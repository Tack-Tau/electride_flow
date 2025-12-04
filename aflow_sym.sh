#!/bin/bash
# Analyze symmetry of electride CIF files using AFLOW with multi-threading

set -e

CIF_DIR="electride_CIF"
ELECTRIDE_CSV="Bin-Ele-HT/electride_analysis.csv"
ELECTRIDE_DB="Bin-Ele-HT/electride_data.db"
OUTPUT_CSV="electride_candidates.csv"
TEMP_DIR="aflow_temp"

# Check if input files exist
if [ ! -d "$CIF_DIR" ]; then
    echo "ERROR: CIF directory not found: $CIF_DIR"
    exit 1
fi

if [ ! -f "$ELECTRIDE_CSV" ]; then
    echo "ERROR: Electride analysis CSV not found: $ELECTRIDE_CSV"
    exit 1
fi

if [ ! -f "$ELECTRIDE_DB" ]; then
    echo "ERROR: Electride database not found: $ELECTRIDE_DB"
    exit 1
fi

# Detect number of cores for multi-threading
if [[ "$OSTYPE" == "darwin"* ]]; then
    NCORES=$(sysctl -n hw.ncpu)
else
    NCORES=$(nproc)
fi

echo "==================================================================="
echo "AFLOW Symmetry Analysis for Electride Candidates"
echo "==================================================================="
echo "CIF directory: $CIF_DIR"
echo "Electride CSV: $ELECTRIDE_CSV"
echo "Electride DB: $ELECTRIDE_DB"
echo "Using $NCORES CPU cores"
echo "==================================================================="

# Query database for true electrides only (before changing directory)
# Extract structure_id from key_value_pairs JSON field
echo "Querying database for structures with is_electride=true..."
electride_list=$(sqlite3 "$ELECTRIDE_DB" "SELECT key_value_pairs FROM systems WHERE key_value_pairs LIKE '%is_electride\": true%' ORDER BY id;" | \
    python3 -c "import sys, json; [print(json.loads(line)['structure_id']) for line in sys.stdin if line.strip()]")
total=$(echo "$electride_list" | wc -l | xargs)

echo "Found $total true electride candidates"
echo "==================================================================="

# Create temporary directory for aflow intermediate files
mkdir -p "$TEMP_DIR"
cd "$TEMP_DIR"

# Create temporary CSV for aflow results
echo "structure_id,aflow_spacegroup" > aflow_results.csv

# Process each true electride
count=0
echo "Analyzing symmetry with AFLOW..."

while IFS= read -r structure_id; do
    # Skip empty lines
    if [ -z "$structure_id" ]; then
        continue
    fi
    
    cif_file="../$CIF_DIR/${structure_id}.cif"
    
    if [ ! -f "$cif_file" ]; then
        echo "WARNING: CIF not found for $structure_id, using original spacegroup"
        echo "$structure_id,NA" >> aflow_results.csv
    else
        # Run aflow with multi-threading and extract space group number
        # Redirect stderr to suppress verbose output and run in temp dir to contain intermediate files
        sg_output=$(aflow --np=$NCORES --aflowSG < "$cif_file" 2>/dev/null || echo "ERROR")
        
        if [ "$sg_output" == "ERROR" ]; then
            echo "WARNING: aflow failed for $structure_id, using original spacegroup"
            sg_num="NA"
        else
            # Extract space group number from output like "C2/m #12"
            sg_num=$(echo "$sg_output" | grep -o '#[0-9]*' | sed 's/#//')
            
            # If no number found, try parsing differently
            if [ -z "$sg_num" ]; then
                sg_num="NA"
            fi
        fi
        
        echo "$structure_id,$sg_num" >> aflow_results.csv
    fi
    
    count=$((count + 1))
    if [ $((count % 100)) -eq 0 ]; then
        echo "  Processed $count/$total structures..."
    fi
done <<< "$electride_list"

echo "  Processed $count/$total structures (complete)"

# Clean up aflow intermediate files
rm -f aflow.*.out.xz 2>/dev/null || true

cd ..

# Now merge the results: use Python for reliable CSV handling
python3 << 'PYTHON_EOF'
import csv
import sys
import sqlite3
import json

# Get list of true electride structure_ids from database
conn = sqlite3.connect('Bin-Ele-HT/electride_data.db')
cursor = conn.cursor()
cursor.execute("SELECT key_value_pairs FROM systems WHERE key_value_pairs LIKE '%is_electride\": true%'")
rows = cursor.fetchall()
conn.close()

# Extract structure_id from JSON
true_electride_ids = set()
for row in rows:
    try:
        kvp = json.loads(row[0])
        structure_id = kvp.get('structure_id')
        if structure_id:
            true_electride_ids.add(structure_id)
    except:
        pass

print(f"Found {len(true_electride_ids)} true electrides in database")

# Read original electride analysis data (all structures)
all_data = {}
with open('Bin-Ele-HT/electride_analysis.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        all_data[row['formula']] = row

# Filter to only true electrides
electride_data = []
for struct_id in true_electride_ids:
    if struct_id in all_data:
        electride_data.append(all_data[struct_id])
    else:
        print(f"WARNING: {struct_id} in database but not in CSV")

print(f"Matched {len(electride_data)} structures from CSV")

# Read aflow results
aflow_results = {}
with open('aflow_temp/aflow_results.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        aflow_results[row['structure_id']] = row['aflow_spacegroup']

# Merge: update spacegroup with aflow results
updated_count = 0
for row in electride_data:
    formula = row['formula']
    if formula in aflow_results:
        aflow_sg = aflow_results[formula]
        if aflow_sg != 'NA':
            row['spacegroup'] = aflow_sg
            updated_count += 1

# Write updated CSV (sorted by e_above_hull as in original)
electride_data.sort(key=lambda x: float(x['e_above_hull']))

with open('electride_candidates.csv', 'w', newline='') as f:
    fieldnames = ['formula', 'composition', 'e0025', 'e05', 'e10', 'band0', 'band1', 'spacegroup', 'e_above_hull']
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(electride_data)

print(f"Updated spacegroup for {updated_count}/{len(electride_data)} structures")

PYTHON_EOF

# Clean up temporary directory
rm -rf "$TEMP_DIR"

echo "==================================================================="
echo "Analysis Complete"
echo "==================================================================="
echo "Output file: $OUTPUT_CSV"
echo ""
echo "Top 10 space groups in electride candidates:"
tail -n +2 "$OUTPUT_CSV" | awk -F',' '{count[$8]++} END {for (sg in count) print count[sg], sg}' | sort -rn | head -10

echo ""
echo "High symmetry electride candidates (cubic, space group >= 195):"
cubic_count=$(awk -F',' '$8 >= 195 {count++} END {print count+0}' "$OUTPUT_CSV")
echo "  Found $cubic_count cubic structures"

if [ $cubic_count -gt 0 ]; then
    echo ""
    echo "  Cubic electrides:"
    awk -F',' 'NR==1 {print $1","$2","$8","$9; next} $8 >= 195 {print $1","$2","$8","$9}' "$OUTPUT_CSV" | column -t -s','
fi

echo "==================================================================="
