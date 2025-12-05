import pandas as pd
from pyxtal.db import database_topology
from pyxtal import pyxtal
from pyxtal.util import new_struc_wo_energy

df = pd.read_csv('electride_analysis.csv')
db = database_topology('electride_data.db')
db1 = database_topology('final.db')
xtals = db1.get_all_xtals()

for row in db.db.select():
    id = row.id
    df_row = df[df['formula'] == row.structure_id]
    #print(df_row)#; import sys; sys.exit()
    if len(df_row) > 0:
        csv_row = df_row.iloc[0]
        if max(csv_row['e0025'], csv_row['e05'], csv_row['e10']) < 20:
            continue
        if max(csv_row['band0'], csv_row['band1']) < 20:
            continue

        xtal = pyxtal()
        for tol in [5e-2, 1e-2, 1e-3, 1e-4]:
            #print(id, tol)
            atoms = row.toatoms()
            xtal.from_seed(atoms, tol=tol)
            if not xtal.valid: continue

            dicts = {}
            for field in ['e0025', 'e05', 'e10', 'band0', 'band1', 'e_above_hull']:
                if field in csv_row.index:
                    dicts[field] = csv_row[field]
            if new_struc_wo_energy(xtal, xtals):
                db1.add_xtal(xtal, dicts)
                xtals.append(xtal)
                print('added:', row.structure_id, xtal.group.number, tol, db1.db.count(), row.id)
            else:
                print('duplicate:', row.structure_id, tol, row.id)
            break
