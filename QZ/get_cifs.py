# ase db final.db -c formula,pearson_symbol,e0025,e05,e10,band0,band1,space_group_number,e_above_hull,age -s space_group_number-
from pyxtal.db import database_topology

db = database_topology('final.db')
count = 0
for row in db.db.select():
    spg = row.space_group_number
    if row.space_group_number > 100 and row['e0025'] > 30 and row['e05'] > 30:
        xtal = db.get_pyxtal(row.id)
        name = f'{spg}-{row.formula}-{row.id}.cif'
        xtal.to_file(name, fmt='cif')
        count += 1
        print(name, count, row['e0025'], row['e05'], row['e_above_hull'], db.db.count())
