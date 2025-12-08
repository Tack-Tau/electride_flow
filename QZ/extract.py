import pandas as pd
from pyxtal.db import database_topology
from pyxtal import pyxtal
from pyxtal.util import new_struc_wo_energy

df = pd.read_csv('electride_analysis.csv')
db = database_topology('electride_data_18k.db')
db1 = database_topology('final.db')
xtals = db1.get_all_xtals()

for row in db.db.select():
    id = row.id
    if id == 2000: break
    df_row = df[df['formula'] == row.structure_id]
    if len(df_row) > 0:
        csv_row = df_row.iloc[0]
        if max(csv_row['e0025'], csv_row['e05'], csv_row['e10']) < 20:
            continue
        if max(csv_row['band0'], csv_row['band1']) < 20:
            continue

        add = False
        xtal = pyxtal()
        for tol in [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]:
            atoms = row.toatoms()#; atoms.write('tmp.vasp', vasp5=True, format='vasp')
            xtal.from_seed(atoms, tol=tol)
            add = True
            if not xtal.valid or xtal.source is None: continue
            #print(xtal, id, tol)
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

        if not add:
            atoms.write('tmp.vasp', vasp5=True, direct=True, format='vasp')
            print('Error in spg determiation, check tmp.vasp', tol); import sys; sys.exit()
