import os
from pyxtal.db import database_topology
from pyxtal import pyxtal

if not os.path.exists('cifs'): os.mkdir('cifs')
if os.path.exists('stable_electrides.db'): os.remove('stable_electrides.db')

db1 = database_topology('Bin-Ele-HT/stable_electrides.db')
db2 = database_topology('Ter-Ele-HT/stable_electrides.db')
db = database_topology('stable_electrides.db')
count = 0
for db0 in [db1, db2]:
    for row in db0.db.select():
        id = row.id
        xtal = pyxtal()
        atoms = row.toatoms()
        xtal.from_seed(atoms, tol=1e-2)

        dicts = {"e0025": row.e0025,
                 "e05": row.e05,
                 "e10": row.e10,
                 "band0": row.band0,
                 "band1": row.band1,
                 "vasp_e_ref_hull": row.dft_e_hull,
                 "vasp_energy": row.vasp_energy_per_atom,
                 "ff_energy": row.mattersim_energy_per_atom,
                 "ff_e_ref_hull": row.mattersim_e_hull,
                 "structure_id": row.structure_id,
                 "N_elements": len(xtal.numIons)}
        status, dup_id = db.check_new_structure(xtal, same_group=True, d_tol=0.05, return_id=True)

        if status:
            count += 1
            db.add_xtal(xtal, dicts)
            xtal.to_file(f"cifs/{row.composition}-spg{row.space_group_number}-ehull-{row.dft_e_hull:.2f}.cif")#, format='cif')
        else:
            row0 = db.get_row(dup_id)
            if row.dft_e_hull < row0.vasp_e_ref_hull:
                db.db.delete([dup_id])
                db.add_xtal(xtal, dicts)
                print("Update", id, xtal.formula)
            print("Duplicate", id, xtal.formula)
    print(db0.db_name, count)

# To do: add new tags e_hull_refined, eletride, N_excess, band_gap, .etc
for row in db0.db.select():
    db.db.update(row.id, e_hull_refined=0.0, electride='Yes')
