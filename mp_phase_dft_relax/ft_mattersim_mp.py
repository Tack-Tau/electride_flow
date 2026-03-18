#!/usr/bin/env python3
"""
Fine-tune MatterSim on VASP-relaxed MP phase structures.

Reads VASP relaxation results from mp_relax_workflow.json (or mp_vasp_comparison.json),
extracts structures with DFT energies/forces/stresses from vasprun.xml ionic trajectories,
builds a train/val dataset in extended XYZ format, and launches MatterSim finetuning.

All ionic steps from each VASP relaxation are included, giving the model
diverse geometries with non-trivial forces.

Train/val/test split is done by mp_id (not by frame) to prevent data leakage.
Default split ratio: 80/10/10.
"""

import os
import sys
import json
import argparse
import subprocess
import numpy as np
from pathlib import Path
from datetime import datetime

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write as ase_write
from ase.calculators.singlepoint import SinglePointCalculator


KBAR_TO_EVA3 = 0.1 / 160.21766208


def _atoms_from_ionic_step(structure, energy, forces, stress_kbar, mp_id, step_idx):
    """Convert one ionic step to ASE Atoms with SinglePointCalculator."""
    stress_eV_A3 = np.array(stress_kbar) * KBAR_TO_EVA3
    atoms = AseAtomsAdaptor().get_atoms(structure)
    stress_voigt = stress_eV_A3.flatten()[[0, 4, 8, 5, 2, 1]]
    calc = SinglePointCalculator(
        atoms,
        energy=float(energy),
        forces=np.array(forces, dtype=np.float64),
        stress=stress_voigt,
    )
    atoms.calc = calc
    atoms.info['mp_id'] = mp_id
    atoms.info['ionic_step'] = int(step_idx)
    return atoms


def _extract_from_vasprun(vr, mp_id, skip_first=0, max_force=50.0):
    """
    Extract Atoms from all ionic steps of a parsed Vasprun.

    Skips the first `skip_first` steps, steps where electronic SCF hit NELM,
    and steps with max |force| > max_force.
    """
    entries = []
    nelm = vr.parameters.get('NELM', 120)

    for step_idx, step in enumerate(vr.ionic_steps):
        if step_idx < skip_first:
            continue
        if len(step['electronic_steps']) >= nelm:
            continue
        forces = step['forces']
        if np.max(np.abs(forces)) > max_force:
            continue
        entries.append(_atoms_from_ionic_step(
            step['structure'], step['e_0_energy'],
            forces, step['stress'], mp_id, step_idx,
        ))
    return entries


def collect_structures_from_workflow_db(db_path, skip_first=0, max_force=50.0):
    """
    Collect structures from mp_relax_workflow.json.
    Reads vasprun.xml for each RELAX_DONE/RELAX_TMOUT entry.
    """
    with open(db_path, 'r') as f:
        db = json.load(f)

    entries = []
    skipped = 0

    for mp_id, sdata in db['structures'].items():
        if sdata['state'] not in ('RELAX_DONE', 'RELAX_TMOUT'):
            continue

        relax_dir = sdata.get('relax_dir')
        if not relax_dir:
            skipped += 1
            continue

        vasprun_path = Path(relax_dir) / 'vasprun.xml'
        if not vasprun_path.exists():
            skipped += 1
            continue

        try:
            vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
            frames = _extract_from_vasprun(vr, mp_id, skip_first, max_force)
            if frames:
                entries.extend(frames)
            else:
                skipped += 1

        except Exception as e:
            print(f"  Error processing {mp_id}: {e}")
            skipped += 1

    n_mp = len(set(a.info['mp_id'] for a in entries))
    print(f"Collected {len(entries)} frames from {n_mp} structures ({skipped} skipped)")
    return entries


def collect_structures_from_comparison(comparison_path, workflow_db_path,
                                       skip_first=0, max_force=50.0):
    """
    Collect structures using mp_vasp_comparison.json for validated mp_ids,
    then read vasprun.xml from relax_dir in workflow DB.
    """
    with open(comparison_path, 'r') as f:
        comparison = json.load(f)

    with open(workflow_db_path, 'r') as f:
        db = json.load(f)

    valid_mp_ids = set()
    for entry in comparison.get('structures', []):
        valid_mp_ids.add(entry['mp_id'])

    entries = []
    skipped = 0

    for mp_id in valid_mp_ids:
        sdata = db['structures'].get(mp_id)
        if not sdata:
            skipped += 1
            continue

        relax_dir = sdata.get('relax_dir')
        if not relax_dir:
            skipped += 1
            continue

        vasprun_path = Path(relax_dir) / 'vasprun.xml'
        if not vasprun_path.exists():
            skipped += 1
            continue

        try:
            vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
            frames = _extract_from_vasprun(vr, mp_id, skip_first, max_force)
            if frames:
                entries.extend(frames)
            else:
                skipped += 1

        except Exception as e:
            print(f"  Error processing {mp_id}: {e}")
            skipped += 1

    n_mp = len(set(a.info['mp_id'] for a in entries))
    print(f"Collected {len(entries)} frames from {n_mp} structures ({skipped} skipped)")
    return entries


def train_val_test_split(entries, val_fraction=0.1, test_fraction=0.1, seed=42):
    """
    Split by mp_id (not by frame) to prevent data leakage.
    All frames from the same trajectory go into the same split.
    Default: 80/10/10 train/val/test.
    """
    if val_fraction + test_fraction >= 1.0:
        raise ValueError("val_fraction + test_fraction must be < 1.0")

    mp_id_to_atoms = {}
    for a in entries:
        mp_id_to_atoms.setdefault(a.info['mp_id'], []).append(a)

    mp_ids = sorted(mp_id_to_atoms.keys())
    rng = np.random.default_rng(seed)
    indices = np.arange(len(mp_ids))
    rng.shuffle(indices)
    mp_ids_shuffled = [mp_ids[i] for i in indices]

    n_test = max(1, int(len(mp_ids_shuffled) * test_fraction))
    n_val = max(1, int(len(mp_ids_shuffled) * val_fraction))

    test_ids = mp_ids_shuffled[:n_test]
    val_ids = mp_ids_shuffled[n_test:n_test + n_val]
    train_ids = mp_ids_shuffled[n_test + n_val:]

    train = [a for mid in train_ids for a in mp_id_to_atoms[mid]]
    val = [a for mid in val_ids for a in mp_id_to_atoms[mid]]
    test = [a for mid in test_ids for a in mp_id_to_atoms[mid]]

    return train, val, test


def write_xyz_dataset(atoms_list, output_path):
    """Write list of ASE Atoms to extended XYZ file."""
    ase_write(str(output_path), atoms_list, format='extxyz')
    print(f"  Written {len(atoms_list)} structures to {output_path}")


def run_finetune(train_path, val_path, save_path, model_name, epochs,
                 batch_size, lr, device, include_stresses,
                 patience, re_normalize):
    """Launch MatterSim finetuning via torchrun."""
    cmd = [
        'torchrun', '--nproc_per_node=1',
        '-m', 'mattersim.training.finetune_mattersim',
        '--load_model_path', model_name,
        '--train_data_path', str(train_path),
        '--valid_data_path', str(val_path),
        '--save_path', str(save_path),
        '--epochs', str(epochs),
        '--batch_size', str(batch_size),
        '--lr', str(lr),
        '--device', device,
        '--save_checkpoint',
        '--ckpt_interval', '20',
        '--early_stop_patience', str(patience),
        '--include_forces',
    ]

    if include_stresses:
        cmd.append('--include_stresses')
    if re_normalize:
        cmd.append('--re_normalize')

    print(f"\nFinetuning command:")
    print(f"  {' '.join(cmd)}\n")
    sys.stdout.flush()

    result = subprocess.run(cmd)
    return result.returncode


def main():
    parser = argparse.ArgumentParser(
        description="Fine-tune MatterSim on VASP-relaxed MP phase structures"
    )
    parser.add_argument(
        '--workflow-db',
        type=str,
        default='./mp_relax_workflow.json',
        help="Path to mp_relax_workflow.json"
    )
    parser.add_argument(
        '--comparison-json',
        type=str,
        default=None,
        help="Path to mp_vasp_comparison.json (uses validated structures only)"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='./ft_mattersim',
        help="Output directory for datasets and finetuned model"
    )
    parser.add_argument(
        '--model',
        type=str,
        default='mattersim-v1.0.0-5m',
        choices=['mattersim-v1.0.0-1m', 'mattersim-v1.0.0-5m'],
        help="Pre-trained MatterSim model to finetune"
    )
    parser.add_argument(
        '--epochs',
        type=int,
        default=300,
        help="Number of training epochs"
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=4,
        help="Training batch size (small for <200 structures)"
    )
    parser.add_argument(
        '--lr',
        type=float,
        default=2e-5,
        help="Learning rate (default: 2e-5)"
    )
    parser.add_argument(
        '--val-fraction',
        type=float,
        default=0.1,
        help="Fraction of structures for validation (default: 0.1)"
    )
    parser.add_argument(
        '--test-fraction',
        type=float,
        default=0.1,
        help="Fraction of structures for test set (default: 0.1)"
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help="Random seed for train/val split"
    )
    parser.add_argument(
        '--device',
        type=str,
        default='cuda',
        help="Device for training (cuda or cpu)"
    )
    parser.add_argument(
        '--include-stresses',
        action='store_true',
        default=True,
        help="Include stress tensor in training (default: True)"
    )
    parser.add_argument(
        '--no-stresses',
        action='store_true',
        help="Disable stress tensor in training"
    )
    parser.add_argument(
        '--patience',
        type=int,
        default=50,
        help="Early stopping patience in epochs (default: 50)"
    )
    parser.add_argument(
        '--re-normalize',
        action='store_true',
        help="Re-normalize energy/force scales to the finetuning data"
    )
    parser.add_argument(
        '--dataset-only',
        action='store_true',
        help="Only prepare dataset (train/val/test XYZ files), do not launch training"
    )
    parser.add_argument(
        '--skip-first',
        type=int,
        default=0,
        help="Skip the first N ionic steps per trajectory (default: 0)"
    )
    parser.add_argument(
        '--max-force',
        type=float,
        default=50.0,
        help="Skip ionic steps with max |force| > this value in eV/A (default: 50.0)"
    )

    args = parser.parse_args()

    if args.no_stresses:
        args.include_stresses = False

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("MatterSim Fine-tuning on VASP-relaxed MP Phases")
    print("=" * 70)
    print(f"Output directory: {output_dir}")
    print(f"Model: {args.model}")
    print(f"Epochs: {args.epochs}")
    print(f"Batch size: {args.batch_size}")
    print(f"Learning rate: {args.lr}")
    print(f"Split: train/{args.val_fraction:.0%} val/{args.test_fraction:.0%} test")
    print(f"Include stresses: {args.include_stresses}")
    print(f"Early stop patience: {args.patience}")
    print(f"Re-normalize: {args.re_normalize}")
    print(f"Device: {args.device}")
    print(f"Trajectory filter: skip_first={args.skip_first}, max_force={args.max_force}")
    print("=" * 70 + "\n")

    traj_kw = dict(skip_first=args.skip_first, max_force=args.max_force)

    # Collect structures
    if args.comparison_json and Path(args.comparison_json).exists():
        print(f"Using validated structures from: {args.comparison_json}")
        entries = collect_structures_from_comparison(
            args.comparison_json, args.workflow_db, **traj_kw
        )
    else:
        print(f"Reading structures from workflow DB: {args.workflow_db}")
        entries = collect_structures_from_workflow_db(args.workflow_db, **traj_kw)

    if len(entries) == 0:
        print("ERROR: No valid structures collected. Exiting.")
        return 1

    n_mp = len(set(a.info['mp_id'] for a in entries))
    print(f"\nTotal frames: {len(entries)} from {n_mp} structures")

    # Train/val/test split (by mp_id to prevent data leakage)
    train_set, val_set, test_set = train_val_test_split(
        entries, args.val_fraction, args.test_fraction, args.seed
    )
    train_ids = sorted(set(a.info['mp_id'] for a in train_set))
    val_ids = sorted(set(a.info['mp_id'] for a in val_set))
    test_ids = sorted(set(a.info['mp_id'] for a in test_set))
    print(f"Train: {len(train_set)} frames ({len(train_ids)} structures)")
    print(f"Val:   {len(val_set)} frames ({len(val_ids)} structures)")
    print(f"Test:  {len(test_set)} frames ({len(test_ids)} structures)")

    # Write XYZ datasets
    train_path = output_dir / 'train.xyz'
    val_path = output_dir / 'val.xyz'
    test_path = output_dir / 'test.xyz'
    write_xyz_dataset(train_set, train_path)
    write_xyz_dataset(val_set, val_path)
    write_xyz_dataset(test_set, test_path)

    # Save split metadata
    meta = {
        'timestamp': datetime.now().isoformat(),
        'workflow_db': str(args.workflow_db),
        'comparison_json': str(args.comparison_json) if args.comparison_json else None,
        'skip_first': args.skip_first,
        'max_force': args.max_force,
        'n_structures': n_mp,
        'n_total_frames': len(entries),
        'n_train_frames': len(train_set),
        'n_val_frames': len(val_set),
        'n_test_frames': len(test_set),
        'n_train_structures': len(train_ids),
        'n_val_structures': len(val_ids),
        'n_test_structures': len(test_ids),
        'val_fraction': args.val_fraction,
        'test_fraction': args.test_fraction,
        'seed': args.seed,
        'model': args.model,
        'train_mp_ids': train_ids,
        'val_mp_ids': val_ids,
        'test_mp_ids': test_ids,
    }
    meta_path = output_dir / 'dataset_meta.json'
    with open(meta_path, 'w') as f:
        json.dump(meta, f, indent=2)
    print(f"  Dataset metadata saved to {meta_path}")

    if args.dataset_only:
        print("\n--dataset-only: skipping training.")
        return 0

    # Run finetuning
    print("\n" + "=" * 70)
    print("Starting MatterSim Finetuning")
    print("=" * 70)
    sys.stdout.flush()

    save_path = output_dir / 'model'
    exit_code = run_finetune(
        train_path=train_path,
        val_path=val_path,
        save_path=save_path,
        model_name=args.model,
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        device=args.device,
        include_stresses=args.include_stresses,
        patience=args.patience,
        re_normalize=args.re_normalize,
    )

    if exit_code == 0:
        best_src = save_path / 'best_model.pth'
        if best_src.exists():
            date_str = datetime.now().strftime('%Y%m%d')
            base_name = args.model.replace('.pth', '')
            ft_name = f"{base_name}_{date_str}_ft.pth"
            ft_dst = save_path / ft_name
            best_src.rename(ft_dst)
            print(f"\nRenamed best_model.pth -> {ft_name}")

        print("\n" + "=" * 70)
        print("Finetuning completed successfully!")
        print(f"Model saved to: {save_path}")
        print("=" * 70)
    else:
        print(f"\nERROR: Finetuning failed with exit code {exit_code}")

    return exit_code


if __name__ == '__main__':
    sys.exit(main())
