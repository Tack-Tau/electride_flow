#!/usr/bin/env python
"""
Manual structure relaxation using MatterSim.

This bypasses mattergen-evaluate's BatchRelaxer and uses direct
MatterSimCalculator with ASE optimizers (which we know works).
"""

import json
import zipfile
from pathlib import Path
import torch
from tqdm import tqdm
from ase.io import read, write
from ase.optimize import FIRE
from mattersim.forcefield import MatterSimCalculator
import sys

def relax_structure(atoms, device='cpu', fmax=0.05, max_steps=100):
    """Relax a single structure using MatterSim + FIRE optimizer."""
    
    # Attach calculator
    calc = MatterSimCalculator(
        load_path="MatterSim-v1.0.0-5M.pth",
        device=device
    )
    atoms.calc = calc
    
    # Get initial energy
    energy_initial = atoms.get_potential_energy()
    
    # Relax with FIRE (more robust than BFGS)
    opt = FIRE(atoms, logfile=None)
    opt.run(fmax=fmax, steps=max_steps)
    
    # Get final energy
    energy_final = atoms.get_potential_energy()
    
    return atoms, energy_initial, energy_final


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Manually relax structures with MatterSim")
    parser.add_argument("--structures_zip", type=str, required=True,
                        help="Path to generated_crystals_cif.zip")
    parser.add_argument("--output_dir", type=str, default=None,
                        help="Output directory (default: same as structures_zip)")
    parser.add_argument("--device", type=str, default="cpu",
                        help="Device: cpu or cuda (default: cpu)")
    parser.add_argument("--fmax", type=float, default=0.05,
                        help="Force convergence criterion (eV/Å)")
    parser.add_argument("--max_steps", type=int, default=100,
                        help="Maximum optimization steps")
    parser.add_argument("--max_structures", type=int, default=None,
                        help="Maximum number of structures to process (default: all)")
    
    args = parser.parse_args()
    
    # Determine output directory
    if args.output_dir is None:
        args.output_dir = str(Path(args.structures_zip).parent)
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("Manual Structure Relaxation with MatterSim")
    print("="*70)
    print(f"Input: {args.structures_zip}")
    print(f"Output: {output_dir}")
    print(f"Device: {args.device}")
    print(f"Convergence: fmax={args.fmax} eV/Å, max_steps={args.max_steps}")
    print("="*70)
    
    # Check device
    if args.device == "cuda" and not torch.cuda.is_available():
        print("  CUDA not available, falling back to CPU")
        args.device = "cpu"
    
    # Load structures from zip
    print(f"\nLoading structures from zip...")
    structures = []
    
    with zipfile.ZipFile(args.structures_zip, 'r') as zf:
        cif_files = sorted([f for f in zf.namelist() if f.endswith('.cif')])
        print(f"Found {len(cif_files)} CIF files")
        
        if args.max_structures:
            cif_files = cif_files[:args.max_structures]
            print(f"Processing first {len(cif_files)} structures")
        
        for cif_file in tqdm(cif_files, desc="Loading CIFs"):
            try:
                with zf.open(cif_file) as f:
                    cif_content = f.read().decode('utf-8')
                    from io import StringIO
                    from pymatgen.io.cif import CifParser
                    
                    parser = CifParser(StringIO(cif_content))
                    # primitive=False keeps the original cell (same as MatterGen output)
                    structure = parser.parse_structures(primitive=False)[0]
                    
                    # Convert to ASE atoms
                    from pymatgen.io.ase import AseAtomsAdaptor
                    atoms = AseAtomsAdaptor.get_atoms(structure)
                    
                    structures.append({
                        'cif_file': cif_file,
                        'atoms': atoms,
                        'formula': structure.composition.reduced_formula
                    })
            except Exception as e:
                print(f"  Error loading {cif_file}: {e}")
    
    print(f"Successfully loaded {len(structures)} structures")
    
    # Check for existing progress (checkpointing)
    stability_json = output_dir / "stability.json"
    checkpoint_json = output_dir / "stability_checkpoint.json"
    extxyz_file = output_dir / "relaxed_structures.extxyz"
    
    # Load existing results if available
    results = []
    processed_cifs = set()
    
    if checkpoint_json.exists():
        print(f"\n  Found checkpoint file: {checkpoint_json}")
        with open(checkpoint_json, 'r') as f:
            results = json.load(f)
        processed_cifs = {r['cif_file'] for r in results}
        print(f"  Resuming from checkpoint: {len(results)}/{len(structures)} already processed")
    
    # Relax structures (skip already processed ones)
    print(f"\nRelaxing structures (with checkpointing)...")
    
    for idx, item in enumerate(tqdm(structures, desc="Relaxing")):
        # Skip if already processed
        if item['cif_file'] in processed_cifs:
            continue
        
        try:
            atoms_relaxed, E_initial, E_final = relax_structure(
                item['atoms'].copy(),
                device=args.device,
                fmax=args.fmax,
                max_steps=args.max_steps
            )
            
            result = {
                'cif_file': item['cif_file'],
                'formula': item['formula'],
                'num_atoms': len(atoms_relaxed),
                'energy_initial': float(E_initial),
                'energy_final': float(E_final),
                'energy_per_atom': float(E_final / len(atoms_relaxed)),
                'energy_change': float(E_final - E_initial),
                'converged': True
            }
            results.append(result)
            
            # Append to extxyz file (incremental save)
            write(str(extxyz_file), atoms_relaxed, format='extxyz', append=True)
            
            # Save checkpoint every 10 structures
            if (idx + 1) % 10 == 0:
                with open(checkpoint_json, 'w') as f:
                    json.dump(results, f, indent=2)
            
            # Free memory
            del atoms_relaxed
            if args.device == 'cuda' and torch.cuda.is_available():
                torch.cuda.empty_cache()
            
        except Exception as e:
            print(f"  Error relaxing {item['cif_file']}: {e}")
            results.append({
                'cif_file': item['cif_file'],
                'formula': item['formula'],
                'error': str(e),
                'converged': False
            })
    
    # Save final results
    with open(stability_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Stability results saved to: {stability_json}")
    
    # Remove checkpoint file (no longer needed)
    if checkpoint_json.exists():
        checkpoint_json.unlink()
    
    # Confirm extxyz file was created
    if extxyz_file.exists():
        print(f"  Relaxed structures saved to: {extxyz_file}")
    
    # Print summary
    converged = [r for r in results if r.get('converged', False)]
    
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Total structures: {len(results)}")
    print(f"Successfully relaxed: {len(converged)}")
    
    if converged:
        energies = [r['energy_per_atom'] for r in converged]
        print(f"\nEnergy per atom statistics:")
        print(f"  Mean: {sum(energies)/len(energies):.4f} eV/atom")
        print(f"  Min:  {min(energies):.4f} eV/atom")
        print(f"  Max:  {max(energies):.4f} eV/atom")
        
        energy_changes = [r['energy_change'] for r in converged]
        print(f"\nEnergy change during relaxation:")
        print(f"  Mean: {sum(energy_changes)/len(energy_changes):.4f} eV")
        print(f"  Min:  {min(energy_changes):.4f} eV")
        print(f"  Max:  {max(energy_changes):.4f} eV")
    
    print("="*70)


if __name__ == "__main__":
    main()

