#!/usr/bin/env python3
"""
Reset failed MP phase DFT relaxation jobs for retry.

Resets RELAX_FAILED -> PENDING so mp_dft_relaxflow.py will resubmit them.

Usage:
    python3 reset_mp_failed.py [--db mp_relax_workflow.json] [--clean] [--dry-run]

Options:
    --db PATH       Path to mp_relax_workflow.json (default: ./mp_relax_workflow.json)
    --clean         Remove entire Relax directory for fresh restart
    --dry-run       Show what would be done without making changes
    --list          List all failed structures and exit

Examples:
    # List failed structures
    python3 reset_mp_failed.py --list

    # Dry run
    python3 reset_mp_failed.py --dry-run --clean

    # Reset and clean directories for fresh restart
    python3 reset_mp_failed.py --clean
"""

import json
import shutil
import argparse
from pathlib import Path


def list_failed(db):
    """List all failed structures."""
    failed = [(mp_id, s) for mp_id, s in db['structures'].items()
              if s['state'] == 'RELAX_FAILED']

    if not failed:
        print("No RELAX_FAILED structures found.")
        return

    print(f"\nRELAX_FAILED: {len(failed)} structures")
    print("-" * 80)
    for mp_id, sdata in failed:
        formula = sdata.get('formula', 'N/A')
        chemsys = sdata.get('chemsys', 'N/A')
        relax_dir = sdata.get('relax_dir', 'N/A')
        print(f"  {mp_id:<15} {formula:<15} {chemsys:<12} {relax_dir}")
    print("-" * 80)


def reset_failed(db_path, clean=False, dry_run=False):
    """Reset RELAX_FAILED -> PENDING."""
    db_path = Path(db_path)

    if not db_path.exists():
        print(f"Error: Database not found at {db_path}")
        return

    print(f"Loading database: {db_path}")
    with open(db_path, 'r') as f:
        db = json.load(f)

    failed = [(mp_id, s) for mp_id, s in db['structures'].items()
              if s['state'] == 'RELAX_FAILED']

    if not failed:
        print("No RELAX_FAILED structures found. Nothing to do.")
        return

    print(f"Found {len(failed)} RELAX_FAILED structures")
    if dry_run:
        print("DRY RUN MODE - no changes will be made\n")
    if not clean:
        print("NOTE: Running without --clean. Use --clean to remove Relax dirs for fresh restart.\n")

    reset_count = 0
    cleaned_dirs = []

    for mp_id, sdata in failed:
        relax_dir = Path(sdata['relax_dir']) if sdata.get('relax_dir') else None

        print(f"  {mp_id} ({sdata.get('formula', '?')}): RELAX_FAILED -> PENDING")
        if relax_dir:
            print(f"    dir: {relax_dir}")

        if not dry_run:
            sdata['state'] = 'PENDING'
            sdata['slurm_id'] = None
            sdata['vasp_energy_per_atom'] = None
            sdata['update_time'] = None
            sdata['submit_time'] = None
            reset_count += 1

        if clean and relax_dir and relax_dir.exists():
            print(f"    Removing Relax directory")
            if not dry_run:
                shutil.rmtree(relax_dir)
                cleaned_dirs.append(str(relax_dir))

    if not dry_run and reset_count > 0:
        backup_path = db_path.with_suffix('.json.bak')
        print(f"\nBacking up database: {backup_path}")
        shutil.copy2(db_path, backup_path)

        print(f"Saving updated database: {db_path}")
        with open(db_path, 'w') as f:
            json.dump(db, f, indent=2)

    print("\n" + "=" * 70)
    print("Reset Summary")
    print("=" * 70)
    print(f"  RELAX_FAILED -> PENDING: {reset_count if not dry_run else len(failed)}")
    if cleaned_dirs:
        print(f"  Cleaned directories: {len(cleaned_dirs)}")
    if dry_run:
        print("\n  DRY RUN - no changes were made")
    else:
        print("\n  Resume mp_dft_relaxflow.py to retry failed jobs")
    print("=" * 70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Reset failed MP phase DFT relaxation jobs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--db', type=str, default='./mp_relax_workflow.json',
        help="Path to mp_relax_workflow.json (default: ./mp_relax_workflow.json)"
    )
    parser.add_argument(
        '--clean', action='store_true',
        help="Remove Relax directories for fresh restart"
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help="Show what would be done without making changes"
    )
    parser.add_argument(
        '--list', action='store_true',
        help="List all failed structures and exit"
    )

    args = parser.parse_args()
    db_path = Path(args.db).expanduser()

    if not db_path.exists():
        print(f"Error: Database not found at {db_path}")
        return 1

    if args.list:
        with open(db_path, 'r') as f:
            db = json.load(f)
        list_failed(db)
        return 0

    reset_failed(db_path, clean=args.clean, dry_run=args.dry_run)
    return 0


if __name__ == '__main__':
    exit(main())
