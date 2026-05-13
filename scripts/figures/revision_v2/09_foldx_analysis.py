#!/usr/bin/env python3
"""
FoldX Analysis - Calculate ΔΔG for key resistance mutations
"""

import subprocess
from pathlib import Path
import pandas as pd
import shutil
import os

# Setup paths
BASE_DIR = Path(__file__).parent.parent.parent
FOLDX_DIR = BASE_DIR / "tools" / "foldxWindows"
FOLDX_EXE = FOLDX_DIR / "foldx_20261231.exe"
STRUCT_DIR = BASE_DIR / "data" / "structures" / "prepared"
OUTPUT_DIR = BASE_DIR / "results" / "revision_v2"
WORK_DIR = BASE_DIR / "data" / "structures" / "foldx_work"
WORK_DIR.mkdir(parents=True, exist_ok=True)

# Key mutations to analyze
MUTATIONS = [
    ('N57H', 'NA57H'),  # FoldX format: ChainResidueNewAA
    ('M66I', 'MA66I'),
    ('Q67H', 'MA67H'),
    ('N74D', 'MA74D'),
]

def run_foldx_repair(pdb_file, work_dir):
    """Run FoldX RepairPDB to optimize structure"""

    print(f"  Running RepairPDB on {pdb_file.name}...")

    # Copy PDB to work directory
    work_pdb = work_dir / pdb_file.name
    shutil.copy(pdb_file, work_pdb)

    # Copy rotabase.txt
    rotabase = FOLDX_DIR / "rotabase.txt"
    work_rotabase = work_dir / "rotabase.txt"
    if not work_rotabase.exists():
        shutil.copy(rotabase, work_rotabase)

    # Run RepairPDB
    cmd = [
        str(FOLDX_EXE),
        "--command=RepairPDB",
        f"--pdb={pdb_file.stem}.pdb"
    ]

    try:
        result = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=300
        )

        # Check for repaired file
        repaired_file = work_dir / f"{pdb_file.stem}_Repair.pdb"
        if repaired_file.exists():
            print(f"  [OK] Repaired structure created")
            return repaired_file
        else:
            print(f"  [WARNING] Repair may have failed")
            print(f"  stdout: {result.stdout[:200]}")
            return work_pdb  # Use original if repair failed

    except Exception as e:
        print(f"  [ERROR] RepairPDB failed: {e}")
        return work_pdb

def create_mutation_list(mutations, output_file):
    """Create individual_list.txt for FoldX BuildModel"""

    with open(output_file, 'w') as f:
        for mut_name, foldx_mut in mutations:
            f.write(f"{foldx_mut};\n")

    print(f"  Created mutation list: {len(mutations)} mutations")

def run_foldx_buildmodel(pdb_file, mutation_list, work_dir):
    """Run FoldX BuildModel to calculate ΔΔG"""

    print(f"\n  Running BuildModel...")

    # Copy files to work directory
    work_pdb = work_dir / pdb_file.name
    if not work_pdb.exists():
        shutil.copy(pdb_file, work_pdb)

    work_mutlist = work_dir / "individual_list.txt"
    if work_mutlist != mutation_list:
        shutil.copy(mutation_list, work_mutlist)

    # Ensure rotabase.txt is present
    rotabase = FOLDX_DIR / "rotabase.txt"
    work_rotabase = work_dir / "rotabase.txt"
    if not work_rotabase.exists():
        shutil.copy(rotabase, work_rotabase)

    # Run BuildModel
    cmd = [
        str(FOLDX_EXE),
        "--command=BuildModel",
        f"--pdb={pdb_file.stem}.pdb",
        "--mutant-file=individual_list.txt",
        "--numberOfRuns=3"
    ]

    try:
        result = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=600
        )

        print(f"  BuildModel completed")
        print(f"  Output preview: {result.stdout[:300]}")

        # Look for output files
        output_files = list(work_dir.glob("*BuildModel*.fxout"))
        if output_files:
            print(f"  [OK] Found {len(output_files)} output files")
            return output_files[0]
        else:
            print(f"  [WARNING] No output files found")
            return None

    except Exception as e:
        print(f"  [ERROR] BuildModel failed: {e}")
        return None

def parse_foldx_output(output_file):
    """Parse FoldX output file to extract ΔΔG values"""

    if not output_file or not output_file.exists():
        return []

    results = []

    with open(output_file, 'r') as f:
        lines = f.readlines()

    # Find data section (after header)
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip().startswith('Pdb'):
            data_start = i + 1
            break

    # Parse data lines
    for line in lines[data_start:]:
        if line.strip() and not line.startswith('#'):
            parts = line.split()
            if len(parts) >= 3:
                try:
                    pdb_name = parts[0]
                    total_energy = float(parts[1])
                    ddg = float(parts[2]) if len(parts) > 2 else 0.0

                    results.append({
                        'pdb': pdb_name,
                        'total_energy': total_energy,
                        'ddG': ddg
                    })
                except:
                    continue

    return results

def main():
    """Main execution"""
    print("="*60)
    print("FoldX ΔΔG Analysis")
    print("="*60)

    # Check FoldX
    if not FOLDX_EXE.exists():
        print(f"ERROR: FoldX not found at {FOLDX_EXE}")
        return

    print(f"FoldX executable: {FOLDX_EXE}")
    print(f"Work directory: {WORK_DIR}")

    # Find prepared structures
    pdb_files = list(STRUCT_DIR.glob("*.pdb"))
    if not pdb_files:
        print("ERROR: No prepared PDB files found")
        return

    print(f"\nFound {len(pdb_files)} PDB files:")
    for pdb in pdb_files:
        print(f"  {pdb.name}")

    # Use 6VKV (WT structure)
    wt_pdb = STRUCT_DIR / "6VKV.pdb"
    if not wt_pdb.exists():
        print("ERROR: 6VKV.pdb not found")
        return

    # Repair structure
    print("\n" + "="*60)
    print("Step 1: Repair Structure")
    print("="*60)
    repaired_pdb = run_foldx_repair(wt_pdb, WORK_DIR)

    # Create mutation list
    print("\n" + "="*60)
    print("Step 2: Create Mutation List")
    print("="*60)
    mutation_list = WORK_DIR / "individual_list.txt"
    create_mutation_list(MUTATIONS, mutation_list)

    # Run BuildModel
    print("\n" + "="*60)
    print("Step 3: Run BuildModel")
    print("="*60)
    output_file = run_foldx_buildmodel(repaired_pdb, mutation_list, WORK_DIR)

    # Parse results
    print("\n" + "="*60)
    print("Step 4: Parse Results")
    print("="*60)

    if output_file:
        results = parse_foldx_output(output_file)

        if results:
            results_df = pd.DataFrame(results)
            results_df.to_csv(OUTPUT_DIR / "foldx_results.csv", index=False)
            print(f"[OK] Saved results: {len(results)} mutations")
            print("\nΔΔG values:")
            for _, row in results_df.iterrows():
                print(f"  {row['pdb']}: ΔΔG = {row['ddG']:.2f} kcal/mol")
        else:
            print("[WARNING] No results parsed")
    else:
        print("[ERROR] No output file to parse")

    print("\n" + "="*60)
    print("FoldX analysis complete!")
    print("="*60)

if __name__ == "__main__":
    main()
