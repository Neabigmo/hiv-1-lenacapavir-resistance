#!/usr/bin/env python3
"""
Structure Preparation - Download and prepare PDB structures for FoldX analysis
"""

import requests
from pathlib import Path
import gzip
import shutil

# Setup paths
BASE_DIR = Path(__file__).parent.parent.parent
STRUCT_DIR = BASE_DIR / "data" / "structures"
STRUCT_DIR.mkdir(parents=True, exist_ok=True)

RAW_DIR = STRUCT_DIR / "raw"
PREPARED_DIR = STRUCT_DIR / "prepared"
RAW_DIR.mkdir(exist_ok=True)
PREPARED_DIR.mkdir(exist_ok=True)

# PDB IDs to download
PDB_IDS = {
    '6VKV': 'WT lenacapavir-CA hexamer complex',
    '7RAO': 'M66I mutant capsid',
}

def download_pdb(pdb_id):
    """Download PDB structure from RCSB"""

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb.gz"
    output_gz = RAW_DIR / f"{pdb_id}.pdb.gz"
    output_pdb = RAW_DIR / f"{pdb_id}.pdb"

    if output_pdb.exists():
        print(f"  {pdb_id}.pdb already exists, skipping download")
        return output_pdb

    print(f"  Downloading {pdb_id} from RCSB...")

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        # Save compressed file
        with open(output_gz, 'wb') as f:
            f.write(response.content)

        # Decompress
        with gzip.open(output_gz, 'rb') as f_in:
            with open(output_pdb, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        # Remove compressed file
        output_gz.unlink()

        print(f"  [OK] Downloaded {pdb_id}.pdb")
        return output_pdb

    except Exception as e:
        print(f"  [ERROR] Failed to download {pdb_id}: {e}")
        return None

def clean_pdb_for_foldx(pdb_file):
    """Clean PDB file for FoldX (remove waters, keep only protein and ligand)"""

    output_file = PREPARED_DIR / pdb_file.name

    if output_file.exists():
        print(f"  {output_file.name} already prepared, skipping")
        return output_file

    print(f"  Cleaning {pdb_file.name} for FoldX...")

    kept_lines = []

    with open(pdb_file, 'r') as f:
        for line in f:
            # Keep ATOM and HETATM lines (protein and ligand)
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Skip water molecules
                if 'HOH' not in line and 'WAT' not in line:
                    kept_lines.append(line)
            # Keep essential header lines
            elif line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'REMARK')):
                kept_lines.append(line)
            # Keep connectivity
            elif line.startswith(('CONECT', 'END')):
                kept_lines.append(line)

    # Write cleaned file
    with open(output_file, 'w') as f:
        f.writelines(kept_lines)

    print(f"  [OK] Cleaned {output_file.name} ({len(kept_lines)} lines)")
    return output_file

def main():
    """Main execution"""
    print("="*60)
    print("Structure Preparation for FoldX")
    print("="*60)

    # Download structures
    print("\nDownloading PDB structures...")
    downloaded = {}

    for pdb_id, description in PDB_IDS.items():
        print(f"\n{pdb_id}: {description}")
        pdb_file = download_pdb(pdb_id)
        if pdb_file:
            downloaded[pdb_id] = pdb_file

    # Clean structures
    print("\n" + "="*60)
    print("Cleaning structures for FoldX...")
    print("="*60)

    prepared = {}
    for pdb_id, pdb_file in downloaded.items():
        print(f"\n{pdb_id}:")
        cleaned_file = clean_pdb_for_foldx(pdb_file)
        if cleaned_file:
            prepared[pdb_id] = cleaned_file

    # Summary
    print("\n" + "="*60)
    print("Structure Preparation Complete")
    print("="*60)
    print(f"\nPrepared structures in: {PREPARED_DIR}")
    for pdb_id, pdb_file in prepared.items():
        print(f"  {pdb_id}: {pdb_file.name}")

    print("\nReady for FoldX analysis!")
    print("="*60)

if __name__ == "__main__":
    main()
