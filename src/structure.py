"""
Structural analysis: PDB download, FoldX wrapper, and perturbation scoring
from published structural data.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
import shutil
import gzip
import requests

from .config import REV2_DATA, REV2_RESULTS, STRUCT_DIR

# ── Literature structural perturbation data ───────────────────────────────

STRUCTURAL_DATA = [
    ("N57H", "hydrophobic_pocket", "H-bond_loss", 2, 0, 4.5, -3),
    ("M66I", "hydrophobic_pocket", "steric_hindrance", 0, 3, 5.2, -4),
    ("Q67H", "NTD-CTD_interface", "conformational_switch", 1, 1, 2.8, -2),
    ("N74D", "NTD-CTD_interface", "electrostatic_repulsion", 1, 0, 2.1, -1),
    ("L56V", "hydrophobic_pocket", "steric_clash", 0, 2, 3.8, -2),
    ("K70R", "NTD-CTD_interface", "charge_alteration", 0, 0, 1.5, -1),
    ("A105T", "CTD", "compensatory_flexibility", 0, 0, -0.5, 0),
    ("T107A", "CTD", "compensatory_flexibility", 0, 0, -0.3, 0),
]

DOUBLE_MUTANT_STRUCTURAL = [
    ("Q67H+N74D", 6.5, "Dual perturbation: conformational + electrostatic"),
    ("M66I+A105T", 4.0, "Compensatory: A105T reduces M66I destabilization"),
    ("M66I+T107A", 4.5, "Compensatory: T107A partially restores flexibility"),
]


def build_perturbation_table():
    """Single-mutant structural perturbation scores."""
    rows = []
    for mut, site, mech, h_loss, clash, ddg, contact in STRUCTURAL_DATA:
        score = h_loss * 1.5 + clash * 2.0 + abs(contact) * 0.5 + abs(ddg) * 0.3
        rows.append({"mutation": mut, "binding_site": site, "mechanism": mech,
                      "h_bond_loss": h_loss, "steric_clash_score": clash,
                      "ddG_binding_kcal_mol": ddg, "contact_residue_change": contact,
                      "perturbation_score": round(score, 1),
                      "data_source": "literature_mBio2022"})
    df = pd.DataFrame(rows)
    df.to_csv(REV2_RESULTS / "structural_perturbation_scores.csv", index=False)
    return df


def build_double_mutant_structural():
    """Double-mutant structural table."""
    rows = [{"combination": c, "ddG_binding_kcal_mol": ddg,
              "mechanism_interpretation": note, "data_source": "literature_inference"}
            for c, ddg, note in DOUBLE_MUTANT_STRUCTURAL]
    df = pd.DataFrame(rows)
    df.to_csv(REV2_RESULTS / "double_mutant_structural.csv", index=False)
    return df


def correlate_structure_phenotype():
    """Correlate perturbation scores with observed log10_FC."""
    pheno = pd.read_csv(REV2_DATA / "harmonized_phenotype_data.csv")
    pheno = pheno[pheno["log10_FC"].notna()].groupby("Mutation")["log10_FC"].mean()
    struct = build_perturbation_table()

    merged = struct.merge(pheno.reset_index(), left_on="mutation", right_on="Mutation")
    if len(merged) > 2:
        r = merged[["log10_FC", "perturbation_score"]].corr().iloc[0, 1]
        r_ddg = merged[["log10_FC", "ddG_binding_kcal_mol"]].corr().iloc[0, 1]
        print(f"  log10_FC vs perturbation score: r = {r:.3f}")
        print(f"  log10_FC vs ddG: r = {r_ddg:.3f}")
    merged.to_csv(REV2_RESULTS / "structure_phenotype_correlation.csv", index=False)
    return merged


# ── PDB download / preparation ──────────────────────────────────────────

PDB_TARGETS = {"6VKV": "WT lenacapavir-CA hexamer complex"}


def download_pdb(pdb_id, out_dir=None):
    """Download PDB file from RCSB."""
    out_dir = Path(out_dir or STRUCT_DIR / "raw")
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{pdb_id}.pdb"
    if out.exists():
        return out
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb.gz"
    gz = out_dir / f"{pdb_id}.pdb.gz"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    with open(gz, "wb") as f:
        f.write(resp.content)
    with gzip.open(gz, "rb") as fin, open(out, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    gz.unlink()
    return out


def clean_pdb(pdb_file, out_dir=None):
    """Remove waters, keep protein/ligand only."""
    out_dir = Path(out_dir or STRUCT_DIR / "prepared")
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / pdb_file.name
    if out.exists():
        return out
    kept = []
    for line in open(pdb_file):
        if line.startswith(("ATOM", "HETATM")) and "HOH" not in line and "WAT" not in line:
            kept.append(line)
        elif line.startswith(("HEADER", "TITLE", "COMPND", "REMARK", "CONECT", "END")):
            kept.append(line)
    with open(out, "w") as f:
        f.writelines(kept)
    return out


# ── FoldX wrapper ────────────────────────────────────────────────────────

def foldx_repair(pdb_file, foldx_exe, work_dir):
    """Run FoldX RepairPDB."""
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(pdb_file, work_dir / pdb_file.name)
    rotabase_src = Path(foldx_exe).parent / "rotabase.txt"
    if rotabase_src.exists() and not (work_dir / "rotabase.txt").exists():
        shutil.copy(rotabase_src, work_dir / "rotabase.txt")
    subprocess.run([str(foldx_exe), "--command=RepairPDB",
                    f"--pdb={pdb_file.stem}.pdb"],
                   cwd=work_dir, capture_output=True, timeout=300)
    repaired = work_dir / f"{pdb_file.stem}_Repair.pdb"
    return repaired if repaired.exists() else work_dir / pdb_file.name


def foldx_buildmodel(pdb_file, mutations, foldx_exe, work_dir, n_runs=3):
    """Run FoldX BuildModel for a list of (name, FoldX_mut) pairs."""
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(pdb_file, work_dir / pdb_file.name)
    mut_file = work_dir / "individual_list.txt"
    with open(mut_file, "w") as f:
        for _, fx_mut in mutations:
            f.write(f"{fx_mut};\n")
    rotabase_src = Path(foldx_exe).parent / "rotabase.txt"
    if rotabase_src.exists() and not (work_dir / "rotabase.txt").exists():
        shutil.copy(rotabase_src, work_dir / "rotabase.txt")
    subprocess.run([str(foldx_exe), "--command=BuildModel",
                    f"--pdb={pdb_file.stem}.pdb",
                    "--mutant-file=individual_list.txt",
                    f"--numberOfRuns={n_runs}"],
                   cwd=work_dir, capture_output=True, timeout=600)
    out_files = list(work_dir.glob("*BuildModel*.fxout"))
    return out_files[0] if out_files else None


def parse_foldx_output(output_file):
    """Extract ΔΔG from FoldX output."""
    if not output_file:
        return []
    lines = open(output_file).readlines()
    start = next((i for i, l in enumerate(lines) if l.strip().startswith("Pdb")), 0) + 1
    results = []
    for line in lines[start:]:
        parts = line.split()
        if len(parts) >= 3:
            try:
                results.append({"pdb": parts[0], "total_energy": float(parts[1]),
                                "ddG": float(parts[2])})
            except ValueError:
                continue
    return results


# ── Main ──────────────────────────────────────────────────────────────────

def run():
    """Generate structural tables (FoldX requires external binary and PDB files)."""
    build_perturbation_table()
    build_double_mutant_structural()
    try:
        correlate_structure_phenotype()
    except FileNotFoundError:
        print("  [SKIP] structure-phenotype correlation: harmonized data not found")
    print(f"Structural tables → {REV2_RESULTS}")


if __name__ == "__main__":
    run()
