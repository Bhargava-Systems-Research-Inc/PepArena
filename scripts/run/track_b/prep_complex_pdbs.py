#!/usr/bin/env python3
import os
import pandas as pd
import gemmi

REPO = "/12TBDrive1/mega_pep_bench"
FRESH = os.environ.get("FRESH", f"{REPO}/runs/track_b/fresh")

def write_chain_pdb(structure, chains, out_pdb):
    sel = gemmi.Structure()
    sel.spacegroup_hm = structure.spacegroup_hm
    model = gemmi.Model("1")
    for ch in structure[0]:
        if ch.name in chains:
            model.add_chain(ch.clone())
    sel.add_model(model)
    sel.setup_entities()
    sel.write_pdb(out_pdb)

def main():
    wl = pd.read_csv(f"{FRESH}/worklist.csv")
    ok = 0
    for r in wl.itertuples(index=False):
        outdir = f"{FRESH}/inputs/{r.complex_id}"
        os.makedirs(outdir, exist_ok=True)
        try:
            st = gemmi.read_structure(r.native_path)
            rec_chains = set(str(r.receptor_chains).split(","))
            pep_chains = {str(r.peptide_chain)}
            write_chain_pdb(st, rec_chains, f"{outdir}/receptor.pdb")
            write_chain_pdb(st, pep_chains, f"{outdir}/peptide.pdb")
            with open(f"{outdir}/peptide.fasta", "w") as fh:
                fh.write(f">{r.complex_id}_peptide\n{r.peptide_seq}\n")
            ok += 1
        except Exception as e:
            print("FAIL", r.complex_id, e)
    print(f"materialized {ok}/{len(wl)} receptor/peptide pairs -> {FRESH}/inputs/")

if __name__ == "__main__":
    main()
