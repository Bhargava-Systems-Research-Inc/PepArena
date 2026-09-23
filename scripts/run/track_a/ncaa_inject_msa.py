from __future__ import annotations
import json, re
from pathlib import Path

ROOT = Path("/12TBDrive1/mega_pep_bench")
NC = ROOT / "ncaa_ccd"
MSA = NC / "_msa"
PROT_MSA_DIRS = ROOT / "runs/track_a/predictions_ncaa/_protenix_msa"

def a3m_text(seq_id):
    p = MSA / f"{seq_id}.a3m"
    return p.read_text() if p.exists() else None

def main():
    cmap = json.load(open(MSA / "receptor_map.json"))
    for m in ("protenix_msa", "boltz2_msa", "af3_msa"):
        (NC / m).mkdir(parents=True, exist_ok=True)
    n = ok = 0
    for cid, rec_ids in cmap.items():
        n += 1
        pj = NC / "protenix" / f"{cid}.json"
        if pj.exists():
            d = json.loads(pj.read_text())
            chains = d[0]["sequences"]
            for idx, ch in enumerate(chains):
                dest = PROT_MSA_DIRS / cid / str(idx)
                dest.mkdir(parents=True, exist_ok=True)
                if idx < len(rec_ids) and a3m_text(rec_ids[idx]):
                    (dest / "non_pairing.a3m").write_text(a3m_text(rec_ids[idx]))
                else:
                    (dest / "non_pairing.a3m").write_text(f">query\n{ch['proteinChain']['sequence']}\n")
                ch["proteinChain"]["msa"] = {"precomputed_msa_dir": str(dest.resolve()),
                                             "pairing_db": "uniref100"}
            (NC / "protenix_msa" / f"{cid}.json").write_text(json.dumps(d, indent=2))
        by = NC / "boltz2" / f"{cid}.yaml"
        if by.exists():
            lines, ri = by.read_text().splitlines(), 0
            out = []
            for ln in lines:
                out.append(ln)
                mobj = re.match(r"^      sequence: (.+)$", ln)
                if mobj:
                    if ri < len(rec_ids) and (MSA / f"{rec_ids[ri]}.a3m").exists():
                        out.append(f"      msa: {(MSA / (rec_ids[ri] + '.a3m')).resolve()}")
                    else:
                        out.append("      msa: empty")
                    ri += 1
            (NC / "boltz2_msa" / f"{cid}.yaml").write_text("\n".join(out) + "\n")
        aj = NC / "af3" / f"{cid}.json"
        if aj.exists():
            d = json.loads(aj.read_text())
            ci = 0
            for s in d["sequences"]:
                if "protein" not in s:
                    continue
                seq = s["protein"]["sequence"]
                a3m = a3m_text(rec_ids[ci]) if ci < len(rec_ids) else None
                if s["protein"].get("modifications"):
                    s["protein"]["unpairedMsa"] = ""
                else:
                    s["protein"]["unpairedMsa"] = a3m if a3m else f">query\n{seq}\n"
                s["protein"]["pairedMsa"] = ""
                s["protein"]["templates"] = []
                ci += 1
            (NC / "af3_msa" / f"{cid}.json").write_text(json.dumps(d, indent=2))
        ok += 1
    print(f"injected MSA into {ok}/{n} complexes x 3 models")
    print(f"  protenix_msa: {len(list((NC/'protenix_msa').glob('*.json')))} | "
          f"boltz2_msa: {len(list((NC/'boltz2_msa').glob('*.yaml')))} | "
          f"af3_msa: {len(list((NC/'af3_msa').glob('*.json')))}")

if __name__ == "__main__":
    main()
