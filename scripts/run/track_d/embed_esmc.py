from __future__ import annotations
import argparse, hashlib, sys, time
from pathlib import Path
import glob
import numpy as np
import pandas as pd
import torch
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig
from esm.tokenization import get_esmc_model_tokenizers
from esm.utils.constants.esm3 import data_root

ROOT = Path(__file__).resolve().parents[3]
MAN = ROOT / "data/curated/interaction_tpeppro.parquet"
sha1 = lambda s: hashlib.sha1(s.encode()).hexdigest()

ESMC_CFG = {"esmc_300m": (960, 15, 30, "esmc-300"), "esmc_600m": (1152, 18, 36, "esmc-600")}

def load_esmc(name, device):
    d, h, n, key = ESMC_CFG[name]
    model = ESMC(d_model=d, n_heads=h, n_layers=n,
                 tokenizer=get_esmc_model_tokenizers(), use_flash_attn=False).eval()
    pth = glob.glob(str(Path(data_root(key)) / "data/weights/*.pth"))[0]
    state = torch.load(pth, map_location="cpu")
    model.load_state_dict(state.get("model", state) if isinstance(state, dict) else state)
    return model.to(device)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="esmc_300m")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--max-len", type=int, default=2046)
    args = ap.parse_args()
    cache = ROOT / f"runs/track_d/cache/emb_{args.model}"
    cache.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(MAN)
    seqs = sorted(set(df.receptor_seq) | set(df.peptide_seq))
    todo = [(sha1(s), s[: args.max_len]) for s in seqs if not (cache / f"{sha1(s)}.npy").exists()]
    print(f"{len(seqs)} unique seqs; {len(todo)} to embed ({len(seqs)-len(todo)} cached)", flush=True)
    if not todo:
        print("all cached"); return

    client = load_esmc(args.model, args.device)
    cfg = LogitsConfig(sequence=True, return_embeddings=True)
    todo.sort(key=lambda x: len(x[1]))
    t0 = time.time()
    for i, (h, s) in enumerate(todo):
        with torch.no_grad():
            t = client.encode(ESMProtein(sequence=s))
            emb = client.logits(t, cfg).embeddings[0]
        vec = emb[1:-1].mean(0).float().cpu().numpy().astype(np.float32)
        np.save(cache / f"{h}.npy", vec)
        if (i + 1) % 1000 == 0:
            r = (i + 1) / (time.time() - t0)
            print(f"  {i+1}/{len(todo)} ({r:.0f}/s, ETA {(len(todo)-i-1)/r/60:.0f} min)", flush=True)
    print(f"done -> {cache}", flush=True)

if __name__ == "__main__":
    main()
