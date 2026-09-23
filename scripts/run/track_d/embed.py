from __future__ import annotations
import argparse, hashlib, sys, time
from pathlib import Path
import urllib.request
import numpy as np
import pandas as pd
import torch
import esm

ROOT = Path(__file__).resolve().parents[3]
MAN = ROOT / "data/curated/interaction_tpeppro.parquet"

def sha1(s):
    return hashlib.sha1(s.encode()).hexdigest()

def load_model(name, device, fp16):
    if not fp16:
        model, alphabet = getattr(esm.pretrained, name)()
        return model.eval().to(device), alphabet
    hub = Path(torch.hub.get_dir()) / "checkpoints"
    hub.mkdir(parents=True, exist_ok=True)
    ckpt = hub / f"{name}.pt"
    if not ckpt.exists():
        url = f"https://dl.fbaipublicfiles.com/fair-esm/models/{name}.pt"
        print(f"  downloading {url} (large) ...", flush=True)
        urllib.request.urlretrieve(url, ckpt)
    md = torch.load(ckpt, map_location="cpu", mmap=True, weights_only=False)
    with torch.device(device):
        model, alphabet = esm.pretrained.load_model_and_alphabet_core(name, md, None)
    return model.half().eval(), alphabet

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="esm2_t12_35M_UR50D")
    ap.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"])
    ap.add_argument("--fp16", action="store_true", help="fp16 + mmap + on-GPU construct (15B)")
    ap.add_argument("--threads", type=int, default=24)
    ap.add_argument("--max-len", type=int, default=1022)
    ap.add_argument("--token-budget", type=int, default=12000)
    args = ap.parse_args()
    torch.set_num_threads(args.threads)
    device = ("cuda" if torch.cuda.is_available() else "cpu") if args.device == "auto" else args.device
    print(f"device: {device}", flush=True)

    cache = ROOT / f"runs/track_d/cache/emb_{args.model}"
    cache.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(MAN)
    seqs = sorted(set(df.receptor_seq) | set(df.peptide_seq))
    todo = [(sha1(s), s[: args.max_len]) for s in seqs if not (cache / f"{sha1(s)}.npy").exists()]
    print(f"{len(seqs)} unique seqs; {len(todo)} to embed ({len(seqs)-len(todo)} cached). model={args.model}", flush=True)
    if not todo:
        print("all cached — nothing to do"); return

    model, alphabet = load_model(args.model, device, args.fp16)
    bc = alphabet.get_batch_converter()
    repr_layer = model.num_layers

    todo.sort(key=lambda x: len(x[1]))
    done, t0, i = 0, time.time(), 0
    while i < len(todo):
        batch, toks = [], 0
        while i < len(todo) and (not batch or toks + len(todo[i][1]) + 2 <= args.token_budget):
            batch.append(todo[i]); toks += len(todo[i][1]) + 2; i += 1
            if len(batch) >= 64:
                break
        _, _, tokens = bc([(h, s) for h, s in batch])
        tokens = tokens.to(device)
        with torch.no_grad():
            out = model(tokens, repr_layers=[repr_layer])["representations"][repr_layer]
        for j, (h, s) in enumerate(batch):
            L = len(s)
            vec = out[j, 1:L + 1].mean(0).float().cpu().numpy().astype(np.float32)
            np.save(cache / f"{h}.npy", vec)
        done += len(batch)
        if done % 512 < len(batch):
            rate = done / (time.time() - t0 + 1e-9)
            print(f"  embedded {done}/{len(todo)}  ({rate:.1f}/s, ETA {(len(todo)-done)/rate/60:.0f} min)", flush=True)
    print(f"done: embedded {done} sequences -> {cache}", flush=True)

if __name__ == "__main__":
    main()
