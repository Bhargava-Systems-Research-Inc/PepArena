from __future__ import annotations
import argparse, hashlib, re, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from transformers import T5Tokenizer, T5EncoderModel

ROOT = Path(__file__).resolve().parents[3]
MAN = ROOT / "data/curated/interaction_tpeppro.parquet"
sha1 = lambda s: hashlib.sha1(s.encode()).hexdigest()
HF = "Rostlab/prot_t5_xl_uniref50"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="prot_t5_xl")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--max-len", type=int, default=1500)
    ap.add_argument("--token-budget", type=int, default=8000)
    args = ap.parse_args()
    cache = ROOT / f"runs/track_d/cache/emb_{args.name}"
    cache.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(MAN)
    seqs = sorted(set(df.receptor_seq) | set(df.peptide_seq))
    todo = [(sha1(s), s[: args.max_len]) for s in seqs if not (cache / f"{sha1(s)}.npy").exists()]
    print(f"{len(seqs)} unique; {len(todo)} to embed", flush=True)
    if not todo:
        print("all cached"); return

    tok = T5Tokenizer.from_pretrained(HF, do_lower_case=False, legacy=True)
    model = T5EncoderModel.from_pretrained(HF).to(args.device).eval()
    if args.device == "cuda":
        model = model.half()
    todo.sort(key=lambda x: len(x[1]))
    t0, i = time.time(), 0
    while i < len(todo):
        batch, toks = [], 0
        while i < len(todo) and (not batch or toks + len(todo[i][1]) <= args.token_budget):
            batch.append(todo[i]); toks += len(todo[i][1]); i += 1
            if len(batch) >= 32:
                break
        proc = [" ".join(re.sub(r"[UZOB]", "X", s)) for _, s in batch]
        enc = tok(proc, add_special_tokens=True, padding="longest", return_tensors="pt")
        ids = enc.input_ids.to(args.device); mask = enc.attention_mask.to(args.device)
        with torch.no_grad():
            rep = model(input_ids=ids, attention_mask=mask).last_hidden_state
        for j, (h, s) in enumerate(batch):
            L = int(mask[j].sum()) - 1
            vec = rep[j, :L].float().mean(0).cpu().numpy().astype(np.float32)
            np.save(cache / f"{h}.npy", vec)
        if i % 2000 < len(batch):
            r = i / (time.time() - t0)
            print(f"  {i}/{len(todo)} ({r:.0f}/s, ETA {(len(todo)-i)/r/60:.0f} min)", flush=True)
    print(f"done -> {cache}", flush=True)

if __name__ == "__main__":
    main()
