#!/usr/bin/env bash
set -uo pipefail
REPO=/12TBDrive1/mega_pep_bench

PY=""
for c in "${PYTHON:-}" "$HOME/anaconda3/bin/python" "$(command -v python3 || true)"; do
  [ -n "$c" ] && [ -x "$c" ] && "$c" -c "import pandas" 2>/dev/null && { PY="$c"; break; }
done
[ -n "$PY" ] || { echo "ABORT: no python with pandas found. Set PYTHON=/path/to/python."; exit 1; }
MAN=$REPO/manuscript
DRY=0; BULK=1; DOCS=0; DEST=""
for a in "$@"; do
  case "$a" in
    --dry-run) DRY=1 ;;
    --no-bulk) BULK=0 ;;
    --docs-only) DOCS=1 ;;
    *) DEST="$a" ;;
  esac
done
DEST=${DEST:-$REPO/dist/PepArena_Proteins_Submission}
Z="$DEST/zenodo"

say()  { printf '%-58s %s\n' "$1" "$2"; }
size() { du -sh "$1" 2>/dev/null | cut -f1; }
copy() { # copy() SRC DEST_DIR  [as-file]
  local src="$1" dst="$2"
  if [ ! -e "$src" ]; then say "  MISSING $src" "-"; return; fi
  say "  $(basename "$src")" "$(size "$src")"
  [ "$DRY" -eq 1 ] && return
  mkdir -p "$dst"
  cp -r "$src" "$dst"/
  find "$dst/$(basename "$src")" -name __pycache__ -type d -prune -exec rm -rf {} + 2>/dev/null
}

echo "destination: $DEST   (dry-run=$DRY, bulk=$BULK)"
[ "$DRY" -eq 0 ] && mkdir -p "$DEST"

echo; echo "01_Manuscript"
for f in PepArena_Manuscript.docx PepArena_Manuscript_with_figures.docx \
         PepArena_Cover_Letter.docx PepArena_References.docx; do
  copy "$MAN/build/$f" "$DEST/01_Manuscript"
done

echo; echo "02_Figures"
for f in "$MAN"/figures/Figure_*.png "$MAN/figures/Graphical_Abstract.png"; do
  copy "$f" "$DEST/02_Figures"
done

echo; echo "03_Supporting_Information"
copy "$MAN/build/PepArena_Supporting_Information.pdf" "$DEST/03_Supporting_Information"
copy "$MAN/build/PepArena_Supporting_Information.docx" "$DEST/03_Supporting_Information"
for f in "$MAN"/supplementary/Table_S*.csv; do
  copy "$f" "$DEST/03_Supporting_Information/tables"
done

echo; echo "04_Single_File"
copy "$MAN/build/PepArena_Manuscript_single_file.pdf" "$DEST/04_Single_File"
copy "$MAN/build/PepArena_Manuscript_single_file.docx" "$DEST/04_Single_File"

if [ "$DOCS" -eq 1 ]; then
  copy "$REPO/scripts/release/templates/SUBMISSION_CHECKLIST.md" "$DEST"
  echo; echo "zenodo/ left untouched (--docs-only); run without the flag when run artifacts change."
  exit 0
fi

echo; echo "zenodo/dataset  (the citable core: an immutable snapshot of data + splits)"
for f in "$REPO"/data/curated/*.parquet; do copy "$f" "$Z/dataset/curated"; done
for f in "$REPO"/data/splits/*.parquet; do copy "$f" "$Z/dataset/splits"; done
for f in DATASET.md STRATIFICATION.md SOURCES.md LICENSES.md QC_REPORT.md; do
  copy "$REPO/data/$f" "$Z/dataset"
done
copy "$REPO/docs/DATASHEET.md" "$Z/dataset"

echo; echo "zenodo/scores  (four flat tables, one per track)"
if [ "$DRY" -eq 0 ]; then
  rm -rf "$Z/scored"
  "$PY" "$REPO/scripts/release/make_master_table.py" >/dev/null || exit 1
  "$PY" "$REPO/scripts/release/make_track_tables.py" >/dev/null || exit 1
fi
for f in PepArena_scores.csv PepArena_docking.csv PepArena_affinity.csv PepArena_interaction.csv; do
  say "  $f" "$(size "$Z/$f")"
done

echo; echo "zenodo/analysis  (the harness that produced all of it)"
copy "$REPO/scripts" "$Z/analysis"
rm -f "$Z/analysis/scripts/run/track_a/run_cycesmfold2_cyclic.py"
copy "$MAN/src/DATA_FACTS.md" "$Z/analysis"
copy "$MAN/src/verify_facts.py" "$Z/analysis"
copy "$MAN/figures/scripts/make_figures.py" "$Z/analysis"

if [ "$BULK" -eq 1 ]; then
  [ "$DRY" -eq 0 ] && rm -f "$Z/UPLOAD_NOTE.txt"

  echo; echo "zenodo/structures  (native reference complexes)"
  copy "$REPO/data/curated/structures/native" "$Z/structures"

  echo; echo "zenodo/predictions  (Track A best-ranked model per complex -- GPU-weeks to regenerate,"
  echo "                     and AF3's weights are licence-gated, so impossible for some readers)"
  PY_BEST="$REPO/scripts/release/deposit_best_predictions.py"
  if [ "$DRY" -eq 0 ]; then
    "$PY" "$PY_BEST" "$REPO" "$Z/predictions" \
      || { echo "ABORT: predictions step failed; refusing to ship the previous build's structures."; exit 1; }
  else
    "$PY" "$PY_BEST" "$REPO" "$Z/predictions" --dry-run
  fi

  if [ "$DRY" -eq 0 ]; then
    [ -d "$Z/predictions/predictions" ] && mv "$Z/predictions/predictions" "$Z/predictions/canonical"
    [ -d "$Z/predictions/predictions_ncaa" ] && mv "$Z/predictions/predictions_ncaa" "$Z/predictions/ncaa"
  fi

  echo; echo "zenodo/msa  (ColabFold is throttled; these are shared across every model)"
  copy "$REPO/ncaa_ccd/_msa" "$Z/msa"
  copy "$REPO/runs/track_c/boltz2_affinity/msa" "$Z/msa/track_c_affinity"
  if [ "$DRY" -eq 0 ]; then
    [ -d "$Z/msa/_msa" ] && mv "$Z/msa/_msa" "$Z/msa/shared"
    [ -d "$Z/msa/track_c_affinity/msa" ] && mv "$Z/msa/track_c_affinity/msa"/* "$Z/msa/track_c_affinity/" \
      && rmdir "$Z/msa/track_c_affinity/msa"
    find "$Z/msa" -type d -name '_tmp_*' -prune -exec rm -rf {} + 2>/dev/null
  fi

  echo; echo "zenodo/embeddings  (lets the split/distance/negative experiments replay with no GPU)"
  copy "$REPO/runs/track_d/cache" "$Z/embeddings/track_d"
else
  echo; echo "zenodo bulk payloads SKIPPED (--no-bulk)"
  [ "$DRY" -eq 0 ] && cat > "$Z/UPLOAD_NOTE.txt" <<'NOTE'
This copy of zenodo/ is the LIGHT part only: the dataset manifests, splits, provenance,
the four per-track score tables, and the analysis harness.

The four bulk payloads (structures/, predictions/, msa/, embeddings/) are NOT here. Build
a full copy on the machine that holds the runs:

    bash scripts/release/make_submission_bundle.sh /path/to/PepArena_Proteins_Submission

then follow UPLOAD_GUIDE.txt.
NOTE
fi

echo; echo "READMEs (templates -- these were orphaned until now and never reached the deposit)"
copy "$REPO/scripts/release/templates/SUBMISSION_CHECKLIST.md" "$DEST"
[ "$DRY" -eq 0 ] && cp "$REPO/scripts/release/templates/ZENODO_UPLOAD_GUIDE.txt" "$DEST/ZENODO_UPLOAD_GUIDE.txt"
say "  ZENODO_UPLOAD_GUIDE.txt (beside zenodo/, not in it)" "$(size "$REPO/scripts/release/templates/ZENODO_UPLOAD_GUIDE.txt")"

if [ "$DRY" -eq 0 ]; then
  echo; echo "README"
  rm -f "$Z/MANIFEST.sha256"
  echo; echo "TOTAL: $(size "$DEST")   (zenodo/: $(size "$Z"))"
else
  echo; echo "(dry run — nothing copied)"
fi
