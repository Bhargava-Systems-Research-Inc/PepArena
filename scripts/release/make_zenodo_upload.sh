#!/usr/bin/env bash
set -euo pipefail
REPO=/12TBDrive1/mega_pep_bench
SRC="${1:-$REPO/dist/PepArena_Proteins_Submission}"
Z="$SRC/zenodo"
OUT="$REPO/dist/zenodo_upload"
TAR=peparena_v1.0_data.tar.gz

[ -d "$Z" ] || { echo "no payload at $Z — run make_submission_bundle.sh first"; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT"

echo "packing $(du -sh "$Z" | cut -f1) from $Z"
tar --format=gnu --sort=name -czf "$OUT/$TAR" -C "$SRC" zenodo

cp "$REPO/scripts/release/templates/zenodo_README.md" "$OUT/README.md"
for c in "$Z"/PepArena_*.csv; do [ -f "$c" ] && cp "$c" "$OUT/"; done
cp "$REPO/scripts/release/templates/ZENODO_UPLOAD_GUIDE.txt" "$OUT/ZENODO_UPLOAD_GUIDE.txt"

echo
echo "dist/zenodo_upload/ ready:"
ls -lh "$OUT"
echo
echo "upload: README.md  PepArena_*.csv (4)  $TAR"
