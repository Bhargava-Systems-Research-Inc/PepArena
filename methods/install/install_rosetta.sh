#!/bin/bash
# Rosetta FlexPepDock / PIPER-FlexPepDock (Track B). LICENSE-GATED: get a free academic license
# at https://www.rosettacommons.org/software/license-and-download, then this builds the binaries.
# Multi-GB source + long compile. Provided as a documented script, NOT run automatically.
set -euo pipefail
SRC=${ROSETTA_SRC:-$HOME/rosetta}      # point at the extracted Rosetta source you downloaded
if [ ! -d "$SRC/source" ]; then
  cat <<MSG
Rosetta source not found at $SRC.
1. Register + download the Rosetta source bundle (academic license, free):
     https://www.rosettacommons.org/software/license-and-download
2. Extract it, then: ROSETTA_SRC=/path/to/rosetta.source.<ver> bash $0
MSG
  exit 1
fi
cd "$SRC/source"
./scons.py -j"$(nproc)" mode=release bin/FlexPepDocking.linuxgccrelease
echo "OK: FlexPepDock built at $SRC/source/bin/"
echo "For PIPER-FlexPepDock also fetch PIPER (separate academic license)."
