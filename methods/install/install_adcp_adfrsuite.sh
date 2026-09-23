#!/bin/bash
# AutoDock CrankPep (ADCP) via ADFRsuite (Track B). ADFRsuite ships as a self-extracting binary
# tarball from CCSB (not conda). Native cyclic-peptide support (-cyc / -cys). CPU.
# After install, export ADFR_HOME so run_adcp.sh finds the binaries.
set -euo pipefail
DEST=${1:-$HOME/ADFRsuite}
VER=1.0
TARBALL=ADFRsuite_x86_64Linux_${VER}.tar.gz
URL=https://ccsb.scripps.edu/adfr/download/1038/${TARBALL}

mkdir -p "$DEST" && cd "$DEST"
if [ ! -f "$TARBALL" ]; then
  echo "downloading ADFRsuite ${VER} (~500 MB)..."
  wget -q --show-progress "$URL" || curl -L -o "$TARBALL" "$URL"
fi
tar xzf "$TARBALL"
cd ADFRsuite_x86_64Linux_${VER}
echo "Y" | ./install.sh -d "$DEST/install" -c 0    # -c 0 = no color, batch
echo ""
echo "OK. Add to shell / run scripts:"
echo "  export ADFR_HOME=$DEST/install/bin"
echo "  export PATH=\$ADFR_HOME:\$PATH"
echo "Test: \$ADFR_HOME/adcp --help"
# pin the CCSB download id if the mirror URL changes; it is versioned there.
