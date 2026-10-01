#!/usr/bin/env bash
# Install the Oxford Nanopore basecaller used for the re-basecalling experiment (task 5).
# Needs an NVIDIA GPU for useful speed (results were produced on an RTX 4060 Ti, 8 GB, via WSL2).
# Models (hac, sup for R10.4.1 5 kHz) are downloaded by dorado itself on first use.
set -euo pipefail
VER=1.4.0
mkdir -p ~/tools && cd ~/tools
[ -d dorado-$VER-linux-x64 ] || curl -sSL https://cdn.oxfordnanoportal.com/software/analysis/dorado-$VER-linux-x64.tar.gz | tar xz
~/tools/dorado-$VER-linux-x64/bin/dorado --version
