#!/usr/bin/env bash
# Installs the APA Python packages from GitHub into the trident environment.
# Run by create_environment.sh.

set -euo pipefail
pip install --no-deps \
    "git+https://github.com/idiap/apatools.git@359bda1995d68f71f5c9c4647eae963dfc0badf8" \
    "git+https://github.com/idiap/bulkanalysis.git@4d6023eed9958a289a851120b671f9b7899a93c5" \
    "git+https://github.com/idiap/genomemanager.git@a911a8fb933600af6c8e06179d78a3b10b89148e" \
    "git+https://github.com/idiap/kallistomanager.git@c91b1d41aa0306ae4812d79ab4236984e98f52ea" \
    "statannotations==0.5.0"
