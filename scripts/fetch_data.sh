#!/usr/bin/env bash
# Download the Levine_13dim gated FCS files (~8 MB) from the HDCytoData
# source mirror at UZH. Raw data stays out of git (data/ is gitignored).
set -euo pipefail
cd "$(dirname "$0")/../data"
curl -fsSL -O "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/Levine_13dim/Levine_13dim_fcs_files.zip"
mkdir -p Levine_13dim
unzip -q -o Levine_13dim_fcs_files.zip -d Levine_13dim -x "__MACOSX/*"
echo "fetched $(ls Levine_13dim/*.fcs | wc -l) fcs files into data/Levine_13dim/"
