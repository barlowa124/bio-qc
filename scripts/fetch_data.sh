#!/usr/bin/env bash
# Download the Levine gated FCS benchmark panels from the HDCytoData
# source mirror at UZH. Raw data stays out of git (data/ is gitignored).
set -euo pipefail
cd "$(dirname "$0")/../data"

for panel in Levine_13dim Levine_32dim; do
    zip="${panel}_fcs_files.zip"
    [ -d "$panel" ] || {
        curl -fsSL -O "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/${panel}/${zip}"
        mkdir -p "$panel"
        unzip -q -o "$zip" -d "$panel" -x "__MACOSX/*"
    }
    echo "${panel}: $(ls "$panel"/*.fcs | wc -l) fcs files"
done
