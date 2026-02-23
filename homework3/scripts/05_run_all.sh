#!/usr/bin/env bash
set -euo pipefail

bash homework3/scripts/01_download_data.sh
bash homework3/scripts/02_verify_integrity.sh
bash homework3/scripts/03_summarize_genome.sh
bash homework3/scripts/04_summarize_annotation.sh
