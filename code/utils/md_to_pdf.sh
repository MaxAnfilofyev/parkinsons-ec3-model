#!/usr/bin/env bash
# Simple helper script to build the manuscript PDF from Markdown
# Usage (from repo root):
#   bash code/md_to_pdf.sh

set -euo pipefail

md-to-pdf --config-file code/utils/md-to-pdf-config.js Manuscript/pd_manuscript.md
