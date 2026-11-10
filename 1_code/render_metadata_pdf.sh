#!/usr/bin/env bash
# Render the metadata README to PDF (same folder as the .md).
# Run from anywhere; the script moves to the project root first.
# --to typst avoids needing a LaTeX install (Quarto's default pdf output does).
set -e
cd "$(dirname "$0")/.."
quarto render 3_output/README_metadata.md --to typst
