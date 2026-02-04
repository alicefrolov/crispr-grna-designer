# CRISPR gRNA Designer

Simple tool for designing CRISPR guide RNAs with basic quality checks.

## Features
- PAM sequence detection (NGG)
- GC content filtering
- Off-target risk assessment (simple homopolymer detection)
- Multiple gRNA candidate generation

## Usage
```python
python grna_designer.py "ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
Background
This tool helps identify potential CRISPR-Cas9 guide RNA sequences in a target DNA sequence. It looks for PAM (NGG) sites and evaluates candidates based on:
GC content (optimal: 40-60%)
Homopolymer runs (potential off-target issues)
Guide RNA length (20 bp)
Author
Alice Frolov - Learning CRISPR gene editing principles