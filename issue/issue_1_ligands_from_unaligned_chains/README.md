# Issue 1: Ligands from Unaligned Chains

## Problem
When querying with 4HJO (monomer, chain A only), Foldseek returned 5CNN (dimer, chains A and B).
The code was extracting ligands from **all chains** of 5CNN, including chain B which was not aligned.

This caused ANP ligands from chain B to appear "floating" in space, far from the binding site.

## Example
- **Query**: 4HJO (EGFR kinase, chain A only)
- **Hit**: 5CNN (EGFR kinase mutant I682Q, chains A and B)
- **Alignment RMSD**: 0.43 Å (excellent)
- **Problem**: ANP from chain B was not aligned, appeared floating

## Root Cause
`extract_ligands_from_pdb()` function did not filter by chain - it extracted all HETATM records.

## Fix
Added `target_chain` parameter to `extract_ligands_from_pdb()` to only extract ligands from the Foldseek-matched chain.

## Files
- `4HJO.pdb` - Query structure
- `5CNN.pdb` - Hit structure (dimer)
- `check_align.pml` - PyMOL script to reproduce
- `alignment_check_4HJO_5CNN.pse` - PyMOL session showing the issue

## Commit
Fixed in commit 4c678c0
