# Pull Request Note

## Summary

This PR adds exact pair-state constraint generation from the S-box lookup table.
It adds support constraints and v2-weighted constraints.
The public API follows the same style as the DDT and LAT APIs.

## Why

These constraints help build exact support and 2-adic MILP models for differential-linear work.
They also make the local S-box model easier to reuse in other projects.

## Code Notes

- The generic pair-state code is in `sboxanalyzer/pair_relations.py`.
- The public methods stay in `sboxanalyzer/__init__.py`.
- The pair-state bases are named `quasidifferential`, `linear`, `subset-difference`, and `ultrametric-difference`.
- The slow PRESENT regression tests are opt-in.