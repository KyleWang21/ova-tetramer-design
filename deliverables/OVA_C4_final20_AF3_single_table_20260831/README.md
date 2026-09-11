# OVA C4 final 20 — single-table package

- `all20_unified_summary.tsv`: the only table in this package. It has one row per candidate and joins sequence/mutations, AF3 aggregate statistics, five seed summaries, all 25 AF3 ipTM/weak-interface/interface-PAE values, mutation-interface recurrence, uniform Protenix-v2 seed31 results, and summarized ProLIF chain-pair/hotspot information.
- `cif/`: exactly one AF3 CIF per candidate. Filenames start with the zero-padded AF3 rank (`01_`–`20_`); each is the highest-ipTM AF3 model among the fixed 5 seeds × 5 samples for that sequence.
- `all20_sequences.fasta`: all 20 subunit sequences.

Column prefixes identify their source: `af3_`, `protenix_`, and `prolif_`. The `af3_cif` field gives the package-relative coordinate path. ProLIF one-to-many details are compacted into semicolon-separated hotspot fields so the delivery remains a single candidate-level table.

The design-to-predicted-interface match is reported in two complementary ways. `design_interface_match_score_25models` is the mean geometric-interface recurrence of all designed mutations across the 25 AF3 models; `design_mutations_interface_ge80_*_25models` reports designed mutations that are interfacial in at least 20/25 models. `best_af3_prolif_design_any_contact_*` and `best_af3_prolif_design_nonvdw_*` report how many designed mutations actually contact another chain in the delivered best AF3 CIF, with the latter restricted to hydrophobic, hydrogen-bond, salt-bridge, pi-stacking, or cation-pi interactions. `best_af3_prolif_interface_designed_position_fraction` is the reverse coverage: the fraction of all ProLIF interface positions that are designed mutation sites.

All CIF coordinates are AF3 predictions, not Protenix structures. Computational predictions do not establish the solution-state oligomer experimentally.
