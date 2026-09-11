# OVA C4 final 20 package

## Contents

- `cif/`: 20 AF3 coordinate files, exactly one CIF per candidate. For each sequence this is the highest-ipTM model among its fixed 5 seeds × 5 samples (25 AF3 models).
- `all20_sequences.fasta`: all 20 homotetramer subunit sequences; each CIF contains four identical copies of the corresponding sequence.
- `tables/01_final_candidates_AF3.tsv`: primary sequence, mutation and aggregate AF3 table.
- `tables/02_AF3_per_seed_summary.tsv`: five-seed summary.
- `tables/03_AF3_per_model_scores.tsv`: all 500 AF3 model scores and source paths.
- `tables/04_mutation_interface_recurrence.tsv`: mutation/interface recurrence over 25 AF3 models per candidate.
- `tables/05_Protenix_v2_seed31_ranked.tsv`: independent Protenix-v2 ranking of all 20 candidates.
- `tables/06_Protenix_v2_seed31_AF3_order.tsv`: the same Protenix-v2 results in AF3 rank order.
- `tables/07_ProLIF_chain_pair_summary.tsv`: AF3 interface chemistry by chain pair.
- `tables/08_ProLIF_interface_residue_hotspots.tsv`: AF3 interface residue hotspots.
- `tables/09_ProLIF_top_crossmodel_consensus.tsv`: AF3/Protenix interface consensus for the original AF3 top candidate.
- `RESULT.md`: concise methods, criteria and interpretation.

All candidate CIF files are AF3 results, not Protenix coordinates. These are computational predictions and do not establish the solution-state oligomer experimentally.
