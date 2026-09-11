# Invalid Volcengine wave `ky-20260901-558`–`565`

The eight jobs recorded in `volc_design_jobs_v1.tsv` are excluded from every
candidate pool, AF3-label set, ranking, and delivery artifact. Their runtime
commands omitted `OVA_IPTM90_EXPERIMENT`, so the generic worker fell back to
E155 and read the 19K/21B/21S/20A parents instead of E169's
`OVA-C4-073-BB1`–`BB4` AF3-medoid parents.

All eight jobs were stopped with the Volcengine `StopJob` API on 2026-09-01 UTC.
The corrected wave must use
`af3_pipeline/run_iptm90_tied_design_volc_e169_worker.sh`, which validates all
four anchor names before starting optimization.
