# tpt-aero-cli

The `tpt-aero` command: run VLM wing sweeps from a TOML case file and get CSV out.

```text
tpt-aero template wing > wing.toml     # commented starting point
tpt-aero wing wing.toml --out polar.csv
tpt-aero isa 0 5000 11000
```

Columns: `alpha_deg,cl,cdi,cm,cl_roll,cn,lift_n`. Unknown or misspelled case-file keys are rejected.
