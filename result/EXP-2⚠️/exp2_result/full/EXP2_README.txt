EXP-2 (Runtime vs t) - plotting artifact

mode=paper
warmup_reps=1
min_repeats=5
t0=1000
error_mode=stdev (requested=auto)
paper_with_sample_phase=1
variants=sampling,enum_sampling,adaptive

raw_ok_rows_used(after warmup)=1450
wall_points_aggregated=290
points_dropped_by_min_repeats=0

Paper plots:
  exp2_paper_runtime_vs_t.(pdf|png)
  exp2_paper_delta_vs_t.(pdf|png)
  exp2_paper_sample_phase_vs_t.(pdf|png)  [if phases_json available]

Tables:
  exp2_ns_per_sample.csv

Paper method set:
  ours, interval_tree, r_tree, pbsm, tlsop, range_tree

