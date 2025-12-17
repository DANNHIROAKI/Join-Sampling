# Spatial-Join-Sampling

```bash
├── apps
│   ├── sjs_convert_realdata.cpp
│   ├── sjs_gen_dataset.cpp
│   ├── sjs_profile.cpp
│   ├── sjs_run.cpp
│   ├── sjs_sweep.cpp
│   └── sjs_verify.cpp
├── CMakeLists.txt
├── config
│   ├── default.json
│   └── sweeps
│       ├── sweep_alpha.json
│       ├── sweep_dim.json
│       ├── sweep_scale.json
│       └── sweep_t.json
├── data
│   ├── real
│   └── synthetic
├── docs
│   ├── 00_overview.md
│   ├── 01_data_format.md
│   ├── 02_experiment_protocol.md
│   ├── 03_baselines.md
│   ├── 04_extending_to_high_dim.md
│   └── 05_adding_real_datasets.md
├── include
│   └── sjs
│       ├── baselines
│       │   ├── aabb
│       │   │   ├── adaptive.h
│       │   │   ├── enum_sampling.h
│       │   │   └── sampling.h
│       │   ├── baseline_api.h
│       │   ├── interval_tree
│       │   │   ├── adaptive.h
│       │   │   ├── enum_sampling.h
│       │   │   └── sampling.h
│       │   ├── kd_tree
│       │   │   ├── adaptive.h
│       │   │   ├── enum_sampling.h
│       │   │   └── sampling.h
│       │   ├── ours
│       │   │   ├── adaptive.h
│       │   │   ├── enum_sampling.h
│       │   │   └── sampling.h
│       │   ├── pbsm
│       │   │   ├── adaptive.h
│       │   │   ├── enum_sampling.h
│       │   │   └── sampling.h
│       │   ├── range_tree
│       │   │   ├── adaptive.h
│       │   │   ├── enum_sampling.h
│       │   │   └── sampling.h
│       │   ├── rejection
│       │   │   ├── adaptive.h
│       │   │   ├── enum_sampling.h
│       │   │   └── sampling.h
│       │   ├── r_tree
│       │   │   ├── adaptive.h
│       │   │   ├── enum_sampling.h
│       │   │   └── sampling.h
│       │   ├── runners
│       │   │   ├── adaptive_runner.h
│       │   │   ├── enum_sampling_runner.h
│       │   │   └── sampling_runner.h
│       │   ├── sirs
│       │   │   ├── adaptive.h
│       │   │   ├── enum_sampling.h
│       │   │   └── sampling.h
│       │   ├── tlsop
│       │   │   ├── adaptive.h
│       │   │   ├── enum_sampling.h
│       │   │   └── sampling.h
│       │   └── tsunami
│       │       ├── adaptive.h
│       │       ├── enum_sampling.h
│       │       └── sampling.h
│       ├── core
│       │   ├── assert.h
│       │   ├── config.h
│       │   ├── logging.h
│       │   ├── rng.h
│       │   ├── stats.h
│       │   ├── timer.h
│       │   └── types.h
│       ├── data
│       │   └── synthetic
│       │       ├── clustered.h
│       │       ├── generator.h
│       │       ├── hetero_sizes.h
│       │       ├── stripe_ctrl_alpha.h
│       │       └── uniform.h
│       ├── geometry
│       │   ├── box.h
│       │   ├── embedding.h
│       │   ├── point.h
│       │   └── predicates.h
│       ├── index
│       │   ├── aabb_tree.h
│       │   ├── interval_tree.h
│       │   ├── kd_tree.h
│       │   ├── range_tree.h
│       │   └── r_tree.h
│       ├── io
│       │   ├── binary_io.h
│       │   ├── csv_io.h
│       │   ├── dataset.h
│       │   └── realdata_stub.h
│       ├── join
│       │   ├── join_enumerator.h
│       │   ├── join_oracle.h
│       │   ├── join_types.h
│       │   └── sweep_events.h
│       └── sampling
│           ├── alias_table.h
│           ├── rank_sampling.h
│           ├── sample_quality.h
│           └── weighted_choice.h
├── LICENSE
├── plot.ipynb
├── README.md
├── results
│   ├── raw
│   └── summary
├── src
│   ├── apps
│   │   ├── sjs_convert.cpp
│   │   ├── sjs_gen.cpp
│   │   ├── sjs_run.cpp
│   │   ├── sjs_sweep.cpp
│   │   └── sjs_verify.cpp
│   ├── baselines
│   │   ├── aabb
│   │   │   ├── adaptive.cpp
│   │   │   ├── enum_sampling.cpp
│   │   │   └── sampling.cpp
│   │   ├── baseline_factory_2d.cpp
│   │   ├── baseline_factory_2d.h
│   │   ├── baseline_names.cpp
│   │   ├── interval_tree
│   │   │   ├── adaptive.cpp
│   │   │   ├── enum_sampling.cpp
│   │   │   └── sampling.cpp
│   │   ├── kd_tree
│   │   │   ├── adaptive.cpp
│   │   │   ├── enum_sampling.cpp
│   │   │   └── sampling.cpp
│   │   ├── ours
│   │   │   ├── adaptive.cpp
│   │   │   ├── enum_sampling.cpp
│   │   │   └── sampling.cpp
│   │   ├── pbsm
│   │   │   ├── adaptive.cpp
│   │   │   ├── enum_sampling.cpp
│   │   │   └── sampling.cpp
│   │   ├── range_tree
│   │   │   ├── adaptive.cpp
│   │   │   ├── enum_sampling.cpp
│   │   │   └── sampling.cpp
│   │   ├── rejection
│   │   │   ├── adaptive.cpp
│   │   │   ├── enum_sampling.cpp
│   │   │   └── sampling.cpp
│   │   ├── r_tree
│   │   │   ├── adaptive.cpp
│   │   │   ├── enum_sampling.cpp
│   │   │   └── sampling.cpp
│   │   ├── runners
│   │   │   └── instantiations_2d.cpp
│   │   ├── sirs
│   │   │   ├── adaptive.cpp
│   │   │   ├── enum_sampling.cpp
│   │   │   └── sampling.cpp
│   │   ├── tlsop
│   │   │   ├── adaptive.cpp
│   │   │   ├── enum_sampling.cpp
│   │   │   └── sampling.cpp
│   │   └── tsunami
│   │       ├── adaptive.cpp
│   │       ├── enum_sampling.cpp
│   │       └── sampling.cpp
│   ├── CMakeLists.txt
│   ├── core
│   │   └── instantiations_2d.cpp
│   ├── data
│   │   ├── instantiations_2d.cpp
│   │   └── synthetic
│   │       ├── generator_factory_2d.cpp
│   │       └── instantiations_2d.cpp
│   ├── geometry
│   │   └── instantiations_2d.cpp
│   ├── index
│   │   ├── index_factory_2d.cpp
│   │   └── instantiations_2d.cpp
│   ├── io
│   │   ├── load_dataset_2d.cpp
│   │   ├── realdata_stub.cpp
│   │   └── write_results.cpp
│   ├── join
│   │   ├── build_events_2d.cpp
│   │   ├── instantiations_2d.cpp
│   │   └── oracle_2d.cpp
│   └── sampling
│       ├── instantiations_2d.cpp
│       └── quality_checks.cpp
├── tests
│   ├── CMakeLists.txt
│   ├── test_baselines_smoke.cpp
│   ├── test_event_sweep.cpp
│   ├── test_geometry.cpp
│   ├── test_join_oracle.cpp
│   └── test_sampling_quality.cpp
└── third_party
    ├── cxxopts
    └── rapidjson
```



## 1. RQs

- **RQ1 (Correctness / Sample Quality):**
  Does the output satisfy the requirement of i.i.d. uniform sampling with replacement over the join results?

- **RQ2 (End-to-end Efficiency):**
  Does our method significantly outperform the sampling versions of baselines in end-to-end efficiency when only $t$ samples are required?

- **RQ3 (Scalability vs. $\alpha_{\text{pair}}$, $N$, $t$):**
  Do the scalability trends with respect to N, \alpha, and t align with expectations (specifically, remaining insensitive to the join size $|J|$)?

- **RQ4 (Robustness):**
  Does the method remain robust under conditions of spatial skew, extreme aspect ratios, hotspot clustering, or the mixture of large objects (which result in **massive active sets**)?

- **RQ5 (Adaptivity):**
  Does the Adaptive variant achieve runtime performance close to the optimal **$\min(\text{Sampling}, \text{Enumerate+Sampling})$**?

- **RQ6 (Practical Utility):**
  Can the generated samples accurately estimate downstream statistics, such as the **distribution of overlap areas**?

## 2. Selected Results

#### Correctness (p-values + independence)

   ![image-20251217072501269](https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20251217072501269.png)

**Correctness verification on the Oracle dataset ($N=10^4, |J|=10^4$).** (a) Distribution of $p$-values from $\chi^2$ goodness-of-fit tests (25 trials). A uniform distribution over $[0,1]$ indicates unbiased sampling. (b) Lag-1 autocorrelation of sampled pair IDs. Points within the gray band (95% CI) indicate independence. (c) Duplicate rate evolution as sample size $t$ increases. The theoretical curve represents the expected collision rate for sampling with replacement.

#### Runtime vs $t$

   <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20251217072850734.png" width=500  />  

**Most scalable performance with respect to sample size $t$.** While enumeration-based methods (gray line) suffer from a high constant cost dictated by $|J|$, and rejection sampling (purple line) exhibits high variance due to low acceptance rates, our method maintains a consistently low latency with minimal slope.

#### Runtime vs $N$

   <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20251217075901519.png" alt="image-20251217075901519" width=500  /> 

**Scalability of end-to-end runtime with increasing dataset cardinality $N$ (Log-Log scale).** Experiments are conducted on SYN-U datasets with fixed sample size $t=10^4$ and fixed join selectivity $\alpha_{\text{pair}} \approx 10^{-6}$. Ours-Sampling (red) scales nearly linearly ($\approx O(N \log N)$), consistently outperforming tree-based baselines due to better memory locality. In contrast, Enum-Stream (gray) exhibits quadratic growth ($\approx O(N^2)$) as $|J|$ increases, becoming computationally prohibitive at $N \ge 3 \cdot 10^6$. Rejection (green) shows high variance and instability due to fluctuating acceptance rates.

#### Runtime vs $\alpha$ + acceptance rate

   <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20251217080157578.png" alt="image-20251217080157578" width=600  /> 

**Effect of join density $\alpha_{\text{pair}}$ on runtime and acceptance rate.** Experiments on SYN-CTRL datasets ($N=10^6, t=10^4$). (Top) Ours-Sampling (red) remains robustly efficient regardless of $\alpha$, whereas Enum-Stream (gray) degrades linearly with increasing $|J|$. Rejection (green) performs poorly at low densities due to extremely low acceptance rates. Ours-Adaptive (purple line) successfully traces the lower envelope, acting as an enumeration method in sparse regions and switching to sampling in dense regions. (Bottom) The acceptance rate $p_{\text{acc}}$ of the Rejection method is directly proportional to $\alpha$, explaining its inverse runtime trend. 

#### Adaptive phase diagram

   <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20251217081606526.png" alt="image-20251217081606526" width=600  /> 

**Phase diagram of the optimal sampling strategy in the $(\alpha_{\text{pair}}, t)$ parameter space.** The heatmap indicates the speedup of Ours-Sampling over Enumerate+Sampling (blue favors Enumeration, red favors Sampling). The dashed line represents the theoretical switching boundary where both methods have equal cost. Region I (low density, large $t$) favors enumeration due to the fixed overhead of sampling structures. Region II (high density or small $t$) heavily favors direct sampling. The hatched area on the right indicates the OOM Zone ($|J| > J^*_{\text{limit}}$), where enumeration becomes infeasible, making sampling the only viable option. Our Adaptive method automatically selects the optimal strategy based on these operating conditions.

#### Adaptive ratio $\rho$ 

   <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20251217082434270.png" alt="image-20251217082434270" width=600  />

**Heatmap of the Adaptive Ratio $\rho$ across the parameter space.** The ratio $\rho$ represents the overhead of our Adaptive strategy compared to the optimal choice between Enumerate+Sampling and Ours-Sampling (an oracle baseline). In most regions, $\rho \approx 1.0$, indicating negligible overhead. A narrow "transition band" (orange) appears where the costs of the two base strategies are similar ($T_{\text{Enum}} \approx T_{\text{Sample}}$), causing slight decision overhead ($\le 15\%$). A secondary vertical ridge near the memory threshold ($|J| \approx J^*$) reflects the cost of partially materializing results before switching to sampling. Overall, the method incurs minimal penalty while automatically capturing the optimal performance.