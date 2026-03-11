# Metropolis-vs-Cluster-Monte-Carlo-for-2D-Potts-Model-Sampling-Efficiency-Across-Phase-Transition
Final Project For Intro to Numerical Methods (18.335, MIT)

## Requirements
- **MATLAB** R2019b or newer (no toolboxes required).
- Place all `.m` files in the **same folder**.

## File Overview

### `Experiment_1.m` — *scan_equal_work.m*
Scans temperature $T$ and compares **energy** $u(T)$ and **heat capacity** $c(T)$ between **Metropolis** and **Wolff** under **equal work** (work ≔ number of spin flips).

- **Key parameters (top of file):**
  - `N_list`, `q_list`: lattice sizes and Potts $q$ (paired by index; scalars allowed).
  - `W_perT_metro_list`, `W_perT_wolff_list`: target sweeps per $T$ (1 sweep = $N^2$ flips).
  - `T_values`: temperature grid.
  - `scan_burn_sweeps`: burn-in sweeps dropped at each $T$.
  - `repeats`: independent repeats per setting.

**Default parameter values are chosen to be good for quick experiments and sanity checks, although results will not be that accurate. For more accurate result, refer to the parameter setting in the experiment section of the paper.**


- **Outputs:**
  - Two figures per setting: **u(T)** and **c(T)** with error bars (mean ± s.e.).
  - Console logs with run settings.
- **Local helpers (end of file):**
  - `pick_k`, `drop_burn`
  - `collect_sweeps_equal_work_metro(...)`: single-spin Metropolis; accumulates flips until reaching the work target; records energy once per sweep.
  - `collect_sweeps_equal_work_wolff(...)`: Wolff cluster updates; **counts cluster size as flips**; records energy once per accumulated $N^2$ flips.



---

### `Experiment_2.m` — *dyn_compare_at_Tc.m*
At $T_c = 1/\log(1+\sqrt q)$, compares **dynamics**: integrated autocorrelation time (**IACT**), effective sample size (**ESS**), and **Cost/ESS** for Metropolis vs. Wolff.

- **Key parameters (top of file):**
  - `N_list`, `q_list` (paired or scalars),
  - `burn_sweeps`, `keep_sweeps`, `maxlag_dyn`.

**Default parameter values are chosen to be good for quick experiments and sanity checks.**

- **Outputs:**
  - Console summary per repeat/setting: `IACT`, `ESS`, total flips, `Cost/ESS`.
- **Local helpers:**
  - `metropolis_sweeps_series(...)`: records energy once per sweep after burn-in.
  - `wolff_sweeps_series(...)`: records energy once per accumulated $N^2$ flips after burn-in.
  - `iact_ips(x, maxlag)`: positive-sequence estimate of IACT; returns truncation lag `K`.

---

### `calculateH.m`
Computes total Hamiltonian with **periodic boundary conditions**:
$$
H = -J\sum_{\langle i,j\rangle}\mathbf 1\{\sigma_i=\sigma_j\} - h\sum_i \sigma_i.
$$
Returns the scalar energy (accounts for edge double counting).

### `calculateDeltaH.m`
Computes local energy **difference** $\Delta H$ when a single spin at $(i,j)$ is changed to `value`.  
Used by Metropolis proposals for efficient acceptance tests.

---

## How to Run
1. Put `Experiment_1.m`, `Experiment_2.m`, `calculateH.m`, `calculateDeltaH.m` in the same directory.
2. In MATLAB:
   ```matlab
   % Temperature scan with equal work:
   run('Experiment_1.m')

   % Dynamics Comparison at Tc (IACT/ESS/Cost per ESS):
   run('Experiment_2.m')
   
3. Adjust the parameter blocks at the top of each script as needed (lattice size, q, work per T, burn/keep sweeps, etc.).

## Important Notes / Conventions
- **One sweep =** $N^2$ **flips.** Metropolis counts 1 flip per proposal; Wolff counts cluster size as flips. `Experiment_1.m` compares the two methods at **matched work**.
- **Critical temperature:** $T_c = 1/\log(1+\sqrt{q})$ with $k_B=J=1$.
