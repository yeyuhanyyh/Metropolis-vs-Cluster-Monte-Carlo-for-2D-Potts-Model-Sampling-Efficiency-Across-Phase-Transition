# Metropolis vs Cluster Monte Carlo for 2D Potts Model Sampling Efficiency Across Phase Transition

This repository contains MATLAB code for comparing **single-spin Metropolis sampling** and **Wolff cluster Monte Carlo** for the **two-dimensional q-state Potts model** on a square lattice.

The project focuses on two related questions:

1. How accurately do Metropolis and Wolff estimate thermodynamic observables across the phase transition when the computational work is matched?
2. How do the two samplers compare in terms of dynamical efficiency at the critical temperature?

This repository accompanies a final project for *Intro to Numerical Methods (18.335, MIT)*.

---

## Repository Structure

```text
Metropolis-vs-Cluster-Monte-Carlo-for-2D-Potts-Model-Sampling-Efficiency-Across-Phase-Transition/
├── codes/
│   ├── Experiment_1.m
│   ├── Experiment_2.m
│   ├── calculateH.m
│   └── calculateDeltaH.m
├── Report.pdf
├── final_project_requirement fall2025.pdf
└── README.md
```

---

## Model

We study the **2D q-state Potts model** on an `N × N` square lattice with periodic boundary conditions.

Each spin takes values in

```math
\sigma_{ij} \in \{1,2,\dots,q\}.
```

In the experiments implemented in this repository, the simulations are run in the **zero-field setting** with `J = 1` and `h = 0`.

The Hamiltonian is

```math
H(\sigma) = -J \sum_{\langle x,y\rangle} \delta_{\sigma_x,\sigma_y},
```

where the sum is over nearest-neighbor pairs and `\delta` is the Kronecker delta.

The Gibbs distribution at temperature `T` is

```math
\pi_T(\sigma) \propto \exp\!\bigl(-\beta H(\sigma)\bigr),
\qquad
\beta = 1/T.
```

For the 2D Potts model on the square lattice with `J = k_B = 1`, the critical temperature is

```math
T_c = \frac{1}{\log(1+\sqrt{q})}.
```

---

## Observables

The main thermodynamic quantities studied in this project are the **energy density** and the **heat capacity**.

The energy density is

```math
u(T) = \frac{1}{N^2}\,\mathbb{E}[H].
```

The heat capacity is estimated from energy fluctuations as

```math
c(T) = \frac{\beta^2}{N^2}\,\operatorname{Var}(H).
```

These are exactly the quantities plotted in the temperature-scan experiment.

---

## Algorithms

### 1. Metropolis single-spin updates

The Metropolis sampler proposes a change at one lattice site at a time.

A proposed move from `\sigma` to `\sigma'` is accepted with probability

```math
a(\sigma \to \sigma')
=
\min\!\left\{1,\exp(-\beta \Delta H)\right\},
```

where `\Delta H = H(\sigma') - H(\sigma)`.

In the implementation:

- one proposal counts as **one flip-equivalent** of work;
- energy is recorded once per sweep, where one sweep means `N^2` flip-equivalents.

### 2. Wolff cluster updates

The Wolff method grows a same-spin cluster and recolors it in one move.

A nearest neighbor of a cluster site is added to the cluster with probability

```math
p_{\mathrm{add}} = 1 - \exp(-\beta J).
```

In the implementation:

- the computational work of one Wolff update is measured by the **cluster size**;
- cluster size is counted in flip-equivalents so that the total work can be fairly compared with Metropolis.

---

## Experiment 1: Temperature Scan Under Matched Work

`Experiment_1.m` scans a temperature grid and compares Metropolis and Wolff under **equal work**.

### Goal

For each temperature, estimate:

- energy density,
- heat capacity,
- uncertainty across repeated runs.

### Work normalization

To make the comparison fair, this script measures work in **spin flips**.

A single sweep is defined as:

```math
1 \text{ sweep} = N^2 \text{ flip-equivalents}.
```

The two methods are normalized as follows:

- **Metropolis:** each single-spin proposal counts as 1 flip;
- **Wolff:** each cluster update counts as the size of the cluster.

Energy is recorded once per accumulated sweep-equivalent work.

### What the script does

- scans temperature from high `T` to low `T`,
- uses an annealed warm start across the temperature grid,
- discards a prescribed number of burn-in sweeps at each temperature,
- repeats the experiment several times,
- plots mean values with standard-error bars.

### Main outputs

For each parameter setting, the script produces:

- a plot of energy density versus temperature,
- a plot of heat capacity versus temperature,
- console logs summarizing the run settings.

---

## Experiment 2: Dynamical Comparison at the Critical Temperature

`Experiment_2.m` compares the dynamical efficiency of Metropolis and Wolff at the critical temperature.

### Goal

At

```math
T = T_c = \frac{1}{\log(1+\sqrt{q})},
```

the script records an energy time series and estimates:

- integrated autocorrelation time,
- effective sample size,
- cost per effective sample.

### Integrated autocorrelation time

Let `\rho(k)` denote the lag-`k` autocorrelation of the recorded energy sequence. The code uses the estimator

```math
\tau_{\mathrm{int}} = 1 + 2 \sum_{k=1}^{K} \rho(k),
```

where `K` is chosen by positive-sequence truncation.

### Effective sample size

If the recorded time series has length `n`, then the code defines

```math
\mathrm{ESS} = \frac{n}{\tau_{\mathrm{int}}}.
```

### Cost per effective sample

The code also reports

```math
\mathrm{Cost/ESS} = \frac{\text{total flip-equivalents}}{\mathrm{ESS}}.
```

This gives a fair efficiency comparison between the two algorithms after normalizing by computational work.

### Main outputs

For each setting, the script prints:

- estimated IACT,
- estimated ESS,
- total work in flip-equivalents,
- Cost/ESS.

---

## File Overview

### `codes/Experiment_1.m`

Temperature scan under matched work.

This script compares Metropolis and Wolff across a range of temperatures and plots:

- `u(T)` for energy density,
- `c(T)` for heat capacity.

### `codes/Experiment_2.m`

Dynamics comparison at the critical temperature.

This script evaluates:

- integrated autocorrelation time,
- effective sample size,
- cost per effective sample.

### `codes/calculateH.m`

Computes the total Hamiltonian of the current Potts configuration with periodic boundary conditions.

### `codes/calculateDeltaH.m`

Computes the local energy change for a proposed single-spin update.  
This is used in the Metropolis acceptance step.

---

## Requirements

- MATLAB R2019b or newer
- No additional toolboxes are required

---

## How to Run

Open MATLAB, enter the `codes` directory, and run the desired script.

### Run Experiment 1

```matlab
cd codes
run('Experiment_1.m')
```

### Run Experiment 2

```matlab
cd codes
run('Experiment_2.m')
```

---

## Parameters You May Want to Change

At the top of each script, you can adjust parameters such as:

- lattice size `N`,
- number of Potts states `q`,
- temperature grid,
- burn-in length,
- number of kept sweeps,
- number of repeats,
- work budget for each method.

The default settings are suitable for quick tests and sanity checks. For more accurate experiments, see the parameter choices discussed in `Report.pdf`.

---

## Notes and Conventions

- Periodic boundary conditions are used throughout.
- All work comparisons are normalized in **flip-equivalents**.
- In this repository, one sweep always means `N^2` flip-equivalents.
- The Wolff method is compared against Metropolis under this shared work convention.
- The experiments in the scripts are run in the zero-field case.

---

## Included Documents

- `Report.pdf` — final project report
- `final_project_requirement fall2025.pdf` — original project handout

---

## Author

Yuhan Ye
