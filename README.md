# Heterogeneity and Global Climate Action — replication code

Python code accompanying the paper *Heterogeneity and Global Climate Action* (evolutionary political-economy integrated assessment model of heterogeneous countries' climate action), submitted to the *Journal of Economic Dynamics and Control*.

Maintainer: Giorgio Ricchiuti (University of Florence) — www.grarchive.net

## Requirements

Python 3 with `numpy`, `scipy`, `pandas`, `matplotlib`, `networkx` and `openpyxl` (the last is needed by `pandas` to read the `.xlsx` file).

```
pip install numpy scipy pandas matplotlib networkx openpyxl
```

Run every script from the repository root, so that the data files below are found. Scripts save figures as `.png` in the working directory.

## The model in brief

Each period, countries choose whether to adopt climate mitigation legislation. With `x_t` the relative participation (`x = 0` means half of the countries act), `E_t` net emissions and `gamma` the degree of homogeneity of preferences:

```
x_{t+1} = tanh( gamma * (beta_x * x_t + beta_e * E_t - mu) / 2 )
E_{t+1} = E_t * (1 + alpha_0 - alpha_1 * x_t)
```

In the IAM version `gamma` evolves with climate damages (DICE and Weitzman damage functions). Baseline calibration used throughout: `beta_x = 1.1`, `beta_e = 0.18`, `mu = 9.77`, `alpha_0 = 0`, `alpha_1 = 0.02`, `E_0 = 36.8`, `x_0 = -0.97`.

## Files

### Section 3 — deterministic (continuum) model

| File | Content |
|---|---|
| `section3_baseline_model.py` | Baseline model (fixed heterogeneity). Imported by the other Section 3 scripts. |
| `section3_baseline_model_alpha.py` | Same model with autonomous emissions growth (`alpha_0 ≠ 0`). |
| `section3_withIAM.py` | Model with carbon cycle and damages (DICE, Weitzman). |
| `section3_scenario_model.py` | Scenario simulations. |
| `section3_scenario_gamma.py` | Scenarios for different values of the degree of homogeneity `gamma`. |
| `section3_scenario_betax.py` | Scenarios for different values of the peer-pressure parameter `beta_x` (e.g. 1.1, 3, 15). |
| `section3_phase_plots.py` | Phase plots for different `gamma` and `beta_x`. |
| `section3_bifurcation_gamma.py` | Bifurcation diagram of `x_t` in `gamma`. |
| `section3_bifurcation_betax.py` | Bifurcation diagram of `x_t` in `beta_x` (Appendix, Figure 18, left). |
| `section3_bifurcation_mu.py` | Bifurcation diagram of `x_t` in `mu` (Appendix, Figure 18, right). |
| `section3_data_comparison.py` | Simulated vs. observed participation and net emissions, 1989–2021. |
| `appendix_apha.py` | Appendix robustness check on autonomous emissions growth (`alpha_0`). |

### Section 4 — agent-based model with Monte Carlo

| File | Content |
|---|---|
| `section4_agent_baseline.py` | Agent-based model with 194 countries (equal size). |
| `section4_montecarlo_baseline.py` | Monte Carlo runs of the baseline agent-based model. |
| `section4_agent_montecarlo_IAM.py` | Agent-based model with damages (IAM), Monte Carlo. |
| `section4_agent_montecarlo_countrysize.py` | Heterogeneous country size (GDP weights). |
| `section4_agent_montecarlo_sizeemissions.py` | Heterogeneous country size and emissions. |
| `network.py` | Builds the interaction network (fully connected) and the adjacency matrix. |
| `dataset.py` | Loads country weights and emissions from `countries_size_emissions.csv`. |

### Data

| File | Content |
|---|---|
| `countries_size_emissions.csv` | Country-level GDP weights (World Bank, 1989–2019 average) and GHG emissions, used in Section 4. |
| `xt_GHG_series.xlsx` | Observed participation `x_t` (from the Climate Change Laws of the World database) and net GHG emissions (Our World in Data), 1989–2021. |

## Data availability

The two data files above are included. The country-level panel used for the logit estimation of `beta_x`, `beta_e` and `mu` is built from third-party sources and is not redistributed here; it can be reconstructed from:

- Climate Change Laws of the World database (Grantham Research Institute and Sabin Center for Climate Change Law) — mitigation laws
- Our World in Data — GHG emissions
- World Bank — GDP
- Worldwide Governance Indicators (World Bank) — control of corruption
- ND-GAIN Country Index (University of Notre Dame) — vulnerability to climate damages

See the data table in the paper's appendix for the variable definitions.

## Contact

Questions about the code: Giorgio Ricchiuti, www.grarchive.net
