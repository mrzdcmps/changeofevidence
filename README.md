![GitHub release](https://img.shields.io/github/r-package/v/mrzdcmps/changeofevidence.svg)
[![GitHub last commit](https://img.shields.io/github/last-commit/mrzdcmps/changeofevidence.svg)](https://github.com/mrzdcmps/changeofevidence/commits/master)


# Change of Evidence
This package provides functions to test for a change of Bayesian evidence over time which can be used to examine volatile effects. 

## Core Functions

### 1. Bayesian Sequential Testing

The package provides sequential Bayesian hypothesis tests that calculate Bayes Factors (BF) as data accumulates. These functions allow you to track evidence evolution over time.

#### `bfttest()` - Bayesian Sequential t-Test
Performs sequential Bayesian t-tests (parametric) or non-parametric tests (Wilcoxon/Mann-Whitney).

**Parameters:**
- `x`, `y`: Numeric vectors or formula for data input
- `formula`, `data`: Alternative formula interface (e.g., `response ~ group`)
- `parametric`: `TRUE` for t-test, `FALSE` for non-parametric (default: `TRUE`)
- `alternative`: `"two.sided"`, `"less"`, or `"greater"`
- `mu`: Null hypothesis value (default: 0)
- `prior.loc`, `prior.r`: Prior distribution parameters (location and scale)
- `nstart`: Minimum observations before first BF calculation (default: `"auto"`)
- `parallel`: Use parallel processing for speed (default: `TRUE`)

**Test types:**
- One-sample: `bfttest(x)`
- Paired samples: `bfttest(x, y)`
- Independent samples: `bfttest(response ~ group, data = df)`

**Returns:** A `seqbf` object containing sequential BF values, test statistics, p-values, and effect size estimates (δ).

#### `bfbinom()` - Bayesian Sequential Binomial Test
Tests binary data (0s and 1s) sequentially.

**Parameters:**
- `data`: Vector of binary data (0s and 1s)
- `p`: Null hypothesis probability of success (default: 0.5)
- `prior.r`: Prior scale parameter (default: 0.1)
- `alternative`: Direction of test
- `nstart`: Starting sample size (default: 5)

**Returns:** A `seqbf` object with sequential BF values and probability of success estimates.

#### `bfcor()` - Bayesian Sequential Correlation Test
Tests correlation between two continuous variables sequentially.

**Parameters:**
- `x`, `y`: Numeric vectors to correlate
- `alternative`: Test direction
- `prior.r`: Prior scale parameter (default: 0.1)
- `nstart`: Starting sample size (default: 5)

**Returns:** A `seqbf` object with sequential BF values and correlation coefficients.

#### `bfRobustness()` - Prior Robustness Analysis
Examines how sensitive your Bayes Factor is to different prior specifications.

**Parameters:**
- `x`: A `seqbf` object from `bfttest()`, `bfbinom()`, or `bfcor()`
- `informed`: `TRUE` tests multiple locations and widths; `FALSE` tests only widths (default: `TRUE`)
- `prior.loc`, `prior.r`: Custom prior ranges (optional)

**Returns:** A `bfRobustness` object showing BF values across different prior specifications.

### 2. Change of Evidence (CoE) Analysis

These functions test whether temporal patterns in Bayesian evidence are unusual compared to simulated null data.

#### `coe()` - Comprehensive Change of Evidence Analysis
Performs all three CoE tests (MaxBF, Energy, FFT) and computes an overall p-value.

**Parameters:**
- `data`: A `seqbf` object or numeric vector of Bayes Factors
- `sims.df`: Dataframe containing simulated null distributions (from `simcreate()`)

**Returns:** A `coe` object containing:
- Maximum BF results
- Energy analysis results
- FFT analysis results
- Combined CoE p-value (harmonic mean of individual p-values)

#### `maxbf()` - Maximum BF Analysis
Tests whether the highest observed BF is unusually extreme.

**Parameters:**
- `data`: A `seqbf` object or BF vector
- `sims.df`: Simulation dataframe

**Returns:** Maximum BF value, its position, and percentage of simulations with higher peaks.

#### `energybf()` - BF Energy Analysis
Calculates the "energy" of evidence (area under the BF curve minus area below BF=1).

**Parameters:**
- `data`: A `seqbf` object or BF vector
- `sims.df`: Simulation dataframe

**Returns:** Energy value and comparison to simulated distributions.

#### `ffttest()` - Frequency Analysis Test
Tests whether oscillatory patterns in the BF trajectory are unusual using Fast Fourier Transform.

**Parameters:**
- `data`: A `seqbf` object or FFT-transformed data (from `fftcreate()`)
- `sims.df`: Simulation dataframe

**Returns:** Amplitude sum and comparison to simulations.

#### `fftcreate()` - Create FFT
Transforms a BF time series into frequency domain using Fast Fourier Transform.

**Parameters:**
- `data`: A `seqbf` object or numeric vector

**Returns:** Vector of spectral densities (N/2 frequencies).

#### `simredo()` - Adjust Simulation Length
Recalculates FFT densities for simulations to match a different data length.

**Parameters:**
- `df`: Simulation dataframe
- `n`: Desired number of trials per simulation
- `rw`: Also recalculate Random Walk FFT (default: `TRUE`)

**Returns:** Adjusted simulation dataframe.

### 3. Simulation Generation

#### `simcreate()` - Generate Null Simulations
Creates Monte Carlo simulations of random data for CoE analysis.

**Parameters:**
- `trials`: Number of trials per simulation
- `n.sims`: Number of simulations to generate (default: 1000)
- `mean.scores`: If specified, sums bits to create normally distributed scores centered at this value
- `method`: `"pseudo"` (software RNG), `"files"` (quantum RNG files), or `"quantis"` (hardware QRNG)
- `filespath`: Path to random bit files (if using `method = "files"`)
- `parallel`: Use parallel processing (default: `TRUE`)
- `nstart`, `alternative`, `prior.loc`, `prior.r`, `p`: Test parameters matching your analysis

**Returns:** Dataframe with columns: `simid`, `index`, `raw`, `rw`, `density.rw`, `bf`, `density.bf`.

### 4. Plotting Functions

#### `plotbf()` - Plot Sequential Bayes Factors
Visualizes BF trajectories over time with evidence strength annotations.

**Parameters:**
- `...`: One or more `seqbf` objects or BF vectors
- `labels`: Names for multiple datasets
- `sims.df`: Optional simulation data for background
- `color`: Line color (default: `"black"`)

**Example:**
```r
plotbf(bf1, bf2, labels = c("Experiment 1", "Experiment 2"))
```

#### `plotrw()` - Plot Random Walk
Plots random walk trajectories with confidence intervals.

**Parameters:**
- `data`: Vector or list of random walk values
- `sims.df`: Optional simulation data
- `p`: Probability parameter for confidence bounds (default: 0.5)
- `n_bits`: Number of bits summed per trial (for deviation-based walks)

#### `plotfft()` - Plot FFT
Visualizes frequency spectra with 95% confidence intervals from simulations.

**Parameters:**
- `data`: FFT-transformed vector
- `sims.df`: Optional simulation data
- `n.hz`: Number of frequencies to display (default: 50)

#### `plotrobust()` - Plot Robustness Analysis
Shows how BF varies across different prior specifications.

**Parameters:**
- `data`: A `bfRobustness` object

## Interpretation

**Bayes Factors:**
- BF > 1: Anecdotal evidence
- BF > 3: Moderate evidence
- BF > 10: Strong evidence
- BF > 30: Very strong evidence
- BF > 100: Extreme evidence

Note: BF$_{01}$ = 1 / BF$_{10}$

**CoE p-value:**
The harmonic mean of p-values from MaxBF, Energy, and FFT tests. Values < 0.05 suggest the temporal pattern of evidence is unusual and may indicate a volatile effect rather than a stable effect or random fluctuation.



## Installation
```
  # Install remotes if necessary
  if (!requireNamespace("remotes")) install.packages("remotes")
  # Get package from Github
  remotes::install_github("mrzdcmps/changeofevidence")
```
    
## How to use the package
[Vignette](https://mrzdcmps.github.io/changeofevidence/vignette.html)
  
## Random Files

You can download 10,000 data files each containing 1,000,000 random bits generated by a quantum random number generator (Quantis by idquantique) as source for the monte carlo simulations here:
https://osf.io/gs42z/files/
