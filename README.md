![GitHub release](https://img.shields.io/github/r-package/v/mrzdcmps/changeofevidence.svg)
[![GitHub last commit](https://img.shields.io/github/last-commit/mrzdcmps/changeofevidence.svg)](https://github.com/mrzdcmps/changeofevidence/commits/master)


# Change of Evidence
This package provides functions to test for a change of Bayesian evidence over time which can be used to examine volatile effects.

## Installation
```
  # Install remotes if necessary
  if (!requireNamespace("remotes")) install.packages("remotes")
  # Get package from Github
  remotes::install_github("mrzdcmps/changeofevidence")
```

## Functionality

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

**Returns:** A `seqbf` object containing sequential BF values, test statistics, p-values, effect size estimates (δ), and raw data.

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
Convenience wrapper that automatically generates the appropriate simulations based on a `seqbf` object. It inspects the test type and extracts parameters (sample size, prior, alternative) automatically.

**Primary interface (recommended):**
```r
result <- bfttest(rnorm(50, 0.5), mu = 0)
sims <- simcreate(result, n.sims = 1000)
```

**Parameters:**
- `x`: A `seqbf` object (recommended) — test type and parameters are inferred automatically
- `n.sims`: Number of simulations to generate (default: 1000)
- `N`: Override the sample size extracted from the `seqbf` object
- `method`: `"pseudo"` (software RNG), `"files"` (quantum RNG), or `"quantis"` (hardware QRNG)
- `filespath`: Path to quantum random data file (if using `method = "files"`)
- `parallel`: Use parallel processing (default: `TRUE`)
- Override parameters: `mu`, `n_bits`, `p`, `rho`, `data_type`, `prior.loc`, `prior.r`, `alternative`, `nstart`

**Returns:** Dataframe with columns: `simid`, `index`, `raw`, `rw`, `density.rw`, `bf`, `density.bf`.

> **Note:** The old interface `simcreate(trials, mean.scores, ...)` is deprecated. Use the type-specific functions below instead.

#### `simcreate_bin()` - Binomial Test Simulations
Generates Monte Carlo simulations for binomial tests.

**Parameters:**
- `N`: Sample size (number of binary observations)
- `n.sims`: Number of simulations (default: 1000)
- `p`: Null hypothesis probability (default: 0.5)
- `method`: `"pseudo"`, `"files"`, or `"quantis"`
- `parallel`: Use parallel processing (default: `TRUE`)
- `nstart`: Minimum observations before first BF (default: 5)
- `prior.r`: Prior scale parameter (default: 0.1)
- `alternative`: `"two.sided"`, `"less"`, or `"greater"`

**Returns:** Dataframe with columns: `simid`, `index`, `raw`, `rw`, `density.rw`, `bf`, `density.bf`.

#### `simcreate_t()` - T-Test Simulations
Generates Monte Carlo simulations for t-tests with support for different data types.

**Parameters:**
- `N`: Sample size
- `n.sims`: Number of simulations (default: 1000)
- `mu`: Null hypothesis mean
- `data_type`: `"summed_bits"`, `"continuous"`, `"integer"`, or `"unknown"`
- `n_bits`: Number of bits to sum per observation (for `data_type = "summed_bits"`)
- `p`: Bit probability (default: 0.5)
- `int_min`, `int_max`: Integer range (for `data_type = "integer"`)
- `method`: `"pseudo"`, `"files"`, or `"quantis"`
- `parallel`: Use parallel processing (default: `TRUE`)
- `nstart`, `alternative`, `prior.loc`, `prior.r`: Test parameters

**Returns:** Dataframe with columns: `simid`, `index`, `raw`, `rw`, `density.rw`, `bf`, `density.bf`.

#### `simcreate_cor()` - Correlation Test Simulations
Generates Monte Carlo simulations for correlation tests.

**Parameters:**
- `N`: Sample size (number of paired observations)
- `n.sims`: Number of simulations (default: 1000)
- `rho`: Null hypothesis correlation (default: 0)
- `method`: `"pseudo"` or `"quantis"` (files not supported)
- `parallel`: Use parallel processing (default: `TRUE`)
- `nstart`: Minimum observations before first BF (default: 5)
- `prior.r`: Prior scale parameter (default: 0.353)

**Returns:** Dataframe with columns: `simid`, `index`, `raw`, `raw2`, `rw`, `bf`, `density.bf`.

#### `download_quantum_data()` - Download Quantum Random Data
Downloads pregenerated quantum random bits for use with `method = "files"`. After downloading once, the file location is remembered automatically.

**Parameters:**
- `path`: Directory to save the file (default: package extdata directory)
- `force`: Re-download even if file exists (default: `FALSE`)

**Example:**
```r
# Download once
download_quantum_data()

# Then use automatically in simulations
result <- bfbinom(rbinom(100, 1, 0.6))
sims <- simcreate(result, n.sims = 1000, method = "files")
```

### 4. Plotting Functions

#### `plotbf()` - Plot Sequential Bayes Factors
Visualizes BF trajectories over time with evidence strength annotations.

**Parameters:**
- `...`: One or more `seqbf` objects, BF vectors, or a single list of vectors
- `labels`: Names for multiple datasets (auto-generated if `NULL`)
- `sims.df`: Optional simulation data for background
- `sims.df.col`: Column in simulation dataframe to plot (default: `"bf"`)
- `color`: A single color or vector of colors for the BF line(s). A vector applies custom colors to multiple lines (overrides the default palette).
- `coordy`: Manual y-axis limits
- `label.x`: Override x-axis label (default: `"N"`)
- `show_annotations`: Show evidence strength labels (default: `TRUE`)

**Examples:**
```r
# Single result
plotbf(result)

# Multiple results with custom labels and colors
plotbf(bf1, bf2, labels = c("Experiment 1", "Experiment 2"), color = c("blue", "red"))

# With simulation background
plotbf(result, sims.df = sims)
```

#### `plotrw()` - Plot Random Walk
Plots random walk trajectories with confidence intervals.

**Parameters:**
- `...`: One or more numeric vectors of random walk values, or a single named list
- `labels`: Names for multiple random walks (auto-generated if `NULL`)
- `sims.df`: Optional simulation data
- `sims.df.col`: Column in simulation dataframe to plot (default: `"rw"`)
- `color`: A single color or vector of colors. When multiple series are plotted, a vector applies custom colors (overrides default palette).
- `coordy`: Manual y-axis limits
- `p`: Probability parameter for confidence bounds (default: 0.5)
- `n_bits`: Number of bits summed per trial (for deviation-based walks)

**Examples:**
```r
# Single random walk
plotrw(my_rw)

# Multiple random walks
plotrw(rw1, rw2, labels = c("Group A", "Group B"), color = c("steelblue", "coral"))

# With simulation comparison
plotrw(my_rw, sims.df = sims, n_bits = 10)
```

#### `plotfft()` - Plot FFT
Visualizes frequency spectra with 95% confidence intervals from simulations.

**Parameters:**
- `data`: FFT-transformed vector or `seqbf` object
- `sims.df`: Optional simulation data
- `sims.df.col`: Column to compare to (default: `"density.bf"`)
- `n.hz`: Number of frequencies to display (default: 50)
- `color`: Line color (default: `"black"`)

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

Note: BF<sub>01</sub> = 1 / BF<sub>10</sub>

**CoE p-value:**
The harmonic mean of p-values from MaxBF, Energy, and FFT tests. Values < 0.05 suggest the temporal pattern of evidence is unusual and may indicate a volatile effect rather than a stable effect or random fluctuation.



## How to use the package
[Vignette](https://mrzdcmps.github.io/changeofevidence/vignette.html)

## Random Files

You can download pregenerated quantum random bits directly from within R using `download_quantum_data()`. Alternatively, 10,000 data files each containing 1,000,000 random bits generated by a quantum random number generator (Quantis by idquantique) are available here:
https://osf.io/gs42z/files/
