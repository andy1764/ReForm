# Reference interval calibration via conFormal prediction (ReForm)

**Author**: Andrew A. Chen, chenandr@musc.edu 

**Maintainer**: Andrew A. Chen, chenandr@musc.edu

**License**: Artistic License 2.0

ReForm calibrates fitted reference intervals to new samples using split 
conformalized quantile regression (Romano et al., 2019). ReForm requires a 
set of observations from the new sample, called the calibration set, which is 
used to separately adjust the lower and upper bounds of the reference interval.
Empirical results (Chen et al., 2025) demonstrate that ReForm performs well so 
long as a sufficiently sized calibration set is available (~40 for 90% 
reference intervals, ~80 for 95%, and ~120 for 99%).

Currently, our package supports standard regression models including:

1. Linear models (`lm`)
2. Generalized additive models for location, scale, and shape (GAMLSS) via [`gamlss2`](https://github.com/gamlss-dev/gamlss2)
3. Quantile regression via [`quantreg`](https://cran.r-project.org/web/packages/quantreg/index.html)
4. Quantile generalized additive models via [`qgam`](https://cran.r-project.org/web/packages/qgam/index.html)

This package is still a work-in-progress and will be continually supported and 
updated to include additional modeling options.

## 1. Installation
The latest version can be installed via `remotes` by running the following code

```
# install.packages("remotes")
remotes::install_github("andy1764/ReForm")
```

Then, you can load this package via

```
library(ReForm)
```

## 2. Usage
A toy example for calibrating a 90% reference interval using the `iris` dataset 
is provided below:

```
library(quantreg)
# construct 90% reference interval for Sepal.Length in iris virginica
ri <- refint(
  rq(Sepal.Length ~ Sepal.Width,
  data = iris[51:100,], tau = 0.05),
  rq(Sepal.Length ~ Sepal.Width,
  data = iris[51:100,], tau = 0.95)
)

# direct application in iris setosa
predict(ri, newdata = iris[1:10,])

# calibrate `ri` for iris setosa using 40 holdout observations
reformed_ri <- reform(ri, caldata = iris[11:50,])

# application of ReFormed reference interval
predict(reformed_ri, newdata = iris[1:10,])
```

For additional examples, please refer to the manual pages

```
help(package = "ReForm")
```

## 3. Citations
Please cite the following papers when using ReForm:

> Romano, Y., Patterson, E., & Candes, E. (2019). Conformalized quantile 
regression. Advances in neural information processing systems, 32.
>
> Chen, A., ...
