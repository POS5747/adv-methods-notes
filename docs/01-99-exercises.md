

## Exercises

From this week, you should be able to...

1. Starting with a given distribution (pdf or pmf), find the log-likihood function.
1. Optimize the log-likelihood function analytically or numerically (i.e., using `optim()`).
1. Use the invariance property of ML estimator to transform parameter estimates into estimates of quantities of interest.
1. Use the parametric bootstrap to compute confidence intervals.
1. Describe a "sampling distribution" and illustrate it with a computer simulation.
1. Describe the concepts of "bias" and "consistency." Are ML estimates unbaised? Consistent? Under what conditions?
1. Use the predictive distribution to evaluate models.

### Questions About the Exponential Distribution

::: {.exercise #ml-exponential}
Suppose we collect $N$ random samples $x = \{x_1, x_2, ..., x_N\}$ and model each draw as a random variable $X \sim \text{exponential}(\lambda)$ with pdf $f(x_n | \lambda) = \lambda e^{-\lambda x_n}$. 

1. Find the maximum likelihood estimator of $\lambda$.
1. Perform a Monte Carlo simulation to assess the bias in the ML estimator of $\lambda$ you found above. Use 100,000 Monte Carlo simulations. Estimate the Monte Carlo error as $\frac{\text{SD of estimates}}{\sqrt{\text{number of MC simulations}}}$. Try a small sample size (e.g., $N = 5$ and a large (e.g., $N = 1,000$). Demonstrate analytically that the estimate is biased.
1. Is the ML estimator of $\lambda$ consistent? Why or why not?
1. We interpret the parameter $\lambda$ as a "rate." Find the ML estimate of the mean, which is the reciprocal of the rate. Is this estimate unbiased? Consistent?
:::

<details><summary>Solution</summary>
1. The math follows the Poisson example closely. However, the solution is the inverse--$\hat{\lambda} = \frac{N}{\sum_{n = 1}^N x_n } = \frac{1}{\text{avg}(x)}$.
</details>

::: {.exercise #unnamed-chunk-1}
Suppose a data set `x <- c(0.306, 0.023, 0.0471, 0.042, 0.227)`. Model this data set using an exponential distribution and estimate the rate $\lambda$ using maximum likelihood. Find the estimates in two ways. First, compute the ML estimates using the analytical solution you found above. Second, derive the log-likelihood function, plot it, and use `optim()` to find the maximum. Show that the two solutions agree.
:::

::: {.exercise #unnamed-chunk-2}
Load the `cancer` data frame from the survival package in R using `data(cancer, package="survival")`. Model the variable `time` using an exponential distribution. Estimate both the rate and the mean. Use the parametric bootstrap to obtain a 95% confidence interval for each. Use the predictive distribution to evaluate the fit of the model. Do the simulated data sets seems to match the observed data sets?
:::

<!--- # remove for exam
::: {.exercise #ml-govt-duration}
Obtain the data for King, Alta, Burns, and Laver (2008) from Dataverse: https://doi.org/10.7910/DVN/CVJPAN. Find the file `coalcold.tab`. The variable `DURAT` gives the number of months that the government lasted. Model these durations as an exponential distributions. 

1. Estimate the rate $\lambda$ and mean (reciprocal of the rate). Use the parametric bootstrap to obtain a 95% confidence interval for each.
1. Simulate fake data from the predictive distribution and graphically compare these draws to the observed values. How do the observed data deviate from the model?
:::
--->

::: {.exercise #unnamed-chunk-3}
Obtain the data for Barrilleaux and Rainey (2014) from GitHub: https://github.com/carlislerainey/aca-opinion/blob/master/Data/mrp_est.csv. Find the file `mrp_est.csv`. The variable `percent_supporting_expansion` gives the the percent of each state that supports expanding Medicaid under the ACA. 

1. Model the distribution of *proportions* (you'll need to transform the variable from a percent to a proportion) as a beta distribution. Estimate the two shape parameters using maximum likelihood. Transform these estimates into estimates of the mean and SD. Compute a 95% confidence interval for the mean and SD using the parametric bootstap.
1. Simulate fake data from the predictive distribution and graphically compare these draws to the observed values. How do the observed data deviate from the model? Compare the the observed and simulated distributions in several ways, such as histograms, ecdf plots, two-number summaries, five-number summaries, etc. Look for ways that the simulated data sets consistently deviate from the observed.
:::

::: {.exercise #unnamed-chunk-4}
Suppose a discrete uniform distribution from 0 to $K$. The pdf is $f(x; K) = \frac{1}{K}$ for $x \in \{0, 1, ..., K\}$. Suppose I have three samples from the distribution: 676, 459, and 912. Find the ML estimate of $K$.
:::

::: {.exercise #unnamed-chunk-5}
**DeGroot and Schervish, q. 9, p. 425.** Suppose a distribution $f(x; \theta) = \theta x^{\theta - 1}$ for $0 < x < 1$ and $\theta > 0$. Find the ML estimator of $\theta$.
:::

::: {.exercise #unnamed-chunk-6}
Remember that the estimate $\hat{\lambda} = \text{avg}(y)$ for the Poisson distribution is unbiased. Remember that $E(y) = \lambda$ and $\text{SD}(y)$ is $\sqrt{\lambda}$. By the invariance properity, the ML estimator of $SD(y) = \sqrt{\hat{\lambda}}$. Use a Monte Carlo simulation to assess whether (and how much) this estimator of the SD is biased. Be sure to experiment with the sample size of y (e.g., N = 5, N = 200, etc.), but use a large number of Monte Carlo simulations (e.g., 100,000).
:::


