
# Week 4: Confidence Intervals

This week, we expand our confidence interval toolkit. We have three core methods:


1. parametric bootstrap, which can be used directly for coefficients or quantities of interest.
1. nonparametric bootstrap, which can be used directly for coefficients or quantities of interest.
1. Wald confidence interval for coefficients, extended to quantities of interest using the delta method.



## Coverage

Before we discuss these three intervals, let's review how we **evaluate** intervals. How do we know if a particular method works well? 

We evaluate confidence intervals in terms of their coverage: a $100(1 - \alpha)\%$ confidence interval, should capture the parameter $100(1 - \alpha)\%$ of the time under repeated sampling. That is, if we imagine repeating the study over-and-over (in the usual frequentist sense), then $100(1 - \alpha)\%$ of the confidence intervals contain the true parameter.

As an example, let's consider the usual 90% confidence interval for the mean: 95% CI = $[\text{avg}(y) - 1.64 \times \hat{\text{SE}}, \text{avg}(y) + 1.64 \times \hat{\text{SE}}]$, where $\hat{\text{SE}} = \frac{\text{SD}(y)}{\sqrt{n}}$. We learned in an earlier class that this interval should capture the population average in about 90% of repeated trials. For out purposes, the "population average" refers to the mean parameter of some probability distribution.

Let's let the unknown distribution be $Y \sim \text{Poisson}(\lambda = 10)$. The "population" mean hear is $E(Y) = \lambda = 10$. Now let's use a Monte Carlo simulation to evaluate this particular interval. For this study, let's use a small sample size of 15 observations.


```r
# number of MC simulations (i.e., repeated trials)
n_mc_sims <- 10000

# contains for lower and upper bounds of 90% cis
lwr <- numeric(n_mc_sims)
upr <- numeric(n_mc_sims)

# mc simulations
for (i in 1:n_mc_sims) {
  y <- rpois(15, lambda = 10)
  se_hat <- sd(y)/sqrt(length(y))
  lwr[i] <- mean(y) - 1.64*se_hat
  upr[i] <- mean(y) + 1.64*se_hat
}

# combine results into a data frame
mc_sims <- tibble(iteration = 1:n_mc_sims,
                  lwr, upr) %>%
  mutate(captured = lwr < 10 & upr > 10)

# compute the proportion of simulations that capture the parameter
mean(mc_sims$captured)
```

```
## [1] 0.8729
```

This simulation demonstrates that this simple $z$-interval captures the parameter $\lambda = 10$ in about 87% of repeated samples. This interval is *slightly* too narrow, because we should really use the $t$-interval here due to the small sample size. 

The simulation below shows that this interval has better coverage (i.e., closer to 90%).


```r
# number of MC simulations (i.e., repeated trials)
n_mc_sims <- 10000

# contains for lower and upper bounds of 90% cis
lwr <- numeric(n_mc_sims)
upr <- numeric(n_mc_sims)

# mc simulations
for (i in 1:n_mc_sims) {
  y <- rpois(15, lambda = 10)
  se_hat <- sd(y)/sqrt(length(y))
  lwr[i] <- mean(y) - qt(.95, df = length(y) - 1)*se_hat
  upr[i] <- mean(y) + qt(.95, df = length(y) - 1)*se_hat
}

# combine results into a data frame
mc_sims <- tibble(iteration = 1:n_mc_sims,
                  lwr, upr) %>%
  mutate(captured = lwr < 10 & upr > 10)

# compute the proportion of simulations that capture the parameter
mean(mc_sims$captured)
```

```
## [1] 0.8985
```

With this criterion in mind, let's consider three types of confidence intervals that we can use in the context of maximum likelihood estimation.

## The Parametric Bootstrap

We've already seen the parametric bootstrap, but a brief review is worthwhile.

To do compute a confidence interval using the parametric bootstrap, do the following: 

1. Approximate $f(y; \theta)$ with $\hat{f} = f(y; \hat{\theta})$. Simulate a new outcome $y^{\text{bs}}$ from the estimated distribution. 
1. Re-compute the estimate of interest $\hat{\theta}^{\text{bs}}$ or $\hat{\tau}^{\text{bs}}$ using the bootstrapped outcome variable $y^{\text{bs}}$ rather than the observed outcome $y$.
1. Repeat 1 and 2 many times (say 2,000) to obtain many bootstrapped estimates. To obtain the 95% confidence interval, take the 2.5th and 97.5th percentiles of the estimates. In general, to obtain a $100(1 - \alpha)\%$ confidence interval, $\frac{\alpha}{2}$th and $(1 - \frac{\alpha}{2})$th percentiles. This is known as the percentile method.

The parametric bootstrap is a powerful, general tool to obtain confidence intervals for estimates from parametric models, *but it relies heavily on the parametric assumptions*.

We can evaluate the parametric bootstrap using a similar Monte Carlo approach. However, this can be **confusing** because we have two types of simulation happening.

1. Monte Carlo simulation to evaluate the CI. We're having our computer conduct the same "study" over and over to compute the long-run, frequentist properties of the CI.
2. Parametric bootstrap to compute each CI. To compute each CI, we're simulating many fake outcome variables $y^{bs}$ from the fitted parametric distribution and refitting the model to obtain new estimates $\hat{\beta}^{bs}$ that we then summarize to find a single CI.

To start, let's create the data-generating process or "probability model" or population. The data are Bernoulli and the model is logit. 


```r
# model parameters
n <- 50
b0 <- 1
b1 <- 1
b2 <- -0.5

# data
x1 <- rnorm(n, mean = 0, sd = 0.5)
x2 <- rbinom(n, size = 1, prob = 0.5)

# probability of success; Pr(y | x)
Xb <- b0 + b1*x1 + b2*x2
p <- plogis(Xb)
```

Now let's simulate *just one* "observed" data set and compute the confidence intervals for that data set using the parametric bootstrap.


```r
# simulate *one* sample and compute the 90% ci using parametric bs
y_obs <- rbinom(n, size = 1, prob = p)
data <- data.frame(y_obs, x1, x2)
fit <- glm(y_obs ~ x1 + x2, data = data, family = binomial)

# parametric bs for coefficients
n_bs <- 100  # should be 2000 or more
coef_bs <- matrix(nrow = n_bs, ncol = length(coef(fit)))
names(coef_bs) <- names(coef(fit))
for (i in 1:n_bs) {
  p_hat <- predict(fit, type = "response")
  y_bs <- rbinom(length(p_hat), size = 1, prob = p_hat)
  fit_bs <- update(fit, formula = y_bs ~ .)
  coef_bs[i, ] <- coef(fit_bs)
}

# compute quantiles
apply(coef_bs, 2, quantile, probs = c(0.05, 0.95))
```

```
##         [,1]     [,2]        [,3]
## 5%  1.344434 0.938166 -3.66211932
## 95% 5.119518 4.598202 -0.01051164
```

Now let's simulate *many* "observed" data sets and compute the confidence intervals for each using the parametric bootstrap.


```r
n_mc_sims <- 25
ci_list <- list()
for (i in 1:n_mc_sims) {
  # simulate the "observed" data for one "study"
  y_obs <- rbinom(n, size = 1, prob = p)
  data <- data.frame(y_obs, x1, x2)
  fit <- glm(y_obs ~ x1 + x2, data = data, family = binomial)
  
  # parametric bs for coefficients
  n_bs <- 25  # should be 2000 or more
  coef_bs <- matrix(nrow = n_bs, ncol = length(coef(fit)))
  colnames(coef_bs) <- names(coef(fit))
  for (j in 1:n_bs) {
    p_hat <- predict(fit, type = "response")
    y_bs <- rbinom(length(p_hat), size = 1, prob = p_hat)
    fit_bs <- update(fit, formula = y_bs ~ .)
    coef_bs[j, ] <- coef(fit_bs)
  }
  
  # compute quantiles
  cis <- apply(coef_bs, 2, quantile, probs = c(0.05, 0.95))
  
  # put results into data frame
  ci_list[[i]] <- tibble(coef_name = colnames(cis),
                         true = c(b0, b1, b2),
                  lwr = cis["5%", ],
                  upr = cis["95%", ], 
                  bs_id = i)
}

ci_df <- bind_rows(ci_list)

ci_df %>%
  mutate(captured = lwr < true & upr > true)
```

```
## # A tibble: 75 × 6
##    coef_name    true    lwr    upr bs_id captured
##    <chr>       <dbl>  <dbl>  <dbl> <int> <lgl>   
##  1 (Intercept)   1    0.422  2.17      1 TRUE    
##  2 x1            1    0.623  4.30      1 TRUE    
##  3 x2           -0.5 -3.29  -0.521     1 FALSE   
##  4 (Intercept)   1    0.658  3.03      2 TRUE    
##  5 x1            1   -0.463  1.62      2 TRUE    
##  6 x2           -0.5 -3.06   0.492     2 TRUE    
##  7 (Intercept)   1    0.219  1.51      3 TRUE    
##  8 x1            1   -0.429  1.52      3 TRUE    
##  9 x2           -0.5 -1.18   0.957     3 TRUE    
## 10 (Intercept)   1    1.25   3.00      4 FALSE   
## # … with 65 more rows
```

### Example: Toothpaste Cap Problm

The code below implements the parametric bootstrap for the toothpaste cap problem. For 2,000 iterations, it draws 150 observations from $Y \sim \text{Bernoulli}(\hat{\pi} = \frac{8}{150})$. For each iteration, it computes the ML estimate of $\pi$ for the bootstrapped data set. Then it computes the percentiles to obtain the confidence interval.


```r
n_bs <- 2000
bs_est <- numeric(n_bs)  # a container for the estimates
for (i in 1:n_bs) {
  bs_y <- rbinom(150, size = 1, prob = 8/150)
  bs_est[i] <- mean(bs_y)
}
print(quantile(bs_est, probs = c(0.025, 0.975)), digits = 2)  # 95% ci
```

```
##  2.5% 97.5% 
## 0.020 0.093
```

We leave an evaluation of this confidence interval (i.e., Does it capture $\theta$ 95% of the time?) to later in the semester.
