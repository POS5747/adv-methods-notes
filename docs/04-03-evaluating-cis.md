## Evaluating Confidence Intervals



### Coverage

Before we discuss these three intervals, let's review how we **evaluate** intervals. How do we know if a particular method works well? 

We evaluate confidence intervals in terms of their coverage: a $100(1 - \alpha)\%$ confidence interval, should capture the parameter $100(1 - \alpha)\%$ of the time under repeated sampling. That is, if we imagine repeating the study over-and-over (in the usual frequentist sense), then $100(1 - \alpha)\%$ of the confidence intervals contain the true parameter.

### Monte Carlo Simulation to Assess Coverage

#### A Simple Example

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
## [1] 0.8705
```

This simulation demonstrates that this simple interval captures the parameter $\lambda = 10$ in about 90% of repeated samples. This interval is *slightly* too narrow, because we should really use the $t$-interval here due to the small sample size. 

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
## [1] 0.9006
```

With this criterion in mind, let's consider three types of confidence intervals that we can use in the context of maximum likelihood estimation.

#### The Parametric Bootstrap

We can evaluate the parametric bootstrap using a similar Monte Carlo approach. However, this can be **confusing** because we have two types of simulation happening.

1. Monte Carlo simulation to evaluate the CI. We're having our computer conduct the same "study" over and over to compute the long-run, frequentist properties of the CI.
2. Parametric bootstrap to compute each CI. To compute each CI, we're simulating many fake outcome variables $y^{bs}$ from the fitted parametric distribution and refitting the model to obtain new estimates $\hat{\beta}^{bs}$ that we then summarize to find a single CI.

To start, let's create the data-generating process or "probability model" or population. The data are Bernoulli and the model is logit. 


```r
# model parameters
n <- 100
b0 <- 0
b1 <- 0.5
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
n_bs <- 2000  # should be 2000 or more
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
##           [,1]      [,2]        [,3]
## 5%  -0.6790475 0.9415808 -1.37164658
## 95%  0.3439898 2.7477924  0.09866435
```

Now let's simulate *many* "observed" data sets and compute the confidence intervals for each using the parametric bootstrap.


```r
# do a monte carlo simulation to compute the coverage for parametric bs
# note: this is not *at all* optimized for speed
n_mc_sims <- 100
ci_list <- list()
for (i in 1:n_mc_sims) {
  # simulate the "observed" data for one "study"
  y_obs <- rbinom(n, size = 1, prob = p)
  data <- data.frame(y_obs, x1, x2)
  fit <- glm(y_obs ~ x1 + x2, data = data, family = binomial)
  # parametric bs for coefficients
  n_bs <- 100  # should be 2000 or more
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
  # put cis into data frame; these are the intervals for a single "study"
  ci_list[[i]] <- tibble(coef_name = colnames(cis),
                         true = c(b0, b1, b2),
                  lwr = cis["5%", ],
                  upr = cis["95%", ], 
                  bs_id = i)
}

# combine the intervals for the many "studies" in a single data frame
ci_df <- bind_rows(ci_list)

# compute the coverage for each parameter
ci_df %>%
  mutate(captured = lwr < true & upr > true) %>%
  group_by(coef_name) %>%
  summarize(coverage = mean(captured), 
            se_hat = sd(captured)/sqrt(n()))
```

```
## # A tibble: 3 × 3
##   coef_name   coverage se_hat
##   <chr>          <dbl>  <dbl>
## 1 (Intercept)     0.83 0.0378
## 2 x1              0.87 0.0338
## 3 x2              0.82 0.0386
```

Thus the parametric bootstrap works well, the coverage is about 90%. But notice that we know the parametric model here, because we created the "observed" data ourselves and matched the GDP exactly with the parametric bootstrap. 

