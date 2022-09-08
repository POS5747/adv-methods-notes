## Sampling Distribution



What's the most important concept in statistical inference? I don't know, but it could be **the sampling distribution**. For effect, let me back off the hedge.

> The most important concept in statistical inference is the **sampling distribution**.

To define a sampling distribution, you need to imagine repeating a study over and over. If each study has a random component (perhaps random sampling or random assignment to treatment and control), then the estimate will differ from study to study. The distribution of the estimates across the studies is called the sampling distribution.

### Example: The Toothpaste Cap Problem

For a given sample of 150 tosses, we recognize the the ML estimate $\hat{\pi}$ does not (usually) exactly equal the parameter $\pi$. Instead, the particular $\hat{\pi}$ that the study produces is draw from a distribution.

Let's illustrate that with a simulation. For these simulations, I suppose that we toss the toothpaste cap 150 times and the chance of a head is 5%.


```r
n_sims <- 10
ml_est <- numeric(n_sims)  # a container for the estimates
for (i in 1:n_sims) {
  y <- rbinom(150, size = 1, prob = 0.05)
  ml_est[i] <- mean(y)
}
print(ml_est, digits = 2)
```

```
##  [1] 0.040 0.040 0.047 0.047 0.080 0.027 0.060 0.053 0.053 0.067
```

As you can see, the ML estimates vary to from sample to sample--different data sets produce different ML estimates. We need a way to create a confidence interval that consistently captures $\theta$. 

If we repeat the simulations a large number of times, we can see an accuracy picture of the sampling distribution via histogram.


```r
n_sims <- 10000
ml_est <- numeric(n_sims)  # a container for the estimates
for (i in 1:n_sims) {
  y <- rbinom(150, size = 1, prob = 0.05)
  ml_est[i] <- mean(y)
}

gg_data <- data.frame(ml_est = ml_est)
ggplot(gg_data, aes(x = ml_est)) + 
  geom_bar()
```

<img src="01-04-bias-and-consistency_files/figure-html/unnamed-chunk-3-1.png" width="384" />

Many of our methods of evaluating an estimator are statements about the sampling distribution of that estimator. In general, we'd like the sampling distribution to be centered over the true parameter of interest and tightly dispersed.

## Bias

Imagine repeatedly sampling and computing the estimate $\hat{\theta}$ of the parameter $\theta$ for each sample. In this thought experiment, $\hat{\theta}$ is a random variable. We say that $\hat{\theta}$ is **biased** if $E(\hat{\theta}) \neq \theta$. We say that $\hat{\theta}$ is **unbiased** if $E(\hat{\theta}) = \theta$. We say that the **bias** of $\hat{\theta}$ is  $E(\hat{\theta}) - \theta$.

Importantly, **ML estimators are not necessarily unbiased**. Of the models we will see in this course, *most* are biased.

### Example: Bernoulli Distribution

For example, we can compute the bias of our ML estimator of $\pi$ in the toothpaste cap problem.

$$
\begin{aligned}
E\left[ \frac{k}{N}\right] &= \frac{1}{N} E(k) = \frac{1}{N} E  \overbrace{ \left( \sum_{n = 1}^N x_n \right) }^{\text{recall } k = \sum_{n = 1}^N x_n } = \frac{1}{N} \sum_{n = 1}^N E(x_n) = \frac{1}{N} \sum_{n = 1}^N \pi = \frac{1}{N}N\pi \\
&= \pi
\end{aligned}
$$

Thus, $\hat{\pi}^{ML}$ is an unbiased estimator of $\pi$ in the toothpaste cap problem.

We can use a Monte Carlo simulation to check this analytical result.


```r
set.seed(1234)
n_mc_sims <- 100000
pi_hat <- numeric(n_mc_sims)
for (i in 1:n_mc_sims) {
  y <- rbinom(150, size = 1, prob = 0.05)
  pi_hat[i] <- mean(y)
}

# expected value of pi-hat
mean(pi_hat)
```

```
## [1] 0.05006227
```

```r
# estimated monte carlo error
sd(pi_hat)/sqrt(n_mc_sims)
```

```
## [1] 5.631271e-05
```

But notice that the property of unbiasedness does not follow the estimate through transformation. Because the sample is relatively large in this case (150 tosses), the bias is small, but detectable with 100,000 Monte Carlo simulations


```r
odds_hat <- pi_hat/(1 - pi_hat)

# actual odds
0.05/(1 - 0.05)
```

```
## [1] 0.05263158
```

```r
# expected value of odds-hat
mean(odds_hat)
```

```
## [1] 0.05307323
```

```r
# estimated monte carlo error
sd(odds_hat)/sqrt(n_mc_sims)
```

```
## [1] 6.288517e-05
```

```r
# the z-statistic
(mean(odds_hat) - 0.05/0.95)/(sd(odds_hat)/sqrt(n_mc_sims))
```

```
## [1] 7.023072
```


### Example: Poisson Distribution

Using math almost identical to the toothpaste cap problem, we can show that the ML estimator $\hat{\lambda} = \text{avg}(x)$ is an unbiased estimator of $\lambda$. 

We can also illustrate the unbiasedness with a computer simulation.


```r
lambda <- 4.0      # the parameter we're trying to estimate
sample_size <- 10  # the sample size we're using in each "study"

n_mc_sims <- 10000  # the number of times we repeat the "study"
lambda_hat <- numeric(n_mc_sims)  # a container 
for (i in 1:n_mc_sims) {
  x <- rpois(sample_size, lambda = lambda)
  lambda_hat[i] <- mean(x)
}

# expected value of lambda-hat
mean(lambda_hat)
```

```
## [1] 3.99397
```

```r
# estimated monte carlo error
sd(lambda_hat)/sqrt(n_mc_sims)
```

```
## [1] 0.006300177
```

## Consistency

Imagine taking a sample of size $N$ and computing the estimate $\hat{\theta}_N$ of the parameter $\theta$. We say that $\hat{\theta}$ is a **consistent** estimator of $\theta$ if $\hat{\theta}$ converges in probability to $\theta$.

Intuitively, this means the following:

1. For a large enough sample, the estimator returns the exact right answer.
1. For a large enough sample, the estimate $\hat{\theta}$ does not vary any more, but collapses onto a single point and that point is $\theta$.

Under weak, but somewhat technical, assumptions that usually hold, ML estimators are consistent. 

Given that we always have finite samples, why is consistency valuable? In short, it's not valuable, directly. However, consistent estimators tend to be decent with small samples. 

But it does not follow that consistent estimators work well in small samples. However, as a rough guideline, consistent estimators work well for small samples. However, whether they actually work well in any particular situation needs a more careful investigation.

### Example: Illustrative

To illustrate the concept of consistency, consider this estimator of the population mean $\hat{\mu}^{\text{silly}} = \frac{\sum_{i = 1}^N x_i}{N + 10}$. While this estimator is biased, it is a consistent estimator.


```r
population <- c(1, 2, 3, 4, 5)
sample_sizes <- c(2, 5, 10, 100, 1000, 10000, 100000)
n_mc_sims <- 30  # for each sample size
results_list <- list()  # grow with each iteration; slow, but easy
for (i in 1:length(sample_sizes)) {
  ml_est_i <- numeric(n_mc_sims)
  for (j in 1:n_mc_sims) {
    x <- sample(population, sample_sizes[i], replace = TRUE)
    ml_est_i[j] <- sum(x)/(sample_sizes[i] + 10)
  }
  results_list[[i]] <- data.frame(sample_size = sample_sizes[i],
                             ml_est = ml_est_i)
}
results <- dplyr::bind_rows(results_list) 

ggplot(results, aes(x = sample_size, y = ml_est)) + 
  geom_hline(yintercept = mean(population)) + 
  geom_jitter() + 
  scale_x_log10()
```

<img src="01-04-bias-and-consistency_files/figure-html/unnamed-chunk-7-1.png" width="672" />



### Example: Bernoulli Odds

There are two ways to see consistency for the Bernoulli. First, unless our sample size is a multiple of 20, it is impossible to obtain an estimated odds of 0.05/(1 - 0.05). Second, in small samples, the ML estimate of the odds is biased. As the sample size increases, the bias shrinks and the estimates collapse toward (and eventually onto) the true value.



```r
sample_sizes <- c(2:100, 250, 400, 500, 750, 1000)
n_mc_sims <- 10  # for each sample size
results_list <- list()  # grow with each iteration; slow, but easy
for (i in 1:length(sample_sizes)) {
  ml_est_i <- numeric(n_mc_sims)
  for (j in 1:n_mc_sims) {
    x <- rbinom(sample_sizes[i], 1, prob = 0.05)
    pi_hat <- mean(x)
    ml_est_i[j] <- pi_hat/(1 - pi_hat)
  }
  results_list[[i]] <- data.frame(sample_size = sample_sizes[i],
                             ml_est = ml_est_i)
}
results <- dplyr::bind_rows(results_list) 

ggplot(results, aes(x = sample_size, y = ml_est)) + 
  geom_hline(yintercept = 0.05/(1 - 0.05)) + 
  geom_jitter(alpha = 0.5, shape = 19) + 
  scale_x_log10() 
```

<img src="01-04-bias-and-consistency_files/figure-html/unnamed-chunk-8-1.png" width="672" />

```r
ggplot(results, aes(x = sample_size, y = ml_est)) + 
  geom_hline(yintercept = 0.05/(1 - 0.05)) + 
  scale_x_log10() + 
  geom_smooth()
```

<img src="01-04-bias-and-consistency_files/figure-html/unnamed-chunk-8-2.png" width="672" />
