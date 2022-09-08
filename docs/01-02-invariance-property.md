## The Invariance Property

The parameter $\pi$ has a nice interpretation--it's a probability or the expected fraction of 1s in the long-run. However, the model parameters might not always have nice interpretation. (See the "shape" parameters of the beta distribution.) Fortunately, it's easy to transform the ML estimates of the model parameters into ML estimates of a quantity of interest.

Suppose we obtain an ML estimate $\hat{\theta}$ of a parameter $\theta$. But we also (or instead) want to estimate a transformation $\tau(\theta)$. The we can estimate $\tau(\theta)$ by applying the transformation $\tau$ to the ML estimate $\hat{\theta}$, so that $\widehat{\tau(\theta)} = \hat{\tau} = \tau(\hat{\theta})$.

### Example: Bernoulli Odds

Suppose that we want an ML estimator of the *odds* of getting a top for the toothpaste cap problem. We already used ML to estimate the *probability* $\pi$ of getting a top and came up with $\frac{8}{150} \approx 0.053$. We can directly transform a probability into odds using $\text{odds} = \frac{\pi}{1 - \pi}$. This has a nice interpretation: odds = 2 means that a top is twice as likely as not; odds = 0.5 means that a top is half as likely as not. 

In our case, we can plug our ML estimate of $\pi$ into the transformation to obtain the ML estimate of the odds.
$$
\begin{aligned}
\widehat{\text{odds}} &= \frac{\hat{\pi}}{1 - \hat{\pi}} \\
& = \frac{\frac{8}{150}}{1 - \frac{8}{150}} \\
& = \frac{\frac{8}{150}}{\frac{150}{150} - \frac{8}{150}} \\
& = \frac{\frac{8}{150}}{\frac{142}{150}} \\
& = \frac{8}{142} \\
& \approx 0.056
\end{aligned}
$$
This means that tops are about 0.06 times as likelihood as not-tops. Inverted, you're about $\frac{142}{8} \approx 18$ times more likely to not get a top than get a top.

### Example: Poisson SD 

In this example, we use real data from Hultman, Kathman, and Shannon (2013). They are interested in civilian casualties during civil wars. They write: 

> To gauge the effectiveness of peacekeeping, we explore all intrastate armed conflicts in sub-Saharan Africa from 1991 to 2008 with monthly observations. Conflicts are identified using the Uppsala Conflict Data Program/Peace Research Institute, Oslo (UCDP/PRIO) Armed Conflict Dataset v.4–2010 (Gleditsch et al. 2002; Harbom and Wallensteen 2009), which employs a threshold of 25 battle deaths per year. The dataset covers 36 conflicts, 12 of which have a PKO present at some time. Consistent with previous research, we add two years of observations to the end of each conflict episode, as the theoretical processes associated with victimization may continue after the cessation of hostilities (Cunningham, Gleditsch, and Salehyan 2009).

Below are a random sample of 250 observations from their 3,972 monthly observations.


```r
civilian_casualties <- c(0, 0, 0, 0, 0, 13, 0, 0, 61, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 19, 0, 0, 12, 0, 0, 4, 147, 0, 934, 0, 0,
                          42, 0, 24, 124, 0, 1, 0, 0, 0, 145844, 0, 0, 44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 2, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          7971, 0, 0, 0, 0, 72, 0, 40, 0, 0, 444, 0, 0, 0, 0, 48, 109, 33, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 15, 41, 0, 0, 0, 0, 84, 0, 34, 0, 0, 0,
                          0, 0, 0, 0, 1, 0, 15, 0, 0, 15, 0, 104, 0, 24, 0, 0, 104, 0, 0, 4, 0, 0, 0, 0,
                          0, 12, 41, 0, 0, 37, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 12, 0, 4, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 8, 21, 0, 0, 0, 0, 25, 0, 0, 0,
                          3, 0, 0, 27, 0, 0, 576, 3, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0,
                          0, 0, 0, 94, 42, 0, 30, 0, 2, 12, 0, 0, 5, 5 ) 
```

We can estimate a single-parameter Poisson model to estimate a mean $\lambda$ and a rate $\frac{1}{\lambda}$. In the case of the Poisson model, the ML estimate $\hat{lambda}$ of $\lambda$ is $\text{avg}(y)$.


```r
ml_est <- mean(civilian_casualties)
print(ml_est, digits = 3)
```

```
## [1] 630
```

The mean is a nice, interpretable parameter, but we might want also want the SD. For the Poisson distribution, the variance equals the mean, so $\text{Var}(y) = \text{E}(y) = \lambda$. Therefore, the SD is $\sqrt{\lambda}$.  


```r
# ML estimate of SD
ml_est <- sqrt(630)
print(ml_est, digits = 2)
```

```
## [1] 25
```

This is the ML estimate of the SD of the data, and it carries all the properties of ML estimators. We're using the invariance property to move from the mean to the SD by a simple transformation.

### Example: Beta Mean and Variance

Now let's see an example of the beta distribution $Y \sim \text{beta}(\alpha, \beta)$. The beta distribution does not have parameters that are easily interpretable in terms of mean and variance. Instead, it has two "shape" parameters $\alpha$ and $\beta$ that are in tension---one pulls the distribution to the left and the other pulls the distribution to the right.

For this example, I use opinion data from the 50 states from Barrilleaux and Rainey (2014). You can find the data here: https://github.com/carlislerainey/aca-opinion/blob/master/Data/mrp_est.csv

To make these data suitable for the beta distribution, I rescaled the observations from a percent to a proportion that ranges from 0 to 1.


```r
br <- tibble::tribble(
  ~state_abbr, ~prop_favorable_aca,
         "AL",   0.382711108911823,
         "AK",   0.374428493677838,
         "AZ",   0.396721609154912,
         "AR",   0.361623814680961,
         "CA",   0.560999240847165,
         "CO",   0.450011650633043,
         "CT",   0.522239143634457,
         "DE",   0.524637037667977,
         "DC",   0.853595690161985,
         "FL",    0.47022917052716,
         "GA",   0.460216990024346,
         "HI",    0.61965456264517,
         "ID",   0.282992730179373,
         "IL",   0.550517975187469,
         "IN",   0.421854785281297,
         "IA",   0.454007062646206,
         "KS",   0.394817640911206,
         "KY",   0.336156662764729,
         "LA",   0.425588396620569,
         "ME",   0.472319257331465,
         "MD",   0.583719023711148,
         "MA",   0.531871146279692,
         "MI",   0.509096426714406,
         "MN",   0.497981331879903,
         "MS",   0.468038078521612,
         "MO",   0.420161837905426,
         "MT",   0.351773944902139,
         "NE",   0.365225584190989,
         "NV",   0.459026605256376,
         "NH",    0.43886275738451,
         "NJ",   0.531656835425683,
         "NM",   0.528461049175538,
         "NY",     0.6010574821094,
         "NC",   0.452240849305449,
         "ND",   0.367690453757597,
         "OH",   0.456298880813516,
         "OK",   0.309578750918355,
         "OR",   0.455832591683007,
         "PA",    0.45819440292365,
         "RI",   0.536978574569609,
         "SC",   0.444870259057071,
         "SD",   0.377170366708612,
         "TN",   0.368615233253355,
         "TX",   0.428407014559672,
         "UT",   0.248496577141183,
         "VT",   0.553042362822573,
         "VA",   0.470739058046787,
         "WA",   0.496133477680592,
         "WV",   0.295062675817918,
         "WI",   0.489912969415965,
         "WY",   0.263567780036879
  )
```

Now let's find the ML estimates of the two shape parameters of the beta distribution.


```r
# obtain ml estimates
log_lik_fn <- function(par = c(2, 2), y) {
  a <- par[1]  # pulling these out makes the code a bit easier to follow
  b <- par[2]
  log_lik_i <- dbeta(y, shape1 = a, shape2 = b, log = TRUE)
  log_lik <- sum(log_lik_i)
  return(log_lik)
}
opt <- optim(par = c(2, 2), fn = log_lik_fn, y = br$prop_favorable_aca,
             control = list(fnscale = -1))
ml_est <- opt$par
print(ml_est, digits = 3)
```

```
## [1]  9.56 11.49
```

The mean is given by $\frac{\alpha}{\alpha + \beta}$ and the variance is given by $\frac{\alpha\beta}{(\alpha + \beta)^2(\alpha + \beta + 1)}$.

We can use the invariance property to obtain ML estimates of the mean and variance using our ML estimates of $\alpha$ and $\beta$.


```r
a <- ml_est[1]
b <- ml_est[2]

a/(a + b)  # mean
```

```
## [1] 0.4542508
```

```r
(a * b)/((a + b)^2 * (a + b + 1))  # var
```

```
## [1] 0.01123986
```

It's worth noting that these correspond closely, *but not exactly* to the observed mean and variance.


```r
mean(br$prop_favorable_aca)
```

```
## [1] 0.4524527
```

```r
var(br$prop_favorable_aca)
```

```
## [1] 0.01073633
```
