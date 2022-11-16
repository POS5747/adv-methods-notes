



## Count Models

### Poisson Regression


```r
# read data
keep <- c("osvAll", "troopLag", "policeLag", "militaryobserversLag",
          "brv_AllLag", "osvAllLagDum", "incomp", "epduration", "lntpop",
          "conflict_id")
hks <- read_csv("data/hks.csv") %>%
  select(keep) %>%
  na.omit() %>%
  glimpse()
```

```
## Rows: 3972 Columns: 10
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## dbl (10): osvAll, troopLag, policeLag, militaryobserversLag, brv_AllLag, osv...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```
## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
## ℹ Please use `all_of()` or `any_of()` instead.
##   # Was:
##   data %>% select(keep)
## 
##   # Now:
##   data %>% select(all_of(keep))
## 
## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
```

```
## Rows: 3,746
## Columns: 10
## $ osvAll               <dbl> 4, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ troopLag             <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ policeLag            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ militaryobserversLag <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ brv_AllLag           <dbl> 0, 138, 2428, 30, 850, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ osvAllLagDum         <dbl> 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ incomp               <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2…
## $ epduration           <dbl> 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1…
## $ lntpop               <dbl> 10.88525, 10.88525, 10.88525, 10.88525, 10.88525,…
## $ conflict_id          <dbl> 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 7…
```

```r
# fit models
f <- osvAll ~ troopLag + policeLag + militaryobserversLag + 
  brv_AllLag + osvAllLagDum + incomp + epduration + 
  lntpop

# poisson model
fit_pois <- glm(f, data = hks, family = poisson)

# predictive distribution
observed_data <- hks %>%
  mutate(type = "observed", 
         linpred_hat = predict(fit_pois, type = "link"))

sim_list <- list()
for (i in 1:5) {
  sim_list[[i]] <- observed_data %>%
    mutate(osvAll = rpois(nrow(observed_data), 
                          lambda = exp(observed_data$linpred_hat)),
           type = paste0("simulated #", i))
}
gg_data <- bind_rows(sim_list) %>%
  bind_rows(observed_data) %>%
  glimpse()
```

```
## Rows: 22,476
## Columns: 12
## $ osvAll               <dbl> 324, 339, 1179, 285, 432, 189, 181, 159, 185, 175…
## $ troopLag             <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ policeLag            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ militaryobserversLag <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ brv_AllLag           <dbl> 0, 138, 2428, 30, 850, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ osvAllLagDum         <dbl> 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ incomp               <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2…
## $ epduration           <dbl> 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1…
## $ lntpop               <dbl> 10.88525, 10.88525, 10.88525, 10.88525, 10.88525,…
## $ conflict_id          <dbl> 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 7…
## $ type                 <chr> "simulated #1", "simulated #1", "simulated #1", "…
## $ linpred_hat          <dbl> 5.701372, 5.756500, 7.017963, 5.651498, 6.088931,…
```

```r
gg1 <- ggplot(gg_data, aes(x = linpred_hat, y = osvAll + 1)) + 
  geom_point(alpha = 0.1, shape = 21, size = 0.3) + 
  facet_wrap(vars(type)) + 
  scale_y_log10()

gg2 <- ggplot(gg_data, aes(x = troopLag, y = osvAll + 1)) + 
  geom_point(alpha = 0.3, shape = 21, size = 0.3) + 
  facet_wrap(vars(type)) + 
  scale_y_log10() + 
  geom_smooth(se = FALSE)

library(patchwork)
```

```
## Warning: package 'patchwork' was built under R version 4.1.2
```

```r
gg1 + gg2
```

```
## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'
```

<img src="07-02-count-models_files/figure-html/unnamed-chunk-2-1.png" width="672" />
### Negative Binomial Regression


```r
# estimate negative binomial model
fit_nb <- MASS::glm.nb(osvAll ~ troopLag + policeLag + militaryobserversLag + 
              brv_AllLag + osvAllLagDum + incomp + epduration + 
              lntpop, 
            data = hks, init.theta = 5, control = glm.control(epsilon = 1e-12, maxit = 2500, trace = FALSE))

# predictive distribution
observed_data <- hks %>%
  mutate(type = "observed", 
         linpred_hat = predict(fit_nb, type = "link"))

sim_list <- list()
for (i in 1:5) {
  sim_list[[i]] <- observed_data %>%
    mutate(osvAll = rnbinom(nrow(observed_data), 
                            mu = predict(fit_nb, type = "response"), 
                            size = fit_nb$theta),
           type = paste0("simulated #", i))
}
gg_data <- bind_rows(sim_list) %>%
  bind_rows(observed_data) %>%
  glimpse()
```

```
## Rows: 22,476
## Columns: 12
## $ osvAll               <dbl> 0, 0, 1457, 265, 174, 0, 0, 0, 0, 2, 0, 0, 0, 0, …
## $ troopLag             <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ policeLag            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ militaryobserversLag <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ brv_AllLag           <dbl> 0, 138, 2428, 30, 850, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ osvAllLagDum         <dbl> 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ incomp               <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2…
## $ epduration           <dbl> 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1…
## $ lntpop               <dbl> 10.88525, 10.88525, 10.88525, 10.88525, 10.88525,…
## $ conflict_id          <dbl> 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 7…
## $ type                 <chr> "simulated #1", "simulated #1", "simulated #1", "…
## $ linpred_hat          <dbl> 5.351759, 5.448656, 7.065306, 5.371268, 5.949796,…
```

```r
gg1 <- ggplot(gg_data, aes(x = linpred_hat, y = osvAll + 1)) + 
  geom_point(alpha = 0.1, shape = 21, size = 0.3) + 
  facet_wrap(vars(type)) + 
  scale_y_log10()

gg2 <- ggplot(gg_data, aes(x = troopLag, y = osvAll + 1)) + 
  geom_point(alpha = 0.3, shape = 21, size = 0.3) + 
  facet_wrap(vars(type)) + 
  scale_y_log10() + 
  geom_smooth(se = FALSE)

library(patchwork)
gg1 + gg2
```

```
## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'
```

<img src="07-02-count-models_files/figure-html/unnamed-chunk-3-1.png" width="672" />


```r
AIC(fit_pois, fit_nb)
```

```
##          df        AIC
## fit_pois  9 2139137.24
## fit_nb   10   12540.31
```

### Zero-Inflated Negative Binomial Regression









