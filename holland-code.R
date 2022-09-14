
```{r}
# # load data from holland
# holland <- haven::read_dta("data/Enforcement.dta") %>%
#   #select(operations, lower, vendors, budget, population, city) %>%
#   glimpse()
# 
# # fit poisson regression model
# f <- operations ~ lower + vendors + budget + population
# fit <- glm(f, data = holland, family = poisson, 
#            subset = city == "santiago")
# 
# # parametric bs for coefficients
# n_bs <- 100  # should be 2000 or more
# coef_bs <- matrix(nrow = n_bs, ncol = length(coef(fit)))
# names(coef_bs) <- names(coef(fit))
# for (i in 1:n_bs) {
#   lambda_hat <- predict(fit, type = "response")
#   y_bs <- rpois(length(lambda_hat), lambda = lambda_hat)
#   fit_bs <- update(fit, formula = y_bs ~ .)
#   coef_bs[i, ] <- coef(fit_bs)
# }
# 
# # compute quantiles
# cis <- apply(coef_bs, 2, quantile, probs = c(0.05, 0.95))
# 
# # put results into a nice little data frame
# ci_df <- tibble(var = names(coef(fit)),
#                 est = coef(fit),
#                 lwr = cis["5%", ],
#                 upr = cis["95%", ])

```

