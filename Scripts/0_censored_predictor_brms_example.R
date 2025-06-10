### example by B Trueman 
### https://bentrueman.github.io/posts/2023-03-13-regression-models-with-censored-predictors/

library("tidyverse")
library("glue")
library("ggdist")
library("withr")
library("brms")
library("rstan")
library("viridis")
options(mc.cores = parallel::detectCores()) # for parallel MCMC chains

theme_set(
  theme_bw(14) + 
    theme(
      legend.position = "bottom",
      strip.background = element_blank()
    )
)

n <- 50 # total number of observations
lcl <- -1 # lower censoring limit

with_seed(1235, {
  data <- tibble(
    x = rnorm(n, 0, 2),
    y = 2 * x + rnorm(n)
  )
})

data$x_star <- data$x
data$x_star[c(3, 29)] <- NA # add missing values
data$cens_x_star <- replace_na(data$x_star < lcl, FALSE) # censoring indicator
data$x_star <- pmax(data$x_star, lcl) # left-censor values

# train/test split:
train <- 1:25 
data_train <- data[train, ]
data_test <- data[-train, ]


p1 <- data %>% 
  rowid_to_column() %>% 
  mutate(
    type = if_else(rowid %in% train, "Training data", "Test data"),
    missing = is.na(data$x_star)
  ) %>% 
  ggplot(aes(x, y, col = type)) + 
  scale_color_viridis_d() +
  geom_point(
    data = . %>% 
      filter(cens_x_star | missing),
    alpha = .4, shape = 16, size = 2
  ) +
  geom_point(
    data = . %>% 
      filter(!cens_x_star & !missing),
    shape = 16, size = 2
  ) +
  geom_hline(
    data = . %>% 
      filter(missing),
    aes(yintercept = y, col = type)
  ) +
  geom_segment(
    data = . %>%
      filter(cens_x_star),
    aes(x = -1, xend = -Inf, yend = y)
  ) + 
  labs(col = NULL)
p1


#The trick at this point is to treat left-censored `x` values as missing values with an upper bound equal to the censoring limit, as described in the [Stan manual](https://mc-stan.org/docs/2_18/stan-users-guide/censored-data.html). But since `brms` doesn't do this we'll have to code it ourselves. We'll need two functions: one to modify the data list and another to modify the code that `brms` passes to Stan:

modify_standata <- function(sdata, data, lcl, var) {
  
  if (length(lcl) != length(var)) stop("lengths of 'var' and 'lcl' must be equal")
  
  varstan <- str_remove_all(var, "_")
  
  for(i in seq_len(length(var))) {
    sdata[[paste0("Ncens_", varstan[i])]] <- sum(data[[paste0("cens_", var[i])]]) # number of left-censored
    # positions of left-censored:
    sdata[[paste0("Jcens_", varstan[i])]] <- as.array(seq_len(nrow(data))[data[[paste0("cens_", var[i])]]]) 
    sdata[[paste0("U_", varstan[i])]] <- lcl[i] # left-censoring limit
  }
  
  sdata
}

modify_stancode <- function(scode, var) {
  
  var <- str_remove_all(var, "_")
  
  for(i in seq_len(length(var))) {
    
    # modifications to data block:
    n_cens <- glue("int<lower=0> Ncens_{var[i]};  // number of left-censored")
    j_cens <- glue("int<lower=1> Jcens_{var[i]}[Ncens_{var[i]}];  // positions of left-censored")
    u <- glue("real U_{var[i]};  // left-censoring limit")
    # modifications to parameters block:
    y_cens <- glue("vector<upper=U_{var[i]}>[Ncens_{var[i]}] Ycens_{var[i]};  // estimated left-censored")
    # modifications to model block:
    yl <- glue("Yl_{var[i]}[Jcens_{var[i]}] = Ycens_{var[i]}; // add imputed left-censored values")
    
    scode <- scode %>%
      # modifications to data block:
      str_replace(
        "(data \\{\n(.|\n)*?)(?=\n\\})",
        paste(c("\\1", n_cens, j_cens, u), collapse = "\n  ")
      ) %>% 
      # modifications to parameters block:
      str_replace(
        "(parameters \\{\n(.|\n)*?)(?=\n\\})",
        paste(c("\\1\n  ", y_cens), collapse = "")
      ) %>% 
      # modifications to model block:
      str_replace(
        "(model \\{\n(.|\n)*?)(?=\n    mu_)",
        paste(c("\\1\n    ", yl), collapse = "")
      )
    
  }
  
  class(scode) <- "brmsmodel"
  
  scode
    
}


#We'll also need a function to impute missing and censored values for prediction using the fitted model:
  
impute <- function(data, model, var, mi = NULL, cens = NULL, id = NULL) {
  
  varstan <- str_remove_all(var, "_")
  
  censored <- as_draws_df(model) %>% 
    as_tibble() %>% 
    select(starts_with(paste0("Ycens_", varstan))) %>% 
    t()
  
  if (!is.null(id)) {
    censored <- censored[id, ]
  }
  
  missing <- as_draws_df(model) %>% 
    as_tibble() %>% 
    select(starts_with(paste0("Ymi_", varstan))) %>% 
    t()
  
  ndraws <- ncol(censored)
  
  map(
    seq_len(ndraws),
    \(x) {
      data_imputed <- data
      if (!is.null(mi)) {data_imputed[mi, var] <- missing[, x]}
      if (!is.null(cens)) {data_imputed[cens, var] <- censored[, x]}
      data_imputed
    }, .progress = TRUE
  )
  
}

#Now we're ready to set up and run the model. First we need a formula that defines a model for both `y` and `x`. Since we don't have anything else to go on, the model for `x` is intercept only, but we should be able to do better than that in many real applications. The calls to `mi()` define missing values and how they're handled during model fitting---see the `brms` [vignette](https://cran.r-project.org/web/packages/brms/vignettes/brms_missings.html) for details.

bform <- bf(y | mi() ~ mi(x_star)) + 
  bf(x_star | mi() ~ 1) +
  set_rescor(FALSE)

#Next, we modify the data and code passed to Stan... 

sdata <- make_standata(bform, data = data_train) %>% 
  modify_standata(data_train, lcl, "x_star")
scode <- make_stancode(bform, data = data_train) %>% 
  modify_stancode("x_star")

#... and fit (and save) the model using `rstan::stan`.

stanseed <- 1257
model_rstan <- stan(
  model_code = scode,
  data = sdata,
  sample_file = "model-censored-x", # save model as CSVs
  seed = stanseed
)


#Now we can generate posterior predictions along a regular sequence of predictor values, using the posterior draws for the slope and intercept:

post_draws <- as_draws_df(model_rstan) # posterior draws:

xnew <- seq(min(data$x), max(data$x), length.out = 25) # regular sequence of x-values

post_imp <- map(
  seq_len(nrow(post_draws)),
  ~ xnew * post_draws$`bsp_y[1]`[.x] + post_draws$Intercept_y[.x]
) %>%
  do.call(rbind, .)

#Here is the fitted model---the shaded region represents a 95% credible interval on the posterior expectation.

preds <- tibble(
  x = xnew,
  yhat = apply(post_imp, 2, median),
  yhat_min = apply(post_imp, 2, \(x) quantile(x, .025)),
  yhat_max = apply(post_imp, 2, \(x) quantile(x, .975)),
  )

p1 + 
  geom_ribbon(
    data = preds,
    aes(x, ymin = yhat_min, ymax = yhat_max),
    alpha = .3,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = preds,
    aes(x, yhat),
    inherit.aes = FALSE
  )

#Let's compare the parameter estimates from this model with those of a simple linear regression model that doesn't account for censoring or missingness:

model_naive <- brm(
  y ~ x_star,
  data = data_train,
  file = "model-naive",
  file_refit = "on_change"
)
post_draws_naive <- as_draws_df(model_naive)

#The posterior draws of the censoring model match the true parameter values closely, while those of the naive model are a poor match:

list(
  Censored = select(post_draws, Slope = `bsp_y[1]`, Intercept = Intercept_y, Sigma = sigma_y),
  Naive = select(post_draws_naive, Slope = b_x_star, Intercept = b_Intercept, Sigma = sigma)
) %>% 
  bind_rows(.id = "model") %>% 
  pivot_longer(c(Slope, Intercept, Sigma)) %>% 
  ggplot(aes(x = value, y = name, fill = model)) + 
  scale_fill_manual(values = pal[c(2,3)]) +
  stat_halfeye(slab_alpha = .5, position = "dodge") + 
  labs(x = NULL, y = NULL, fill = NULL)

#To estimate prediction performance for the training data, we first impute the censored and missing values. This generates a list of imputed datasets as long as the number of posterior samples.

data_train_imputed <- impute(data_train, model_rstan, "x_star", sdata$Jmi_xstar, sdata$Jcens_xstar)

#Then, we can calculate RMSE for the training data by generating predictions using the imputed datasets and the posterior draws of the model parameters:

rmse <- function(y, yhat) {
  sqrt(mean((y - yhat) ^ 2))
}
rmse_train <- map2_dbl(
  seq_len(nrow(post_draws)),
  data_train_imputed,
  \(x, y) {
    yhat <- y$x_star * post_draws$`bsp_y[1]`[x] + post_draws$Intercept_y[x]
    rmse(y$y, yhat)
  }
)
median_qi(rmse_train) # summarize

#We can also estimate out-of-sample error by imputing missing values and generating posterior predictions. Imputing censored values in the test data, though, is slightly more complicated: since there is no guarantee that model predictions will fall below the censoring limit, we impute by refitting the model to the training set, augmented by the censored observations from the test set:

data_combined <- bind_rows(data_train, data_test[data_test$cens_x_star, c("x_star", "cens_x_star")])
sdata_combined <- make_standata(bform, data = data_combined) %>% 
  modify_standata(data_combined, lcl, "x_star")

#We fit the same model as we did on the training data. The augmented model excludes all of the response values from the test set, though: these are imputed during model fitting but not used further.

model_rstan_combined <- stan(
  model_code = scode,
  data = sdata_combined,
  sample_file = "model-censored-x-combined",
  seed = stanseed,
  control = list(adapt_delta = .99)
)


#After fitting, we impute the censored values, and fill in the missings by posterior prediction.

data_test_imputed <- impute(
  data = data_test,
  model = model_rstan_combined,
  var = "x_star",
  mi = NULL,
  cens = seq_len(nrow(data_test))[data_test$cens_x_star],
  id = 8:13
) %>% 
  # impute missings via posterior prediction:
  map2(
    seq_len(nrow(post_draws)),
    ~ .x %>% 
        mutate(
          x_star = if_else(
            is.na(x_star), 
            post_draws$Intercept_xstar[.y], 
            x_star
          )
        )
  )

#Then, we generate test predictions... 

preds_test <- map2(
  seq_len(nrow(post_draws)),
  data_test_imputed,
  ~ .y$x_star * post_draws$`bsp_y[1]`[.x] + post_draws$Intercept_y[.x]
)

#... and use them to calculate the RMSE. Since the naive model uses complete cases only, we omit the row with a missing `x` value from the RMSE calculation for the censored model.

rmse_test <- apply(
  do.call(rbind, preds_test)[,-4], 
  1, 
  \(x) rmse(x, data_test$y[-4])
) %>% 
  median_qi()

rmse_naive <- apply(
  posterior_epred(model_naive, newdata = data_test)[ , -4],
  1,
  \(x) rmse(x, data_test$y[-4])
) %>% 
  median_qi()


#Finally, let's compare the predictions:
  
preds_naive <- fitted(model_naive, newdata = data_test)
rmse_labels <- bind_rows(list(cens = rmse_test, naive = rmse_naive), .id = "model") %>% 
  mutate(label = glue("RMSE = {round(y, 2)} (95% CI: {round(ymin, 2)} - {round(ymax, 2)})"))

data_test %>% 
  mutate(
    # predictions from censored model:
    yhat_cens = apply(do.call(rbind, preds_test), 2, mean),
    ymin_cens = apply(do.call(rbind, preds_test), 2, \(x) quantile(x, .025)),
    ymax_cens = apply(do.call(rbind, preds_test), 2, \(x) quantile(x, .975)),
    # predictions from naive model:
    yhat_naive = preds_naive[, "Estimate"],
    ymin_naive = preds_naive[, "Q2.5"],
    ymax_naive = preds_naive[, "Q97.5"]
  ) %>% 
  filter(!is.na(x_star)) %>%  # remove so models are comparable
  pivot_longer(
    matches("y.+"),
    names_to = c(".value", "model"),
    names_pattern = "([^_]+)_([^_]+)"
  ) %>% 
  ggplot(aes(y, yhat)) + 
  facet_wrap(
    vars(model),
    labeller = as_labeller(c(cens = "Censored model", naive = "Naive model"))
  ) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
  geom_abline() +
  geom_point() +
  geom_label(
    data = rmse_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust = "inward", vjust = "inward",
    label.padding = unit(0.35, "lines"),
    label.size = 0
  ) +
  labs(x = "Observations", y = "Predictions")

#And that's it! The linear regression model with censored predictors recovers the true parameter values well and yields better prediction performance than the simple linear regression model.