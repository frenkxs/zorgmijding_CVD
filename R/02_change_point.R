# .--------------------------------------------------------------------------------------------
# ----------------------------- Zorgmijding CVD: change point detection------------------------
# .--------------------------------------------------------------------------------------------
# .--------------------------------------------------------------------------------------------
# .--------------------------------------------------------------------------------------------
# 

# This scripts estimates the change point(s) in GP contacts for cardiovascular complaints in the
# period 2017-2020, with a particular emphasis on changes related to the COVID-19 pandemic (from
# March 2020 onwards). It's done separately for each region/database available.
#
# The workflow is as follow:
#
# (1) Determine the ARIMA structure for the time series with the GP contact rate (based on the last
# three years before the pandemic)
# (2) Test for structural change in the whole series, assuming the ARIMA structure determined
# (CUSUM, MUSUM, F test)

# Author: Premysl Velek, p.velek@erasmusmc.nl 




# .-------------------------------------------------------------------------
# Rotterdam ---------------------------------------------------------------
# .-------------------------------------------------------------------------
rm(list = ls())
source(here::here("R", "01_load.R"))
setwd(here::here())


# visualise ------------------------------------------------------
rot <- data$rotterdam$full$n_visits_w_total_40

rot2 <- rot %>%
    dplyr::mutate(week = lubridate::isoweek(date)) %>%
    dplyr::filter(week <= 52) %>%
    dplyr::mutate(year = year(date)) %>%
    mutate(week = yearweek(date))%>%
    as_tsibble(index = week) %>%
    select(week, rate)

autoplot(rot2, rate) +
    labs(title = "GP consultation rate for CVD",
         y = "Number of consultations per 100,000 patients")

# box cox transformation
lambda <- rot2 |> 
    features(rate, features = guerrero) |>
    pull(lambda_guerrero)

# lambda is close to 1 => no transformation is needed

# seasonal difference
rot2 |>
    dplyr::filter(year(week) < 2020) %>%
    gg_tsdisplay(difference(rate, 52),
                 plot_type = 'partial', lag = 52) +
    labs(title = "Seasonally differenced", y = "")



# ARIMA modelling -----------------------------------------------------------------------------

# Get the ARIMA structure
fit_rot <- rot2 |>
    dplyr::filter(year(week) < 2020) |>
    model(
        arima000010 = ARIMA(rate ~ pdq(0,0,0) + PDQ(0,1,0)),
        arima000011d = ARIMA(rate ~ 1 + pdq(0,0,0) + PDQ(0,1,1)),
        arima000011 = ARIMA(rate ~ 0 + pdq(0,0,0) + PDQ(0,1,1)),
        arima000110 = ARIMA(rate ~ pdq(0,0,0) + PDQ(1,1,0)),
        auto = ARIMA(rate, stepwise = FALSE, approx = FALSE)
    )

# print off the models
fit_rot |> pivot_longer(everything(), names_to = "mod_rotel",
                    values_to = "Orders")

# get the best mod_rotels base on AICc
glance(fit_rot) |> arrange(AICc) |> select(.mod_rotel:BIC)


# ARIMA (0,0,0)(1,1,0)[52] w/ drift is the selected model

# get the coefficients and residuals
fit_rot |> select(auto) %>% report
fit_rot |> select(auto) |> gg_tsresiduals(lag = 52)


# Ljung-Box test  
augment(fit_rot) |>
    filter(.mod_rotel == "auto") |>
    features(.innov, ljung_box, lag = 52, dof = 2)



# change point detection --------------------------------------------------


# Convert into OLS regression mod_rotel
# convert to ts object as this is what strucchange expects
rate_rot <- ts(rot2$rate, freq = 52, start = c(2017, 1))
rate_rot <- cbind(rate_rot, 
              stats::lag(rate_rot, k = -52),
              diff(rate_rot, lag = 52)
              )

colnames(rate_rot) <- c("y", "s_lag1", "s_delta_y")
rate_rot <- window(rate_rot, start = c(2017, 1), end = c(2020, 52), frequency = 52)

# model selected above
mod_rot_rot <- s_delta_y ~ s_lag1


# test for structural change(s) 

# estimates based process (fluctuation test) 
re_rate_rot <- efp(mod_rot, data = rate_rot, type = "RE")
plot(re_rate_rot)
plot(re_rate_rot, functional = NULL)

# F-stat test
fs_rate_rot <- Fstats(mod_rot, data = rate_rot, from = 0.1)
plot(fs_rate_rot, main = "supF test")

# CUSUM 
ocus_rate_rot <- efp(mod_rot, data = rate_rot, type = "OLS-CUSUM")
plot(ocus_rate_rot)
rocus_rate_rot <- efp(mod_rot, data = rate_rot, type = "Rec-CUSUM")
plot(rocus_rate_rot)

# MOSUM 
omos_rate_rot <- efp(mod_rot, data = rate_rot, type = "OLS-MOSUM")
plot(omos_rate_rot)
romos_rate_rot <- efp(mod_rot, data = rate_rot, type = "Rec-MOSUM")
plot(romos_rate_rot)


# estimate the date of the change points ----------------------------------

bp_rate_rot <- breakpoints(mod_rot, data = rate_rot, h = 0.1, breaks = 5)
plot(bp_rate_rot)

# plot the results
plot(rate_rot[, "y"])
lines(bp_rate_rot)
lines(confint(bp_rate_rot))



# .-------------------------------------------------------------------------
# Utrecht ---------------------------------------------------------------
# .-------------------------------------------------------------------------
source(here::here("R", "01_load.R"))
setwd(here::here())

# visualise ------------------------------------------------------
rm(list = ls())

utr <- data$utrecht$full$n_visits_w_total_40


utr2 <- utr %>%
    dplyr::mutate(week = lubridate::isoweek(date)) %>%
    dplyr::filter(week <= 52) %>%
    dplyr::mutate(year = year(date)) %>%
    mutate(week = yearweek(date))%>%
    as_tsibble(index = week) %>%
    select(week, rate)

autoplot(utr2, rate) +
    labs(title = "GP consultation rate for CVD",
         y = "Number of consultations per 100,000 patients")


# seasonal difference
utr2 |>
    dplyr::filter(year(week) < 2020) %>%
    gg_tsdisplay(difference(rate, 52),
                 plot_type = 'partial', lag = 52) +
    labs(title = "Seasonally differenced", y = "")

# double differenced
utr2 |>
    dplyr::filter(year(week) < 2020) %>%
    gg_tsdisplay(difference(rate, 52) |> difference(),
                 plot_type = 'partial', lag = 52) +
    labs(title = "Double differenced", y = "")

# ARIMA modelling -----------------------------------------------------------------------------

# Get the ARIMA structure
fit_utr <- utr2 |>
    dplyr::filter(year(week) < 2020) |>
    model(
        arima200010 = ARIMA(rate ~ 1 + pdq(2,0,0) + PDQ(0,1,0)),
        arima002010 = ARIMA(rate ~ 1 + pdq(0,0,2) + PDQ(0,1,0)),
        arima202010 = ARIMA(rate ~ 1 + pdq(2,0,2) + PDQ(0,1,0)),
        arima202011 = ARIMA(rate ~ 1 + pdq(2,0,2) + PDQ(0,1,1)),
        arima202110 = ARIMA(rate ~ 1 + pdq(2,0,2) + PDQ(1,1,0)),
        auto = ARIMA(rate, stepwise = FALSE, approx = FALSE)
    )

# print off the mod_rotels
fit_utr |> pivot_longer(everything(), names_to = "model",
                    values_to = "Orders")

# get the best mod_rotels base on AICc
glance(fit_utr) |> arrange(AICc) |> select(.model:BIC)

# ARIMA (2,0,3)(0,1,1)[52] w/ drift is the best model
# We may also consider ARIMA(2,0,2)(0,1,1)[52] w/ drift as the AICc and BIC are similar


# get the coefficients and residuals
fit_utr |> select(auto) %>% report
fit_utr |> select(arima202011) %>% report

fit_utr |> select(auto) |> gg_tsresiduals(lag = 52) + 
    labs(title = "ARIMA(2,0,3)(0,1,1)[52] w/ drift")
fit_utr |> select(arima202011) |> gg_tsresiduals(lag = 52) + 
    labs(title = "ARIMA(2,0,2)(0,1,1)[52] w/ drift")

# residuals look better in ARIMA(2,0,3)(0,1,1)[52] => this is our model

# Ljung-Box test  
augment(fit_utr) |>
    filter(.model == "auto") |>
    features(.innov, ljung_box, lag = 52, dof = 6)



# change point detection --------------------------------------------------


# Convert into OLS regression mod_rotel
# convert to ts object as this is what strucchange expects
rate_utr <- ts(utr2$rate, freq = 52, start = c(2017, 1))

# estimate the error terms
err <- utr2 |>
    model(errors = ARIMA(rate ~ 1 + pdq(2, 0, 0) + PDQ(0, 1, 0))) |> 
    residuals() |>
    pull(.resid)


err <- ts(err, freq = 52, start = c(2017, 1))

rate_utr <- cbind(
            # y
            rate_utr, 
            
            # p1
            stats::lag(rate_utr, k = -1),
            
            # p2
            stats::lag(rate_utr, k = -2),
            
            # q1
            stats::lag(err, k = -1),
            
            # q2
            stats::lag(err, k = -2),
            
            # q3
            stats::lag(err, k = -3),
            
            # s_p1
            stats::lag(err, k = -52),
            
            # s_d1
            diff(rate_utr, lag = 52)
)

colnames(rate_utr) <- c("y", "lag1", "lag2", "err1", "err2", "err3", "s_err1", "s_delta_y")

rate_utr <- window(rate_utr, start = c(2018, 1), end = c(2020, 52), frequency = 52)

# model selected above
mod_utr <- s_delta_y ~ 1 + lag1 + lag2 + err1 + err2 + err3 + s_err1


# test for structural change(s) ----------------------------------------------


# estimates based process (fluctuation test) 
re_rate_utr <- efp(mod_utr, data = rate_utr, type = "RE")
plot(re_rate_utr)
# plot(re_rate_utr, functional = NULL)

# F-stat test
fs_rate_utr <- Fstats(mod_utr, data = rate_utr, from = 0.1)
plot(fs_rate_utr, main = "supF test")

# CUSUM 
ocus_rate_utr <- efp(mod_utr, data = rate_utr, type = "OLS-CUSUM")
plot(ocus_rate_utr)
rocus_rate_utr <- efp(mod_utr, data = rate_utr, type = "Rec-CUSUM")
plot(rocus_rate_utr)

# MOSUM 
omos_rate_utr <- efp(mod_utr, data = rate_utr, type = "OLS-MOSUM")
plot(omos_rate_utr)
romos_rate_utr <- efp(mod_utr, data = rate_utr, type = "Rec-MOSUM")
plot(romos_rate_utr)


# estimate the date of the change points ----------------------------------

bp_rate_utr <- breakpoints(mod_utr, data = rate_utr, h = 0.1)
summary(bp_rate_utr)
lines(bp_rate_utr, breaks = 2)
bp_rate_utr2 <- breakpoints(bp_rate_utr, breaks = 2)

# plot the results
plot(rate_utr[, "y"])
lines(bp_rate_utr, breaks = 2)
lines(confint(bp_rate_utr, breaks = 2))


# .-------------------------------------------------------------------------
# Amsterdam ---------------------------------------------------------------
# .-------------------------------------------------------------------------
rm(list = ls())
source(here::here("R", "01_load.R"))
setwd(here::here())

# visualise ------------------------------------------------------
ams <- data$amsterdam$full$n_visits_w_total_40


ams2 <- ams %>%
    dplyr::mutate(week = tsibble::yearweek(date),
                  week_n = lubridate::isoweek(date),
                  year = year(date)) %>%
    dplyr::filter(week_n <= 52,
                  year < 2021) %>%
    as_tsibble(index = week) %>%
    select(week, rate)

autoplot(ams2, rate) +
    labs(title = "GP consultation rate for CVD",
         y = "Number of consultations per 100,000 patients")


# seasonal difference
ams2 |>
    dplyr::filter(year(week) < 2020) %>%
    gg_tsdisplay(difference(rate, 52),
                 plot_type = 'partial', lag = 52) +
    labs(title = "Seasonally differenced", y = "")

# double differenced
ams2 |>
    dplyr::filter(year(week) < 2020) %>%
    gg_tsdisplay(difference(rate, 52) |> difference(),
                 plot_type = 'partial', lag = 52) +
    labs(title = "Double differenced", y = "")

# ARIMA modelling -----------------------------------------------------------------------------


# Get the ARIMA structure
fit_ams <- ams2 |>
    dplyr::filter(year(week) < 2020) |>
    model(
        arima200011 = ARIMA(rate ~ 1 + pdq(2,0,0) + PDQ(0,1,1)),
        arima002011 = ARIMA(rate ~ 1 + pdq(0,0,2) + PDQ(0,1,1)),
        arima210011 = ARIMA(rate ~ 1 + pdq(2,1,0) + PDQ(0,1,1)),
        auto = ARIMA(rate, stepwise = FALSE, approx = FALSE)
    )

# print off the mod_rotels
fit_ams |> pivot_longer(everything(), names_to = "model",
                        values_to = "Orders")

# get the best mod_rotels base on AICc
glance(fit_ams) |> arrange(AICc) |> select(.model:BIC)

# ARIMA (2,0,0)(1,1,0)[52] w/ drift is the best model
# We may also consider ARIMA(2,0,0)(0,1,1)[52] w/ drift as the AICc and BIC are similar


# get the coefficients and residuals
fit_ams |> select(auto) %>% report
fit_ams |> select(arima200011) %>% report

fit_ams |> select(auto) |> gg_tsresiduals(lag = 52) + 
    labs(title = "ARIMA(2,0,0)(1,1,0)[52] w/ drift")
fit_ams |> select(arima200011) |> gg_tsresiduals(lag = 52) + 
    labs(title = "ARIMA(2,0,0)(0,1,1)[52] w/ drift")

# residuals look slightly better in ARIMA(2,0,0)(1,1,0)[52] => this is our model

# Ljung-Box test  
augment(fit_ams) |>
    filter(.model == "auto") |>
    features(.innov, ljung_box, lag = 52, dof = 3)

augment(fit_ams) |>
    filter(.model == "arima200011") |>
    features(.innov, ljung_box, lag = 52, dof = 3)




# change point detection --------------------------------------------------


# Convert into OLS regression model
# convert to ts object as this is what strucchange expects
rate_ams <- ts(ams2$rate, frequency = 52, start = c(2017, 1))


rate_ams <- cbind(
    # y
    rate_ams, 
    
    # p1
    stats::lag(rate_ams, k = -1),
    
    # p2
    stats::lag(rate_ams, k = -2),
    
    # s_p1
    stats::lag(rate_ams, k = -52),
    
    # s_d1
    diff(rate_ams, lag = 52)
)

colnames(rate_ams) <- c("y", "lag1", "lag2", "s_lag1", "s_delta_y")

rate_ams <- window(rate_ams, start = c(2018, 1), end = c(2020, 52), frequency = 52)

# model selected above
mod_ams <- s_delta_y ~ 1 + lag1 + lag2 + s_lag1


# test for structural change(s) ----------------------------------------------


# estimates based process (fluctuation test) 
re_rate_ams <- efp(mod_ams, data = rate_ams, type = "RE")
plot(re_rate_ams)
plot(re_rate_ams, functional = NULL)

# F-stat test
fs_rate_ams <- Fstats(mod_ams, data = rate_ams, from = 0.1)
plot(fs_rate_ams, main = "supF test")

# CUSUM 
ocus_rate_ams <- efp(mod_ams, data = rate_ams, type = "OLS-CUSUM")
plot(ocus_rate_ams)
rocus_rate_ams <- efp(mod_ams, data = rate_ams, type = "Rec-CUSUM")
plot(rocus_rate_ams)

# MOSUM 
omos_rate_ams <- efp(mod_ams, data = rate_ams, type = "OLS-MOSUM")
plot(omos_rate_ams)
romos_rate_ams <- efp(mod_ams, data = rate_ams, type = "Rec-MOSUM")
plot(romos_rate_ams)


# estimate the date of the change points ----------------------------------

bp_rate_ams <- breakpoints(mod_ams, data = rate_ams, h = 0.1)
bp_rate_ams

# plot the results
plot(rate_ams[, "y"])
lines(bp_rate_ams)
lines(confint(bp_rate_ams))

