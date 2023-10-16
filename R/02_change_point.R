source(here::here("R", "01_load.R"))


# EDA - visual check ------------------------------------------------------



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

rot2 |>
    dplyr::filter(year(week) < 2020) %>%
    gg_tsdisplay(difference(rate, 52),
                 plot_type = 'partial') +
    labs(title = "Seasonally differenced", y = "")



# Remove the outlying data in the first and last weeks of the year
mu <- mean(rot$rate)
rot3 <- rot |>
    dplyr::mutate(week = lubridate::isoweek(date)) %>%
    dplyr::mutate(rate = case_when(week == 1 ~  mu,
                                   week == 52 ~ mu,
                                   week == 53 ~ mu,
                                   .default = rate)) %>%
    dplyr::mutate(year = year(date)) %>%
    mutate(week = yearweek(date))%>%
    as_tsibble(index = week) %>%
    select(week, rate)

rot3 |> 
    dplyr::filter(year(week) < 2020) %>%
    autoplot(rate) +
    labs(title = "GP consultation rate for CVD",
         y = "Number of consultations per 100,000 patients")
    
rot3 |>
    dplyr::filter(year(week) < 2020) %>%
    gg_tsdisplay(difference(rate, 50),
                 plot_type = 'partial') +
    labs(title = "Seasonally differenced", y = "")

rot3 |>
    dplyr::filter(year(week) < 2020) %>%
    pull(rate) |> 
    difference(lag = 50) |> 
    ggAcf()
  

fit <- rot3 |>
    dplyr::filter(year(week) < 2020) |>
    model(
        arima000010 = ARIMA(rate ~ pdq(0,0,0) + PDQ(0,1,0)),
        arima100010 = ARIMA(rate ~ pdq(1,0,0) + PDQ(0,1,0)),
        arima001010 = ARIMA(rate ~ pdq(0,0,1) + PDQ(0,1,0)),
        arima101010 = ARIMA(rate ~ pdq(1,0,1) + PDQ(0,1,0)),
        auto = ARIMA(rate, stepwise = FALSE, approx = FALSE)
    )

fit |> pivot_longer(everything(), names_to = "Model name",
                    values_to = "Orders")

glance(fit) |> arrange(AICc) |> select(.model:BIC)
fit |> select(auto) %>% report
fit |> select(auto) |> gg_tsresiduals(lag = 52)

augment(fit) |>
    filter(.model == "auto") |>
    features(.innov, ljung_box, lag = 52, dof = 1)

# NOw convert into OLS regression model
# convert to ts object as this is what strucchange expects
rate <- ts(rot3$rate, freq = 52, start = c(2017, 1))

rate <- cbind(rate, 
              stats::lag(rate, k = -52),
              diff(rate, lag = 52)
              )
colnames(rate) <- c("y", "s_lag1", "delta_y")


plot(consults)
mod <- delta_y ~ s_lag1


# estimates based process (fluctuation test) 
re_rate <- efp(mod, data = rate, type = "RE")
plot(re_rate)
plot(re.rate, functional = NULL)

# F-stat test
fs_rate <- Fstats(mod, data = consults, from = 0.1)
plot(fs_rate, main = "supF test")

# CUSUM 
ocus_rate <- efp(mod, data = consults, type = "OLS-CUSUM")
plot(ocus_rate)
rocus_rate <- efp(mod, data = consults, type = "Rec-CUSUM")
plot(rocus_rate)

# MOSUM 
omos_rate <- efp(mod, data = consults, type = "OLS-MOSUM")
plot(omos_rate)
romos_rate <- efp(mod, data = consults, type = "Rec-MOSUM")
plot(romos_rate)

# dating the point change
bp_rate <- breakpoints(mod, data = consults, h = 0.1, breaks = 5)
bp_rate
plot(bp_rate)

# plot the results
plot(consults[, "y"])
lines(bp_rate)
lines(confint(bp_rate))