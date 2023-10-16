source(here::here("R", "01_load.R"))


# EDA - visual check ------------------------------------------------------

ams <- data$amsterdam$full$n_visits_w_total_40
utr <- data$utrecht$full$n_visits_w_total_40
rot <- data$rotterdam$full$n_visits_w_total_40

# all three rates plotted together
weekly_total <- bind_rows(list(rotterdam = rot,
                               amsterdam = ams,
                               utrecht = utr), 
                          .id = "Database")

ggplot(weekly_total, aes(x = date, y = rate, colour = Database)) +
    geom_line() +
    geom_point() +
    scale_x_date(limits = c(ymd("2016-12-27"), ymd("2021-01-03")))


# individual plots
ggplot(rot, aes(x = date, y = rate)) +
    geom_line() +
    geom_point() +
    scale_x_date(limits = c(ymd("2016-12-27"), ymd("2021-01-03")))


# add week numbers to it
rot <- rot %>%
    dplyr::mutate(week = lubridate::isoweek(date))


ggplot(rot, aes(x = date, y = rate, label = as.character(week))) +
    geom_line() +
    geom_point() +
    geom_text(nudge_y = - 100) +
    scale_x_date(limits = c(ymd("2016-12-27"), ymd("2021-01-03")))

#' we remove weeks 1 and 52 (and 53 if exist) as they have very low rates and
#' may not be full weeks. This will also give us a periodicity of exactly 50 weeks
rot2 <- rot %>%
    dplyr::filter((week > 1) & (week < 52)) %>%
    dplyr::mutate(year = year(date))

table(rot2$year)


# check visually
ggplot(rot2, aes(x = date, y = rate, label = as.character(week))) +
    geom_line(alpha = 0.5, colour = "grey") +
    geom_point(alpha = 0.5, colour = "grey") +
    geom_text(nudge_y = - 100) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_date(limits = c(ymd("2016-12-27"), ymd("2021-01-03")))




