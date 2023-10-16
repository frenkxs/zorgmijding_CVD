library(zorgmijding)

library(tidyverse)
library(colorspace)
library(strucchange)
library(forecast)
library(fpp3)
library(fable)

here::here()


# Import data -------------------------------------------------------------

data <- new.env()


local({
    amsterdam <- new.env(parent = data)
    rotterdam <- new.env(parent = data)
    utrecht <- new.env(parent = data)
    
    # load amsterdam data
    local({
        load(here::here("umc_data", "amsterdam", "contact_rates_long.RData"),
             full <- new.env(parent = amsterdam))
        load(here::here("umc_data", "amsterdam", "contact_rates.RData"),
             averages <- new.env(parent = amsterdam))
    }, envir = amsterdam)

    # utrecht data
    local({
        load(here::here("umc_data", "utrecht", "contact_rates_long.RData"),
             full <- new.env(parent = utrecht))
        load(here::here("umc_data", "utrecht", "contact_rates.RData"),
             averages <- new.env(parent = utrecht))
    }, envir  = utrecht)
    
    # rotterdam data
    local({
        load(here::here("umc_data", "rotterdam", "contact_rates_long.RData"),
             full <- new.env(parent = rotterdam))
        load(here::here("umc_data", "rotterdam", "contact_rates.RData"),
             averages <- new.env(parent = rotterdam))
    }, envir = rotterdam)

}, envir = data)



    
    