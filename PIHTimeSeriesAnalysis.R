# Load packages
library(tidyverse)
library(dynlm)
library(stargazer)
library(quantmod)
library(urca)
library(forecast)
library(AER)

# Clear the environment
rm(list = ls())

# Data Collection ---------------

# Import data
DSPIC96 <- read.csv("DSPIC96.csv")
PCEC96 <- read.csv(("PCEC96.csv"))

# Create dataframe for easier data handling
data <- na.omit(data.frame(date = as.Date(DSPIC96$DATE),
                           income = as.numeric(DSPIC96$DSPIC96),
                           consumpt = as.numeric(PCEC96$PCEC96)))
data$quarterly <- as.yearqtr(data$date)
head(data)
summary(data)

# Data as xts objects, as well as first difference
income <- xts(data$income, data$quarterly)
consumpt <- xts(data$consumpt, data$quarterly)
inc_diff <- 400 * log(income / lag(income)) # Annualised percent
con_diff <- 400 * log(consumpt / lag(consumpt))

# Exploratory Data Analysis ---------------

# Plot data as time series
ggplot(data, aes(x = date)) +
  geom_line(mapping = aes(y = income), color = "blue") +
  geom_line(mapping = aes(y = consumpt), color = "red") +
  labs(title = "Real Disposable Personal Income (Blue)\nReal Personal Consumption Expenditures (Red)",
       x = "Year", y = "Billions of Chained 2017 Dollars")

# Plot log values and first differences
par(mfrow = c(2,2)) # Split plot frame into four

plot(as.zoo(log(consumpt)), xaxt = "n", 
     xlab = "Date",
     ylab = "log(Real Consumption)",
     main = "Real Consumption",
     col = "red",
     lwd = 2)
axis(1, at = seq(2007, 2024.25, 1)) # Remove overcrowded x-axis and then replace

plot(as.zoo(log(income)), xaxt = "n",
     xlab = "Date",
     ylab = "log(Real Income)",
     main = "Real Income",
     col = "blue",
     lwd = 2)
axis(1, at = seq(2007, 2024.25, 1))

plot(as.zoo(con_diff), xaxt = "n", 
     xlab = "Date",
     ylab = "Percentage Change",
     col = "red",
     lwd = 2)
axis(1, at = seq(2007, 2024.25, 1))

plot(as.zoo(inc_diff), xaxt = "n", 
     xlab = "Date",
     ylab = "Percentage Change (Income)",
     col = "blue",
     lwd = 2)
axis(1, at = seq(2007, 2024.25, 1))

# Autocorrelation functions
acf(na.omit(log(consumpt)), lag.max = 24,
    plot = TRUE, main = "Real Consumption")
acf(na.omit(log(income)), lag.max = 24,
    plot = TRUE, main = "Real Income")
acf(na.omit(con_diff), lag.max = 24,
    plot = TRUE, main = "")
acf(na.omit(inc_diff), lag.max = 24,
    plot = TRUE, main = "")

# Transforming to stationarity and order of integration ---------------

#ADF t-stat function for AR models of class 'dynlm'
ADF <- function(model) {
  
  ssr <- sum(model$residuals^2)
  t <- length(model$residuals)
  npar<- length(model$coef)
  tstat <- summary(model)$coef[2, "t value"]
  
  return(
    round(c("p" = npar - 1,
            "AIC" = log(ssr/t) + npar * (2/t),
            "BIC" = log(ssr/t) + npar * log(t)/t,
            "t" = tstat), 4)
  )
}

# Test for unit root and order of integration using ADF function and ur.df
# For log(income)
ADF_income_trend <- sapply(1:8, function(x)
  "UR" = ADF(dynlm(diff(log(ts(income))) ~ L(log(ts(income))) + 
                     trend(log(ts(income)), scale = F) + 
                     diff(L(log(ts(income)),1:x)))))
stargazer(ADF_income_trend, type = "text", title = "Income ADF")
summary(ur.df(log(income),
              type = "trend",
              lags = 3,
              selectlags = "Fixed"))

# For first difference of log(income)
ADF_income_const <- sapply(1:8, function(x)
  "UR" = ADF(dynlm(diff(ts(inc_diff)) ~ L(ts(inc_diff)) +
                     diff(L(ts(inc_diff),1:x)))))
stargazer(ADF_income_const, type = "text", title = "Inc Diff ADF")
summary(ur.df(na.omit(inc_diff),
              type = "drift",
              lags = 2,
              selectlags = "Fixed"))

# For log(consumpt)
ADF_consumpt_trend <- sapply(1:8, function(x)
  "UR" = ADF(dynlm(diff(log(ts(consumpt))) ~ L(log(ts(consumpt))) + 
                     trend(log(ts(consumpt)), scale = F) + 
                     diff(L(log(ts(consumpt)),1:x)))))
stargazer(ADF_consumpt_trend, type = "text", title = "Consumpt ADF")
summary(ur.df(log(consumpt),
              type = "trend",
              lags = 1,
              selectlags = "Fixed"))

# For first difference of log(consumpt)
ADF_consumpt_const <- sapply(1:8, function(x)
  "UR" = ADF(dynlm(diff(ts(con_diff)) ~ L(ts(con_diff)) +
                     diff(L(ts(con_diff),1:x)))))
stargazer(ADF_consumpt_const, type = "text", title = "Con Diff ADF")
summary(ur.df(na.omit(con_diff),
              type = "drift",
              lags = 1,
              selectlags = "Fixed"))

# Lag length selection ---------------

# Information criteria function for AR models of class 'dynlm'
IC <- function(model) {
  
  ssr <- sum(model$residuals^2)
  t <- length(model$residuals)
  npar <- length(model$coef)
  ser <- summary(model)$sigma
  tstat <- abs(summary(model)$coef[length(model$coef), "t value"])
  
  return(
    round(c("p" = npar - 1,
            "SER*" = 2 * log(ser),
            "AIC" = log(ssr/t) + npar * (2/t),
            "BIC" = log(ssr/t) + npar * log(t)/t,
            "|t stat|" = tstat), 4)
  )
}

# Lag length selection using IC function
IC_con <- sapply(1:10, function(x)
  "AR" = IC(dynlm(ts(con_diff) ~ L(ts(con_diff), 1:x))))
IC_inc <- sapply(1:10, function(x)
  "AR" = IC(dynlm(ts(inc_diff) ~ L(ts(inc_diff), 1:x))))
stargazer(IC_con, type = "text", title = "IC Con")
stargazer(IC_inc, type = "text", title = "IC Inc")

IC_con[, which.min(IC_con[2, ])] # min SER*
IC_con[, which.min(IC_con[3, ])] # min AIC
IC_con[, which.min(IC_con[4, ])] # min BIC

IC_inc[, which.min(IC_inc[2, ])]
IC_inc[, which.min(IC_inc[3, ])]
IC_inc[, which.min(IC_inc[4, ])]

con_intercept <-dynlm(ts(con_diff) ~ 1)
con_ar1 <- dynlm(ts(con_diff) ~ L(ts(con_diff)))
inc_ar3 <- dynlm(ts(inc_diff) ~ L(ts(inc_diff),1:3))

# ADL model selection using IC function
IC_con_adl <- sapply(1:10, function(x)
  "ADL" = IC(dynlm(ts(con_diff) ~ L(ts(con_diff),1:x) + L(ts(inc_diff),1:x))))
stargazer(IC_con_adl, type = "text", title = "IC Con ADL")
con_adl <- dynlm(ts(con_diff) ~ L(ts(con_diff)) + L(ts(inc_diff)))
summary(con_adl)

IC_inc_adl <- sapply(1:10, function(x)
  "ADL" = IC(dynlm(ts(inc_diff) ~ L(ts(inc_diff),1:x) + L(ts(con_diff),1:x))))
stargazer(IC_inc_adl, type = "text", title = "IC Inc ADL")
inc_adl <- dynlm(ts(inc_diff) ~ L(ts(inc_diff),1:3) + L(ts(con_diff),1:3))
summary(inc_adl)

# Granger Causality ---------------

# Granger Causality F test function
Granger <- function(model) {
  
  npar <- length(model$coefficients)-1
  startindex <- (length(model$coef) + 1) / 2
  coefnames <- names(model$coef)[-(1:startindex)] # Parameter names excluded from restricted model
  Fstat <- linearHypothesis(model,
                            c(coefnames),
                            vcov. = sandwich)$F[2]
  pval <- linearHypothesis(model,
                           c(coefnames),
                           vcov. = sandwich)$`Pr(>F)`[2]
  
  return(round(
    c("p" = npar/2,
      "F" = Fstat,
      "p-value" = pval),4)
  )
}

# F test for income Granger-causing consumption
granger_con <- sapply(1:10, function(x)
  "GT" = Granger(dynlm(ts(con_diff) ~ L(ts(con_diff), 1:x) + L(ts(inc_diff), 1:x))))
stargazer(granger_con, type = "text", title = "Granger Con")

# F test for consumption Granger-causing income
granger_inc <- sapply(1:10, function(x)
  "GT" = Granger(dynlm(ts(inc_diff) ~ L(ts(inc_diff), 1:x) + L(ts(con_diff), 1:x))))
stargazer(granger_inc, type = "text", title = "Granger Inc")

# Structural Breaks ---------------

# Find start and end index for a 15% trimming parameter
0.15 * length(con_diff) # 10.5
0.85 * length(con_diff) # 59.5

con_diff[11] # 2009 Q3
con_diff[60] # 2021 Q4

# Define tau for 2009 Q3 to 2021 Q4
tau <- seq(2009.50, 2021.75, 0.25)

# Initialise vector of F-stats for QLR tests
Fstats_break_con_int <- numeric(length(tau))
Fstats_break_con_ar1 <- numeric(length(tau))
Fstats_break_inc <- numeric(length(tau))

# QLR test for consumption using both AR1 and intercept-only models
for (i in 1:length(tau)) {
  
  D <- time(con_diff) > tau[i]
  
  breaktest <- dynlm(ts(con_diff) ~ 1 + D*1)
  
  Fstats_break_con_int[i] <- linearHypothesis(breaktest,
                                              c("DTRUE=0"),
                                              vcov. = sandwich)$F[2]
}

for (i in 1:length(tau)) {
  
  D <- time(con_diff) > tau[i]
  
  breaktest <- dynlm(ts(con_diff) ~ L(ts(con_diff)) + D*L(ts(con_diff)))
  
  Fstats_break_con_ar1[i] <- linearHypothesis(breaktest,
                                              c("DTRUE=0",
                                                "L(ts(con_diff)):DTRUE"),
                                              vcov. = sandwich)$F[2]
}

QLR_con_int <- max(Fstats_break_con_int) # QLR statistic
QLR_con_int
as.yearqtr(tau[which.max(Fstats_break_con_int)]) # Corresponding break date estimator

QLR_con_ar1 <- max(Fstats_break_con_ar1)
QLR_con_ar1
as.yearqtr(tau[which.max(Fstats_break_con_ar1)])

# Plot QLR statistics for consumption
par(mfrow = c(1,1)) # Reset plot frame

series_break_con_int <- ts(Fstats_break_con_int, 
                           start = tau[1],
                           end = tau[length(tau)],
                           frequency = 4)
plot(series_break_con_int, 
     xlim = c(2007, 2024), 
     ylim = c(0, 15), 
     col = "blue", 
     ylab = "")
abline(h = 8.68, lty = 2, col = "orange")
abline(h = 12.16, lty = 2, col = "red")
legend("topleft", 
       lty = c(2, 2),
       col = c("red", "orange"),
       legend = c("1% Critical Value", "5% Critical Value"))

series_break_con_ar1 <- ts(Fstats_break_con_ar1,
                           start = tau[1],
                           end = tau[length(tau)],
                           frequency = 4)
plot(series_break_con_ar1,
     xlim = c(2007, 2024),
     ylim = c(0, 15),
     col = "blue",
     ylab = "")
abline(h = 5.86, lty = 2, col = "orange")
abline(h = 7.78, lty = 2, col = "red")
legend("topleft", 
       lty = c(2, 2),
       col = c("red", "orange"),
       legend = c("1% Critical Value", "5% Critical Value"))

# QLR test for income AR3 model
for (i in 1:length(tau)) {
  
  D <- time(inc_diff) > tau[i]
  
  breaktest <- dynlm(ts(inc_diff) ~ L(ts(inc_diff)) + L(ts(inc_diff), 2) + L(ts(inc_diff), 3) + 
                       D*L(ts(inc_diff)) + D*L(ts(inc_diff), 2) + D*L(ts(inc_diff), 3))
  
  Fstats_break_inc[i] <- linearHypothesis(breaktest,
                                          c("DTRUE=0",
                                            "L(ts(inc_diff)):DTRUE",
                                            "L(ts(inc_diff), 2):DTRUE",
                                            "L(ts(inc_diff), 3):DTRUE"),
                                          vcov. = sandwich)$F[2]
}

QLR_inc <- max(Fstats_break_inc)
QLR_inc
as.yearqtr(tau[which.max(Fstats_break_inc)]) # 2021 Q2

series_break_inc <- ts(Fstats_break_inc,
                       start = tau[1],
                       end = tau[length(tau)],
                       frequency = 4)
plot(series_break_inc,
     xlim = c(2007, 2024),
     col = "blue",
     ylab = "")
abline(h = 4.09, lty = 2, col = "orange")
abline(h = 5.12, lty = 2, col = "red")
legend("topleft", 
       lty = c(2, 2),
       col = c("red", "orange"),
       legend = c("1% Critical Value", "5% Critical Value"))

# Rerun QLR test on period before 2021 Q2
0.15 * length(inc_diff["2007/2021-1"]) # 8.55
0.85 * length(inc_diff["2007/2021-1"]) # 48.45

inc_diff[9] # 2009 Q1
inc_diff[48] # 2018 Q4

tau2 <- seq(2009.00, 2018.75, 0.25)

Fstats_break_inc2 <- numeric(length(tau2))

for (i in 1:length(tau2)) {
  
  D <- time(inc_diff) > tau2[i]
  
  breaktest <- dynlm(ts(inc_diff) ~ L(ts(inc_diff)) + L(ts(inc_diff), 2) + L(ts(inc_diff), 3) + 
                       D*L(ts(inc_diff)) + D*L(ts(inc_diff), 2) + D*L(ts(inc_diff), 3))
  
  Fstats_break_inc2[i] <- linearHypothesis(breaktest,
                                           c("DTRUE=0",
                                             "L(ts(inc_diff)):DTRUE",
                                             "L(ts(inc_diff), 2):DTRUE",
                                             "L(ts(inc_diff), 3):DTRUE"),
                                           vcov. = sandwich)$F[2]
}

QLR_inc2 <- max(Fstats_break_inc2)
QLR_inc2
as.yearqtr(tau2[which.max(Fstats_break_inc2)]) # 2009 Q1

series_break_inc2 <- ts(Fstats_break_inc2,
                        start = tau2[1],
                        end = tau2[length(tau2)],
                        frequency = 4)
plot(series_break_inc2,
     xlim = c(2007, 2024),
     col = "blue",
     ylab = "")
abline(h = 4.09, lty = 2, col = "orange")
abline(h = 5.12, lty = 2, col = "red")
legend("topright", 
       lty = c(2, 2),
       col = c("red", "orange"),
       legend = c("1% Critical Value", "5% Critical Value"))

# Rerun QLR excluding data before 2009 Q2
income_sub <- inc_diff["2009-3/2021-1"]
0.15 * length(income_sub) # 7.2
0.85 * length(income_sub) # 40.8

income_sub[7] # 2011 Q2
income_sub[41] # 2019 Q3

tau3 <- seq(2011.25, 2019.50, 0.25)

Fstats_break_inc3 <- numeric(length(tau3))

for (i in 1:length(tau3)) {
  
  D <- time(inc_diff) > tau3[i]
  
  breaktest <- dynlm(ts(inc_diff) ~ L(ts(inc_diff)) + L(ts(inc_diff), 2) + L(ts(inc_diff), 3) + 
                       D*L(ts(inc_diff)) + D*L(ts(inc_diff), 2) + D*L(ts(inc_diff), 3))
  
  Fstats_break_inc3[i] <- linearHypothesis(breaktest,
                                           c("DTRUE=0",
                                             "L(ts(inc_diff)):DTRUE",
                                             "L(ts(inc_diff), 2):DTRUE",
                                             "L(ts(inc_diff), 3):DTRUE"),
                                           vcov. = sandwich)$F[2]
}

QLR_inc3 <- max(Fstats_break_inc3)
QLR_inc3
as.yearqtr(tau3[which.max(Fstats_break_inc3)]) # 2010 Q4

series_break_inc3 <- ts(Fstats_break_inc3,
                        start = tau3[1],
                        end = tau3[length(tau3)],
                        frequency = 4)
plot(series_break_inc3,
     xlim = c(2007, 2024),
     ylim = c(0, 8),
     col = "blue",
     ylab = "")
abline(h = 4.09, lty = 2, col = "orange")
abline(h = 5.12, lty = 2, col = "red")
legend("topright", 
       lty = c(2, 2),
       col = c("red", "orange"),
       legend = c("1% Critical Value", "5% Critical Value"))

# Forecasting ---------------

# Re-check lag length using stable epoch income_sub
IC_inc <- sapply(1:10, function(x)
  "AR" = IC(dynlm(ts(income_sub) ~ L(ts(income_sub), 1:x))))
stargazer(IC_inc, type = "text", title = "Inc AR Check")

IC_inc_adl <- sapply(1:10, function(x)
  "ADL" = IC(dynlm(ts(income_sub) ~ L(ts(income_sub),1:x) + L(ts(con_diff),1:x))))
stargazer(IC_inc_adl, type = "text", title = "Inc ADL Check")

# Re-check whether consumption Granger-causes income
granger_inc <- sapply(1:10, function(x)
  "GT" = Granger(dynlm(ts(income_sub) ~ L(ts(income_sub),1:x) + L(ts(con_diff),1:x))))
stargazer(granger_inc, type = "text", title = "Inc Granger Check")

# Model using dynlm
inc_ar3 <- dynlm(ts(income_sub) ~ L(ts(income_sub),1:3))
inc_adl <- dynlm(ts(income_sub) ~ L(ts(income_sub),1:3) + L(ts(con_diff),1:3))

names(inc_ar3$coefficients) <- c("(Intercept)", 
                                 "Income Lag 1", 
                                 "Income Lag 2", 
                                 "Income Lag 3")

names(inc_adl$coefficients) <- c("(Intercept)", 
                                 "Income Lag 1", 
                                 "Income Lag 2", 
                                 "Income Lag 3", 
                                 "Consumpt Lag 1", 
                                 "Consumpt Lag 2", 
                                 "Consumpt Lag 3")

stargazer(inc_ar3, inc_adl, type = "text")

# Check residuals
acf(inc_ar3$residuals, lag.max = 24, plot = TRUE, main = "")
acf(inc_adl$residuals, lag.max = 24, plot = TRUE, main = "")

# Alternative modelling using ar.ols
ar.ols(income_sub,
       order.max = 10,
       aic = TRUE,
       demean = FALSE,
       intercept = TRUE)

# Using arima
arima(income_sub, order = c(3, 0, 0), method = "CSS") # Conditional sum-of-squares so model is same

# Compute forecasts and prediction intervals for the next 12 periods
inc_arima <- arima(income_sub["2009/2020"], order = c(3, 0, 0), method = "CSS")
fc <- forecast(inc_arima, h = 12, level = seq(5, 99, 10))

# Plot fan chart and compare
par(mfrow = c(2,1))

plot(fc, 
     main = "Forecast Fan Chart for AR(3) Model of Simulated Data",
     showgap = F, 
     fcol = "red",
     flty = 2)

plot(as.zoo(inc_diff["2009-3/2024-1"]),
     ylab = "",
     xlab = "Date")

# Cointegration ---------------

# Estimate cointegrating relationship
coin_model <- lm(log(consumpt) ~ log(income))
summary(coin_model) # int = 1.0950, theta-hat = 0.8748

# Test for stationarity in residuals
coin_res <- ts(coin_model$residuals, start = 2007, end = 2024.25, frequency = 4)

par(mfrow = c(1,1))
plot(coin_res,
     ylab = "")

ur.df(coin_res, 
      type = "drift",
      lags = 6,
      selectlags = "BIC")


# Re-estimate and test up to 2020
coin_model2 <- lm(log(consumpt["2007/2020"]) ~ log(income["2007/2020"]))
summary(coin_model2)

coin_res2 <- ts(coin_model2$residuals, start = 2007, end = 2020.75, frequency = 4)

par(mfrow = c(1,1))
plot(coin_res2,
     ylab = "")

ur.df(coin_res2, 
      type = "drift",
      lags = 6,
      selectlags = "BIC")