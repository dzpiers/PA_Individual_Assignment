# Relevant packages
library(ggplot2)
library(lubridate)
library(stringr)
library(forecast)
library(gridExtra)
library(seasonal)
library(tidyverse)
library(TSstudio)
library(stargazer)
library(qpcR)


# Loading in the data
temps <- read.csv("MonthTemp.csv")

# Merging year and month columns
temps$Date <- str_c(temps$Year, '-', temps$Month, '-1')

# Changing year and month to times
temps$Date <- ymd(temps$Date)

# Adding difference variable
temps$Difference <- temps$Highest.temperature - temps$Lowest.temperature

# Making monthly datasets
temps_month <- list()
for (i in seq(12)) {
  temps_month[[i]] <- subset(temps, Month==i)
}

# Making seasonal datasets
temps_summer <- subset(temps, Month %in% c(12, 1, 2))
temps_winter <- subset(temps, Month %in% c(6, 7, 8))


# Plots
# Overall plots
plot1 <- ggplot(data = temps, aes(x=Date)) + 
  geom_line(aes(y=Highest.temperature, color = "Highest temperature"), color = "red") +
  labs(title = "Highest Temperature", y = "Temperature (degrees Celsius)")
plot1

plot2 <- ggplot(data = temps, aes(x=Date)) + 
  geom_line(aes(y=Lowest.temperature, color = "Lowest temperature"), color = "cyan")+
  labs(title = "Lowest Temperature", y = "Temperature (degrees Celsius)")
plot2

plot3 <- ggplot(data = temps, aes(x=Date)) + 
  geom_line(aes(y=Mean.maximum, color = "Mean maximum"), color = "orange")+
  labs(title = "Mean Maximum Temperature", y = "Temperature (degrees Celsius)")
plot3

plot4 <- ggplot(data = temps, aes(x=Date)) + 
  geom_line(aes(y=Mean.minimum, color = "Mean minimum"), color = "blue")+
  labs(title = "Mean Minimum Temperature", y = "Temperature (degrees Celsius)")
plot4

plot5 <- ggplot(data = temps, aes(x=Date)) + 
  geom_line(aes(y=Difference, color = "Difference"), color = "purple")+
  labs(title = "Temperature Difference", y = "Temperature (degrees Celsius)")
plot5

png(filename = "Overall.png", width = 960, height = 480)
grid.arrange(plot1, plot2, plot3, plot4, plot5)
dev.off()

# Monthly plots
plots <- list()
for (i in seq(12)) {
  plots[[i]] <- ggplot(data = temps_month[[i]], aes(x=Date)) + 
    geom_line(aes(y=Highest.temperature, color = "Highest temperature")) +
    geom_line(aes(y=Lowest.temperature, color = "Lowest temperature")) + 
    geom_line(aes(y=Mean.maximum, color = "Mean maximum")) +
    geom_line(aes(y=Mean.minimum, color = "Mean minimum")) +
    labs(title = as.character(month(temps_month[[i]][1,2], label=TRUE)), y = "Temperature") + 
    theme(legend.position="right") +
    scale_color_manual(values = c("red", "cyan", "orange", "blue")) +
    theme(legend.title = element_blank())
}

png(filename = "Monthly.png", width = 960, height = 960)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]], plots[[11]], plots[[12]], ncol = 2)
dev.off()

# Summer plot
plot_summer <- ggplot(data = temps_summer, aes(x=Date)) + 
  geom_line(aes(y=Highest.temperature, color = "Highest.temperature")) +
  geom_line(aes(y=Lowest.temperature, color = "Lowest.temperature")) + 
  geom_line(aes(y=Mean.maximum, color = "Mean.maximum")) +
  geom_line(aes(y=Mean.minimum, color = "Mean.minimum"))
plot_summer

# Winter plot
plot_winter <- ggplot(data = temps_winter, aes(x=Date)) + 
  geom_line(aes(y=Highest.temperature, color = "Highest.temperature")) +
  geom_line(aes(y=Lowest.temperature, color = "Lowest.temperature")) + 
  geom_line(aes(y=Mean.maximum, color = "Mean.maximum")) +
  geom_line(aes(y=Mean.minimum, color = "Mean.minimum"))
plot_winter

# Making a time series object for mean.maximum and difference
temps_ts <- ts(temps$Mean.maximum, frequency = 12, start = c(1970, 7), end = c(2021, 6))
temps_tsd <- ts(temps$Difference, frequency = 12, start = c(1970, 7), end = c(2021, 6))


# Time series decomposition (mean maximum)
# Additive decomposition
fitadd <- decompose(temps_ts, type="additive")
autoplot(fitadd)
fitadd$seasonal
fitadd$trend

#STL with various combinations of window choices
fitstl<-stl(temps_ts,t.window=61, s.window=13)
pd=autoplot(fitstl)
ps=ggsubseriesplot(seasonal(fitstl)) + 
  ggtitle("Seasonal Variations from STL (Global Warming)") + ylab("Seasonal")
png(filename = "GW STL.png", width = 960, height = 480)
grid.arrange(pd,ps,nrow=1)
dev.off()

# Time series decomposition (difference)
# Additive decomposition
fitaddd <- decompose(temps_tsd, type="additive")
autoplot(fitaddd)
fitaddd$seasonal
fitaddd$trend

#STL with various combinations of window choices
fitstld<-stl(temps_tsd,t.window=61, s.window=13)
pdd=autoplot(fitstld)
psd=ggsubseriesplot(seasonal(fitstld)) + 
  ggtitle("Seasonal Variations from STL (Climate Change)") + ylab("Seasonal")
png(filename = "CC STL.png", width = 960, height = 480)
grid.arrange(pdd,psd,nrow=1)
dev.off()


# Regression (mean.maximum)
# Linear trend
fit1 <- tslm(temps_ts~trend+season, data=temps_ts)
summary(fit1)

trend<-1:length(temps_ts)
TT_ts<-trend*fit1$coeff[2] #linear
SS_ts=fit1$fitted.values-TT_ts
par(mfrow=c(3,1))
plot(TT_ts, type="l", main="Trend")
plot(SS_ts, type="l", main="Seasonality")
plot(fit1$resid, main="Residuals")

# Quadratic  trend
fit2 <- tslm(temps_ts~trend+I(trend^2)+season, data=temps_ts)
summary(fit2)

trend<-1:length(temps_ts)
TT_ts2<-trend*fit2$coeff[2]+(trend^2)*fit2$coeff[3] #quadratic
SS_ts2<-fit2$fitted.values-TT_ts2
par(mfrow=c(3,1))
plot(TT_ts2, type="l", main="Trend")
plot(SS_ts2, type="l", main="Seasonality")
plot(fit2$resid, main="Residuals")

# Model selection
CV1 <- rbind(CV(fit1), CV(fit2))
rownames(CV1) <- c("Linear", "Quadratic")
stargazer(t(CV1))


# Forecasts (mean.maximum)
plot(temps_ts)

temps_sets  <- ts_split(temps_ts, sample.out=(2021-2011)*12 + 6)
temps_sets$train
temps_sets$test

fitf1 <- tslm(temps_sets$train~trend+season)
fitf2 <- tslm(temps_sets$train~trend+I(trend^2)+season)

f1 <- forecast(fitf1,h=length(temps_sets$test))
f2 <- forecast(fitf2,h=length(temps_sets$test))
png(filename = "Linear GW.png", width = 960, height = 480)
plot(f1, main = "Linear Forecast (Global Warming)", ylab = "Temperature")
dev.off()
png(filename = "Quadratic GW.png", width = 960, height = 480)
plot(f2, main = "Quadratic Forecast (Global Warming)", ylab = "Temperature")
dev.off()

# Accuracy measures
accuracy(f1, temps_sets$test)
accuracy(f2, temps_sets$test)


# Regression (difference)
# Linear trend
fit1d <- tslm(temps_tsd~trend+season, data=temps_tsd)
summary(fit1d)

trend<-1:length(temps_tsd)
TT_tsd<-trend*fit1d$coeff[2] #linear
SS_tsd=fit1d$fitted.values-TT_tsd
par(mfrow=c(3,1))
plot(TT_tsd, type="l", main="Trend")
plot(SS_tsd, type="l", main="Seasonality")
plot(fit1d$resid, main="Residuals")

# Quadratic  trend
fit2d <- tslm(temps_tsd~trend+I(trend^2)+season, data=temps_tsd)
summary(fit2d)

trend<-1:length(temps_tsd)
TT_ts2d<-trend*fit2d$coeff[2]+(trend^2)*fit2d$coeff[3] #quadratic
SS_ts2d<-fit2d$fitted.values-TT_ts2d
par(mfrow=c(3,1))
plot(TT_ts2d, type="l", main="Trend")
plot(SS_ts2d, type="l", main="Seasonality")
plot(fit2d$resid, main="Residuals")


# Model selection
CV2 <- rbind(CV(fit1d), CV(fit2d))
rownames(CV2) <- c("Linear", "Quadratic")
stargazer(t(CV2))


# Forecasts (difference)
plot(temps_tsd)

temps_setsd  <- ts_split(temps_tsd, sample.out=(2021-2011)*12 + 6)
temps_setsd$train
temps_setsd$test

fitf1d <- tslm(temps_setsd$train~trend+season)
fitf2d <- tslm(temps_setsd$train~trend+I(trend^2)+season)

f1d <- forecast(fitf1d,h=length(temps_setsd$test))
f2d <- forecast(fitf2d,h=length(temps_setsd$test))

png(filename = "Linear CC.png", width = 960, height = 480)
plot(f1d, main = "Linear Forecast (Climate Change)", ylab = "Temperature")
dev.off()
png(filename = "Quadratic CC.png", width = 960, height = 480)
plot(f2d, main = "Quadratic Forecast (Climate Change)", ylab = "Temperature")
dev.off()

accuracy(f1d, temps_setsd$test)
accuracy(f2d, temps_setsd$test)


# Accuracy measures (mean maximum)
f1a <- accuracy(f1, temps_sets$test)[2,]
f2a <- accuracy(f2, temps_sets$test)[2,]
acc_table1 <- rbind(f1a, f2a)
rownames(acc_table1) <- c("Linear", "Quadratic")
stargazer(t(acc_table1))

# DM Test
dm.test((temps_sets$test - f1$mean), (temps_setsd$test - f2$mean), power = 2, alternative = "two.sided")
dm.test((temps_sets$test - f1$mean), (temps_setsd$test - f2$mean), power = 2, alternative = "less")

# Accuracy measures (difference)
f1da <- accuracy(f1d, temps_setsd$test)[2,]
f2da <- accuracy(f2d, temps_setsd$test)[2,]
acc_table2 <- rbind(f1da, f2da)
rownames(acc_table2) <- c("Linear", "Quadratic")
stargazer(t(acc_table2))