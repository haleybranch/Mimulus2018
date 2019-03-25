library(tidyverse) # for data manipulation
library(lme4) # for mixed models
library(lmtest) # for likelihood ratio tests
library(visreg) # for visualizing effects

mydata_i <- na.omit(data.frame(read.csv("point_measure.csv")))


# Add in covariates #
# Long-term climate for hypothesis about propensity to evolve or be plastic
wna <- read_csv("timeseries_lat_Normal_1981_2010Y.csv") %>% 
  select(Site=ID, MAT.clim=MAT,MAP.clim=MAP,CMD.clim=CMD)
wna$Site <- as.factor(wna$Site)

#These covariates are only for the years 2010-2016, use this to calculate anomalies
wna1 <- read_csv("timeseries_lat_2010-2016.csv")
wna2 <- wna1 %>% 
  select(ID_Year1,Latitude,Longitude,Elevation,MAT.weath=MAT,MAP.weath=MAP,CMD.weath=CMD) %>% 
  separate(ID_Year1, into = c("Site", "Year"), sep = "_")
wna2$Site <- as.factor(wna2$Site)
wna2$Year <- as.numeric(wna2$Year)

# join climate and weather 
wna_all <- left_join(wna2, wna, by="Site") %>% 
  mutate(CMD.anom = CMD.clim-CMD.weath,
         MAT.anom = MAT.clim-MAT.weath,
         MAP.anom = log(MAP.clim)-log(MAP.weath))

#join them onto response variables
mydata <- left_join(mydata_i, wna_all, by=c("Site", "Year"))

# Scale variables before running models
mydata <- mydata %>% 
  mutate(#Year.scaled = scale(Year), #treat year as category instead?
    Latitude.scaled = scale(Latitude),
    MAT.clim.scaled = scale(MAT.clim),
    MAP.clim.scaled = scale(MAP.clim),
    CMD.clim.scaled = scale(CMD.clim),    
    MAT.weath.scaled = scale(MAT.weath),
    MAP.weath.scaled = scale(MAP.weath),
    CMD.weath.scaled = scale(CMD.weath),    
    MAT.anom.scaled = scale(MAT.anom),
    MAP.anom.scaled = scale(MAP.anom),
    CMD.anom.scaled = scale(CMD.anom))

plot(CMD.anom.scaled~CMD.clim.scaled, data=mydata)
plot(CMD.weath.scaled~CMD.clim.scaled, data=mydata)

# Make sure factors are set correctly
mydata$Year <- as.factor(mydata$Year)
mydata$Site <- as.factor(mydata$Site)
mydata$Block <- as.factor(mydata$Block)
mydata$Plant.ID <- as.factor(mydata$Plant.ID)

## Full models for fixed and random effects 
# General model structure: fixed effects = treatment*climate*anomaly, random effects = year, family nested within site, block 

# CMD
mod1.cmd= lmer(A ~ Treatment*CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
summary(mod1.cmd)
anova(mod1.cmd)

# drop 3-way
mod2.cmd <- lmer(A ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
summary(mod2.cmd)
lrtest(mod1.cmd, mod2.cmd) # this model is not better or worse than more complicated mod1; simplify to this one

# Drop 2-ways singly
##Drop Trt*anomaly
mod3.cmd <- lmer(A ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(mod2.cmd,mod3.cmd) #this model is not better or worse than more complicated mod2, no support for retaining Trt*anom
## Drop Trt*climate
mod4.cmd <- lmer(A ~ Treatment*CMD.anom.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(mod2.cmd,mod4.cmd) # mod4 is significantly better than mod2, definitely drop Trt*clim
mod4b.cmd <- lmer(A ~ CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(mod3.cmd,mod4b.cmd) # mod4b is significantly better than mod3, definitely drop Trt*clim

## Drop anom*climate
mod5.cmd <- lmer(A ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(mod2.cmd,mod5.cmd) #mod5 is no better or worse than mod2, no support for retaining clim*anom

# Go down to main effects only
mod6.cmd <- lmer(A ~ Treatment + CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
# Drop Treatment
mod7.cmd <- lmer(A ~  CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(mod7.cmd, mod6.cmd) #Treatment is significant; retain it
# Drop Climate
mod8.cmd <- lmer(A ~ Treatment + CMD.anom.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(mod8.cmd, mod6.cmd) #Climate is not better or worse than simpler model, so drop it
# Drop Anomaly
mod9.cmd <- lmer(A ~ Treatment + CMD.clim.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
mod9b.cmd <- lmer(A ~ Treatment + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(mod9.cmd, mod6.cmd) #Anomaly is not better or worse than simpler model, so drop it
lrtest(mod9b.cmd, mod8.cmd) #Anomaly is not better or worse than simpler model, so drop it

visreg(mod9b.cmd, xvar="Treatment") #Dry treatment significantly higher assimilation rate than wet across plants from years

# p value - Scatterthwaite's approx p value and CI
library(lmerTest)
# CIs <- confint(model, method = "profile") # gives confidence intervals 
CIs <- confint(mod9b.cmd, method = "profile")
summary(mod9b.cmd) ## pvalue is <2e-16 for mod9b.cmd --> treatment significant



#gsw linear models

gsw1.cmd= lmer(gsw ~ Treatment*CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
summary(gsw1.cmd)
anova(gsw1.cmd)

# drop 3-way
gsw2.cmd <- lmer(gsw ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
summary(gsw2.cmd)
lrtest(gsw1.cmd, gsw2.cmd) # this model is better than more complicated mod1; simplify to this one

# Drop 2-ways singly
##Drop Trt*anomaly
gsw3.cmd <- lmer(gsw ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(gsw2.cmd,gsw3.cmd) #this model is better than more complicated gsw2, remove Trt*anom
## Drop Trt*climate
gsw4.cmd <- lmer(gsw ~ Treatment*CMD.anom.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(gsw2.cmd,gsw4.cmd) # gsw4 is significantly better than gsw2, drop Trt*clim
gsw4b.cmd <- lmer(gsw ~ CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(gsw3.cmd,gsw4b.cmd) # gsw4b is significantly better than gsw3, definitely drop Trt*clim

## Drop anom*climate
gsw5.cmd <- lmer(gsw ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(gsw2.cmd,gsw5.cmd) #gsw5 is better than gsw2, drop clim*anom

# Go down to main effects only
gsw6.cmd <- lmer(gsw ~ Treatment + CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
# Drop Treatment
gsw7.cmd <- lmer(gsw ~  CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(gsw7.cmd, gsw6.cmd) #Treatment is significant; retain it
# Drop Climate
gsw8.cmd <- lmer(gsw ~ Treatment + CMD.anom.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(gsw8.cmd, gsw6.cmd) #Climate is not better than simpler model, keep it
# Drop Anomaly
gsw9.cmd <- lmer(gsw ~ Treatment + CMD.clim.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(mod9.cmd, mod6.cmd) #Anomaly is not better or worse than simpler model, so drop it
gsw9b.cmd <- lmer(gsw ~ Treatment + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(gsw9b.cmd, gsw9.cmd) #Anomaly is not better or worse than simpler model, so drop it

visreg(gsw9.cmd, xvar="Treatment") #Dry treatment significantly higher stomatal conductance rate than wet across plants from years

# p value - Scatterthwaite's approx p value and CI
library(lmerTest)
# CIs <- confint(model, method = "profile") # gives confidence intervals 
CIs <- confint(gsw9.cmd, method = "profile")
summary(gsw9.cmd) ## pvalue is <2e-16 for gsw9.cmd --> treatment significant


#### Water potential
wpotential <- na.omit(data.frame(read.csv("waterpotentialJuly2018.csv"))) # red sensor 
head(wpotential)
tail(wpotential)

# orient data longwise
wpotential <- wpotential %>%
  gather(Bench, WP, WP1:WP4)
head(wpotential)

# add AM/PM to time 
wpotential <- wpotential %>% 
  mutate(Time2 = paste(Time, AMPM),
         DateTime = paste(Measurement, Time2))
head(wpotential)

# change format to 24 hour, convert format for ggplot
wpotential$Time3 <- as.POSIXct(strptime(wpotential$DateTime, "%Y-%m-%d %I:%M:%S %p"))
wpotential$Day <- as.POSIXct(wpotential$Measurement)
head(wpotential, 20)
str(wpotential)

# remove old date and time columns, 
wpotential <- wpotential %>% 
  select(Day, DayTime = Time3, Bench, WP)
head(wpotential)
tail(wpotential)

# graph of raw WP by bench
ggplot(data=wpotential, aes(x=DayTime, y=WP)) +
  geom_point(aes(color=Bench)) 

# mean per day 
wpotential.mean <- wpotential %>%
  group_by(Bench, Day) %>%
  summarise(WP.mean = mean(WP)) %>% 
  ungroup() %>% 
  select(Bench, Day, WP.mean)
head(wpotential.mean, 10)

# graph of mean wp per day
ggplot(data=wpotential.mean, aes(x=Day, y=WP.mean)) +
  geom_point(aes(color=Bench))

## amy stopped here


##### creating new subsets for days of interest
#attempt 1
wpotential.start <- wpotential.mean %>%
  group_by(Measurement)
wpotential.start

#attempt 2
wpotential.mean$Measurement <- as.numeric(wpotential.mean$Measurement)
drought.start <- wpotential.mean[wpotential$Measurement == c(2018-07-01, 2018-07-14),]
drought.start <- data.frame(Measurement = c(2018-07-01, 2018-07-14))

#attempt3
plot(wpotential.mean$Time[wpotential.mean$Measurement < 2018-07-14],
     wpotential.mean$WP[wpotential.mean$Measurement < 2018-07-14])

#attempt 4
drought.start <- subset(wpotential2.mean, wpotential2.mean$Measurement > 2018-07-01 & wpotential2.mean$Measurement < 2018-07-15)

#attempt 5
replace.value(wpotential.mean, Measurement, from = 2018-07-01, to=as.factor(1), verbose = FALSE)
wpotential.mean

#attempt 6
wpotential.mean$Measurement <- as.character(wpotential.mean$Measurement)

#attempt 7
# create a new column using mutate with an if else statement 
# or split the measurement column into year, month, day -- and use if the day is less than 14 ...
# splitstring function?

#step 1 - separate YYYY-MM-DD into columns
wpotential.mean <- separate(wpotential.mean, "Measurement", c("Year", "Month", "Day"), sep = "-")

#step 2 - create new column with 3 stages: 1 = 3/4 water, 2 = 1/2 water, 3 = 0 water
wpot <- as.data.frame(wpotential.mean)
wpot$Day <- as.numeric(wpot$Day)

wpot$Stage[wpot$Day<=14]<-1
wpot$Stage[wpot$Day>14 & wpot$population<30] <-2
wpot$Stage[wpot$Day>=30] <-3
wpot$Stage <- as.character(wpot$Stage)

#graph the daily water potential during each stage 
ggplot(data=wpot, aes(x=Time, y=WP.mean)) +
  geom_point(aes(color=Bench)) +
  facet_wrap(aes(group=Stage))

