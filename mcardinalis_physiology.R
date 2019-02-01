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

anova(mod9b.cmd) #no significant interaction 

visreg(mod9b.cmd, xvar="Treatment") #Dry treatment significantly higher assimilation rate than wet across plants from years





#gsw linear models
gsw1 <- lmer(gsw ~ Treatment*Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata)
gsw2 <- lmer(gsw ~ Treatment*Year + Treatment*CMD.scaled + Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata)
gsw3 <- lmer(gsw ~ Treatment*Year + (1|Site/Plant.ID) + (1|Block), mydata)#5.031e-06
gsw4 <- lmer(gsw ~ Treatment*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata) #1.376e-05 
gsw5 <- lmer(gsw ~ Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata)#< 2.2e-16 
lrtest(gsw1,gsw2,gsw3,gsw4,gsw5)
anova(gsw2)#treatment is highly sig < 2.2e-16 in all models 

visreg(gsw2, xvar="Year", by="Treatment")
visreg(gsw2, xvar="CMD.scaled", by="Treatment") #plants from all sites have higher stomatal conductance under the dry treatment than wet, plants from dry sites do worse in the wet treatment
visreg(gsw2, xvar="CMD.scaled", by = "Year") # This shows the shift from wetter climate (2010) to drier (2011-2014) and then back to wetter (2015-2016)


#### Water potential
wpotential <- na.omit(data.frame(read.csv("waterpotentialJuly2018.csv"))) # red sensor 
wpotential
# orient data longwise
wpotential <- wpotential %>%
  gather(Bench, WP, WP1:WP4)

#add AM/PM to time 
wpotential$Time <- paste(wpotential$Time, wpotential$AMPM)

wpotential$AMPM <- NULL
str(wpotential)

# change format to 24 hour
###wpotential2$Time <- format(strptime(wpotential2$Time, "%I:%M %p"), format="%H:%M:%S")

head(wpotential, 20)

## group wet and dry tables

ggplot(data=wpotential, aes(x=Time, y=WP)) +
         geom_point(aes(color=Bench))


#### Mean per day -> This didn't do what I hoped it would - check again

wpotential.mean <- wpotential %>%
  group_by(Measurement, Time, Bench) %>%
  summarise(WP.mean= mean(WP))

head(wpotential.mean, 10)

ggplot(data=wpotential.mean, aes(x=Time, y=WP.mean)) +
  geom_point(aes(color=Bench)) +
  facet_wrap(aes(group=Measurement))


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

