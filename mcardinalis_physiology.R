library(tidyverse) # for data manipulation
library(lme4) # for mixed models
library(lmtest) # for likelihood ratio tests
library(visreg) # for visualizing effects

mydata_i <- na.omit(data.frame(read.csv("point_measure.csv")))


# Add in covariates 
#These covariates are only for the years 2010-2016, historical data is not included
wna1 <- read_csv("timeseries_lat_2010-2016.csv")
wna2 <- wna1 %>% 
  select(ID_Year1,Latitude,Longitude,Elevation,MAT,MAP,CMD) %>% 
  separate(ID_Year1, into = c("Site", "Year"), sep = "_")
wna2$Site <- as.factor(wna2$Site)
wna2$Year <- as.numeric(wna2$Year)
mydata <- left_join(mydata_i, wna2, by=c("Site", "Year"))

# Scale variables before running models
mydata <- mydata %>% 
  mutate(#Year.scaled = scale(Year), #treat year as category instead?
    Latitude.scaled = scale(Latitude),
    MAT.scaled = scale(MAT),
    MAP.scaled = scale(MAP),
    CMD.scaled = scale(CMD))
# Make sure factors are set correctly
mydata$Year <- as.factor(mydata$Year)
mydata$Site <- as.factor(mydata$Site)
mydata$Block <- as.factor(mydata$Block)
mydata$Plant.ID <- as.factor(mydata$Plant.ID)

## Full models for fixed and random effects 
# General model structure: fixed effects = year*treatment*climate, random effects = family nested within site, block 

# CMD
mod1.cmd= lmer(A ~ Treatment*Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata)
summary(mod1.cmd)

# drop 3-way
mod2 <- lmer(A ~ Treatment*Year + Treatment*CMD.scaled + Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata)
summary(mod1.cmd.no3way)
lrtest(mod1.cmd, mod1.cmd.no3way) # 3-way interaction is better fit

mod3 <- lmer(A ~ Treatment*Year + (1|Site/Plant.ID) + (1|Block), mydata)
mod4 <- lmer(A ~ Treatment*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata)
mod5 <- lmer(A ~ Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(mod1.cmd,mod2,mod3,mod4,mod5) #mod5 best fit, but mod2 also significant 
anova(mod5) #no significant interaction 
anova(mod2) # treatment significant < 2.2e-16 and Treatment*CMD.scaled significant 0.005664

visreg(mod2, xvar="Year", by="Treatment") #Dry treatment significantly higher assimilation rate than wet across plants from years
visreg(mod2, xvar= "CMD.scaled", by = "Treatment") # decreased assimilation rate in plants from dry sites in both treatments, but more pronounced in wet treatment


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
visreg(gsw2, xvar="CMD.scaled", by = "Year") # This sows the shift from wetter climate (2010) to drier (2011-2014) and then back to wetter (2015-2016)






