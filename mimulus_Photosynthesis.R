# setwd("/Users/haleybranch/Desktop/Branchgithub/Mimulus2018")
# don't need to setwd if we are all working within our local copies of the same R project folder

library(tidyverse) # for data manipulation
library(lme4) # for mixed models
library(lmtest) # for likelihood ratio tests
library(visreg) # for visualizing effects

mydata <- na.omit(data.frame(read.csv("mimulusjuly2018.csv")))

# Add in covariates
wna1 <- read_csv("timeseries_lat_2010-2016.csv")
wna2 <- wna1 %>% 
  select(ID_Year1,Latitude,Longitude,Elevation,MAT,MAP,CMD) %>% 
  separate(ID_Year1, into = c("Site", "Year"), sep = "_")
wna2$Site <- as.factor(wna2$Site)
wna2$Year <- as.numeric(wna2$Year)
mydata <- left_join(mydata, wna2, by=c("Site", "Year"))

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
mod1.cmd = lmer(A ~ Treatment*Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata)
summary(mod1.cmd)

# drop 3-way
mod1.cmd.no3way <- lmer(A ~ Treatment*Year + Treatment*CMD.scaled + Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata)
summary(mod1.cmd.no3way)
lrtest(mod1.cmd, mod1.cmd.no3way) # 3-way interaction is highly significant

visreg(mod1.cmd, xvar="Year", by="Treatment")
visreg(mod1.cmd, xvar="Year", by="Treatment", overlay=T)
visreg(mod1.cmd, by="Year", xvar="Treatment")
visreg(mod1.cmd, by="Year", xvar="Treatment", overlay=T) # more variation among years apparent under wet than dry treatment, but no simple temporal progression
visreg(mod1.cmd, xvar="CMD.scaled", by="Treatment", overlay=TRUE) # this is a neat interaction: plants sampled from wet sites x years do better in wet treatment, plants sampled from dry sites x years do better in dry treatment
visreg(mod1.cmd, xvar="CMD.scaled", by="Year") # cline is positive in some early years, switches to more consistently negative in recent years
visreg(mod1.cmd, xvar="CMD.scaled", by="Year", overlay=T) # cline is positive in some early years, switches to more consistently negative in recent years

# MAP
mod1.map = lmer(A ~ Treatment*Year*MAP.scaled + (1|Site/Plant.ID) + (1|Block), mydata)
summary(mod1.map)

# drop 3-way
mod1.map.no3way <- lmer(A ~ Treatment*Year + Treatment*MAP.scaled + Year*MAP.scaled + (1|Site/Plant.ID) + (1|Block), mydata)
summary(mod1.map.no3way)
lrtest(mod1.map, mod1.map.no3way) # 3-way interaction is highly significant

visreg(mod1.map, xvar="Year", by="Treatment")
visreg(mod1.map, xvar="Year", by="Treatment", overlay=T)
visreg(mod1.map, by="Year", xvar="Treatment")
visreg(mod1.map, by="Year", xvar="Treatment", overlay=T) # more variation among years apparent under wet than dry treatment, but no simple temporal progression
visreg(mod1.map, xvar="MAP.scaled", by="Treatment", overlay=TRUE) # this is a neat interaction: plants sampled from wet sites x years do better in wet treatment, plants sampled from dry sites x years do better in dry treatment
visreg(mod1.map, xvar="MAP.scaled", by="Year") # 2013 is poorly sampled, should probably exclude
visreg(mod1.map, xvar="MAP.scaled", by="Year", overlay=T) # heterogeneity among years but no simple temporal progression



### simplify to only 2010 and 2016 (i.e. pre- and post selection)?

# Subset by year
mydata.subYear <- subset(mydata, subset = Year %in% c(2010,2016))

# Rescale
mydata.subYear <- mydata.subYear %>% 
  mutate(Latitude.scaled = scale(Latitude),
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
mod2.cmd = lmer(A ~ Treatment*Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.subYear)
summary(mod2.cmd)

# drop 3-way
mod2.cmd.no3way <- lmer(A ~ Treatment*Year + Treatment*CMD.scaled + Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.subYear)
summary(mod2.cmd.no3way)
lrtest(mod2.cmd, mod2.cmd.no3way) # 3way NS

# drop 2-ways
mod2.cmd.noTbyY <- lmer(A ~ Treatment*CMD.scaled + Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.subYear)
lrtest(mod2.cmd.no3way, mod2.cmd.noTbyY) # Treatment x Year NS
mod2.cmd.noTbyC <- lmer(A ~ Treatment*Year + Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.subYear)
lrtest(mod2.cmd.no3way, mod2.cmd.noTbyC) # Treatment x CMD NS
mod2.cmd.noYbyC <- lmer(A ~ Treatment*Year + Treatment*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.subYear)
lrtest(mod2.cmd.no3way, mod2.cmd.noYbyC) # Year x CMD NS

# drop main effects
mod2.cmd.mains <- lmer(A ~ Treatment + CMD.scaled + Year + (1|Site/Plant.ID) + (1|Block), mydata.subYear)
mod2.cmd.noT <- lmer(A ~ CMD.scaled + Year + (1|Site/Plant.ID) + (1|Block), mydata.subYear)
lrtest(mod2.cmd.mains, mod2.cmd.noT) #  Treatment significant
mod2.cmd.noC <- lmer(A ~ Treatment + Year + (1|Site/Plant.ID) + (1|Block), mydata.subYear)
lrtest(mod2.cmd.mains, mod2.cmd.noC) #  CMD not significant
mod2.cmd.noY <- lmer(A ~ Treatment + CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.subYear)
lrtest(mod2.cmd.mains, mod2.cmd.noY) #  Year not significant

visreg(mod2.cmd.mains, xvar="Treatment") #lower photo in wet?? maybe we need to include measurement day as a covariate?




### haley's initial work
# I want S11, S15, S32, S07 for years 2010 and 2016
#Call names for subsetting
names(mydata)
#Subset by year
mydata.subYear <- subset(mydata, subset = Year %in% c(2010,2016))
#Subset year by site
mydata.finalsub <- subset(mydata.subYear, subset = Site %in% c("S11","S15","S32","S07"))
#print check
mydata.finalsub
mydata.finalsub <- na.omit(data.frame(mydata.finalsub))
str(mydata.finalsub)



mydata.finalsub$Plant.ID <- as.factor(mydata.finalsub$Plant.ID)
mydata.finalsub$Block <- as.factor(mydata.finalsub$Block) 
mydata.finalsub$Site <- as.factor(mydata.finalsub$Site)
mydata.finalsub$Year <- as.factor(mydata.finalsub$Year)


## ggplot

library(tidyverse)
library(Hmisc)
library(lme4)
library(nlme)

ggplot(data=mydata.finalsub, aes(x=Site, y=A)) +
  geom_boxplot(aes(color=Treatment)) + ### potentially add a column that puts site and block together
  facet_wrap(~Year+Block) 

ggplot(data=mydata.finalsub, aes(x=Plant.ID, y=A)) +
  geom_point(aes(color=Treatment)) +
  geom_line(aes(group=Year))


