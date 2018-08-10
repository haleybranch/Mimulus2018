# setwd("/Users/haleybranch/Desktop/Branchgithub/Mimulus2018")

mydata <- na.omit(data.frame(read.csv("mimulusjuly2018.csv")))

#Add in yearly weather variables
library(tidyverse)
wna1 <- read_csv("timeseries_lat_2010-2016.csv")
wna2 <- wna1 %>% 
  select(ID_Year1,Latitude,Longitude,Elevation,MAT,MAP,CMD) %>% 
  separate(ID_Year1, into = c("Site", "Year"), sep = "_")
wna2$Site <- as.factor(wna2$Site)
wna2$Year <- as.numeric(wna2$Year)
mydata <- left_join(mydata, wna2, by=c("Site", "Year"))

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

# get variables scaled before running models
mydata.scaled <- mydata %>% 
  mutate(#Year.scaled = scale(Year),
        Latitude.scaled = scale(Latitude),
        MAT.scaled = scale(MAT),
        MAP.scaled = scale(MAP),
        CMD.scaled = scale(CMD))
mydata.scaled$Year <- as.factor(mydata.scaled$Year)

library(lme4) # for mixed models
library(lmtest) # for likelihood ratio tests
library(visreg) # for visualizing effects

## Full model for fixed and random effects 
# fixed effects = year*treatment*climate
# random effects = family nested within site, block 
mod1 = lmer(A ~ Treatment*Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.scaled)
summary(mod1)

# drop 3-way
mod1.no3way <- lmer(A ~ Treatment*Year + Treatment*CMD.scaled + Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.scaled)
summary(mod1.no3way)
lrtest(mod1, mod1.no3way) # 3way significant

visreg(mod1, xvar="Year", by="Treatment")
visreg(mod1, by="Year", xvar="Treatment")
visreg(mod1, xvar="CMD.scaled", by="Treatment", overlay=TRUE) # this is a neat interaction: plants from historically wet sites do better in wet treatment, plants from historically dry sites do better in dry treatment
visreg(mod1, xvar="CMD.scaled", by="Year") # this makes me worry that year effects could be confounded by differences in population sampling and vice versa


## with only 2010 and 2016

#Subset by year
mydata.subYear <- subset(mydata, subset = Year %in% c(2010,2016))

# this time treat year as factor
mydata.subYear.scaled <- mydata.subYear %>% 
  mutate(Latitude.scaled = scale(Latitude),
         MAT.scaled = scale(MAT),
         MAP.scaled = scale(MAP),
         CMD.scaled = scale(CMD))
mydata.subYear.scaled$Year <- as.factor(mydata.subYear.scaled$Year)

## Full model for fixed and random effects 
# fixed effects = year*treatment*climate
# random effects = family nested within site, block 
mod1 = lmer(A ~ Treatment*Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.subYear.scaled)
summary(mod1)

# drop 3-way
mod1.no3way <- lmer(A ~ Treatment*Year + Treatment*CMD.scaled + Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.subYear.scaled)
summary(mod1.no3way)
lrtest(mod1, mod1.no3way) # 3way NS

# drop 2-ways
mod1.noTbyY <- lmer(A ~ Treatment*CMD.scaled + Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.subYear.scaled)
lrtest(mod1.no3way, mod1.noTbyY) # Treatment x Year NS
mod1.noTbyC <- lmer(A ~ Treatment*Year + Year*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.subYear.scaled)
lrtest(mod1.no3way, mod1.noTbyC) # Treatment x CMD NS
mod1.noYbyC <- lmer(A ~ Treatment*Year + Treatment*CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.subYear.scaled)
lrtest(mod1.no3way, mod1.noYbyC) # Year x CMD NS

# drop main effects
mod1.mains <- lmer(A ~ Treatment + CMD.scaled + Year + (1|Site/Plant.ID) + (1|Block), mydata.subYear.scaled)
mod1.noT <- lmer(A ~ CMD.scaled + Year + (1|Site/Plant.ID) + (1|Block), mydata.subYear.scaled)
lrtest(mod1.mains, mod1.noT) #  Treatment significant
mod1.noC <- lmer(A ~ Treatment + Year + (1|Site/Plant.ID) + (1|Block), mydata.subYear.scaled)
lrtest(mod1.mains, mod1.noC) #  CMD not significant
mod1.noY <- lmer(A ~ Treatment + CMD.scaled + (1|Site/Plant.ID) + (1|Block), mydata.subYear.scaled)
lrtest(mod1.mains, mod1.noY) #  Year not significant

visreg(mod1.mains, xvar="Treatment") #lower photo in wet??



