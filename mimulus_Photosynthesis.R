setwd("/Users/haleybranch/Desktop/LICOR 6800")
# I want S11, S15, S32, S07 for years 2010 and 2016
mydata <- na.omit(data.frame(read.csv("mimulusjuly2018.csv")))
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

## Model for fixed and random effects 
# fixed effect = year*treatment
# random effect = block 
mod1 = lm(A ~ Treatment+Plant.ID*Block, mydata.finalsub)
summary(mod1)
histogram(resid(mod1), breaks=20)
qqnorm(resid(mod1))
Anova(mod1, type=3)








