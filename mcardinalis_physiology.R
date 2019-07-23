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
########### LRTEST = higher number is the better model 
########### AIC = lower number is the better model 
mod1.cmd= lmer(A ~ Treatment*CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
summary(mod1.cmd)
anova(mod1.cmd) # 2 way Treatment*CMD.clim.scaled significant 

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
anova(gsw1.cmd) #3-way 0.05

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
lrtest(gsw2.cmd,gsw4.cmd) # gsw4 is not better or worse than gsw2, drop Trt*clim
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
lrtest(gsw8.cmd, gsw6.cmd) #Climate is  better than simpler model, keep it
# Drop Anomaly
gsw9.cmd <- lmer(gsw ~ Treatment + CMD.clim.scaled + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(mod9.cmd, mod6.cmd) #Anomaly is not better or worse than simpler model, so drop it
gsw9b.cmd <- lmer(gsw ~ Treatment + (1|Year) + (1|Site/Plant.ID) + (1|Block), mydata)
lrtest(gsw9b.cmd, gsw9.cmd) #keep climate

visreg(gsw9.cmd, xvar="Treatment") #Dry treatment significantly higher stomatal conductance rate than wet across plants from years

# p value - Scatterthwaite's approx p value and CI
library(lmerTest)
# CIs <- confint(model, method = "profile") # gives confidence intervals 
CIs <- confint(gsw9.cmd, method = "profile")
summary(gsw9.cmd) ## pvalue is <2e-16 for gsw9.cmd --> treatment significant

# Adding this section to see what happens when we group years before,during,and post drought
## make new column grouping pre-drought (2010 and 2011) and post drought (2015-2016)
mydata <- mydata %>%
  mutate(PrePost = ifelse(Year == "2010", "1", 
                          ifelse(Year == "2011", "1", 
                                 ifelse(Year == "2012", "2",
                                 ifelse(Year == "2013", "2",
                                 ifelse(Year == "2014", "2",
                                 ifelse(Year == "2015", "3",
                                        ifelse(Year == "2016", "3", NA))))))))
mydata$PrePost <- as.factor(mydata$PrePost)


## are there any group differences, regardless of climate/weather?
gsw1.group= lmer(gsw ~ Treatment*Site*PrePost + (1|Year), mydata)
summary(gsw1.group)
anova(gsw1.group) #treatment, site:prepost significant, trt:site

# drop 3-way
gsw2.group= lmer(gsw ~ Treatment*Site + Site*PrePost + Treatment*PrePost + (1|Year), mydata)
summary(gsw2.group)
anova(gsw2.group) 

model.sel(gsw1.group, gsw2.group)

#reduced model favored 
visreg(gsw2.group, xvar="PrePost", by="Site", cond=list(Treatment="D"))
visreg(gsw2.group, xvar="PrePost", by="Site", cond=list(Treatment="W"))
visreg(gsw2.group, xvar="PrePost", by="Site", cond=list(Treatment=c("W", "D")))

# Can we sub Latitude for Site to make interpretation easier? 
# Don't use this one
gsw1.lat= lmer(gsw ~ Treatment*Latitude*PrePost + (1|Year), mydata)
summary(gsw1.lat)
anova(gsw1.lat) #3-way not significant 
visreg(gsw1.lat, xvar="Latitude", by="PrePost", overlay=T, cond=list(Treatment="D"))
visreg(gsw1.lat, xvar="Latitude", by="PrePost", overlay=T, cond=list(Treatment="W"))

## drop 3-way
gsw2.lat= lmer(gsw ~ Treatment*Latitude + Latitude*PrePost + Treatment*PrePost + (1|Year), mydata)
summary(gsw2.lat)
anova(gsw2.lat) 

model.sel(gsw1.lat, gsw2.lat) #about 50:50 model weight, AIC almost tied 

# Can we sub CMD for Site to make interpretation easier?
gsw1.cmd= lmer(gsw~ Treatment*CMD.clim.scaled*PrePost + (1|Year), mydata)
summary(gsw1.cmd)
anova(gsw1.cmd) #3-way not significant 
visreg(gsw1.cmd, xvar="CMD.clim.scaled", by="PrePost", overlay=T, cond=list(Treatment="D"))
visreg(gsw1.lat, xvar="CMD.clim.scaled", by="PrePost", overlay=T, cond=list(Treatment="W"))

# drop 3-way
gsw2.cmd= lmer(gsw ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*PrePost + Treatment*PrePost + (1|Year), mydata)
summary(gsw2.cmd)
anova(gsw2.cmd) 

model.sel(gsw1.cmd, gsw2.cmd) #strong support for 3-way interaction


