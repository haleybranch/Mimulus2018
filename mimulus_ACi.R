mimulus_aci <- data.frame(read.csv("mimulus_aci.csv"))
head(mimulus_aci, 5)
mimulus_aci<- mimulus_aci[1:2880,]
library(plantecophys)
library(dplyr)
library(tidyr)
library(tidyverse) # for data manipulation
library(lme4) # for mixed models
library(lmtest) # for likelihood ratio tests
library(visreg) # for visualizing effects
library(lmerTest)
library(MuMIn)

# combining columns into single string
aci.data <- unite_(mimulus_aci, "Curve", c("Site", "Plant", "Year", "Treatment"))

# new dataframe with Curve group, Ci, A, Tleaf, PARi
aci.sub <- aci.data %>%
  select(Curve, Ci, A, Tleaf, PARi)

# change name of A to Photo
aci.sub$Photo <- aci.sub$A
aci.sub$A = NULL
aci.fit <- fitacis(aci.sub, "Curve")

# remove S32_335_2016_D from dataset -- weird reading where Ci is -2040 and A was 451
# skipped S36 340, 341, 342 because the data was incorrect 
aci.fit[["S32_335_2016_W"]] <- NULL
aci.fit[["S32_335_2016_D"]] <- NULL
aci.fit[["S36_340_2011_W"]] <- NULL
aci.fit[["S36_340_2011_D"]] <- NULL
aci.fit[["S36_341_2011_W"]] <- NULL
aci.fit[["S36_341_2011_D"]] <- NULL
aci.fit[["S36_342_2011_W"]] <- NULL
aci.fit[["S36_342_2011_D"]] <- NULL

# also remove these from unique Curves list for proper indexing below
aci.sub <- aci.sub %>% 
  filter(Curve != "S32_335_2016_W") %>% 
  filter(Curve!="S32_335_2016_D") %>% 
  filter(Curve!="S36_340_2011_W") %>% 
  filter(Curve!="S36_340_2011_D") %>% 
  filter(Curve!="S36_341_2011_W") %>% 
  filter(Curve!="S36_341_2011_D") %>% 
  filter(Curve!="S36_342_2011_W") %>% 
  filter(Curve!="S36_342_2011_D")

# plotting
plot(aci.fit[["S11_119_2010_D"]], what = c("data"))
plot(aci.fit, "oneplot")


#To look at specific plants

aci.fit[["S11_119_2010_D"]]
plot(aci.fit[["S11_119_2010_D"]])
aci.fit[["S11_119_2010_W"]]
plot(aci.fit[["S11_119_2010_W"]])

plot(aci.fit[["S18_245_2010_D"]])


# getting the data from each individual aci curve 
x <- aci.fit[["S18_278_2016_D"]]
x2 <- x$pars[,-2]
x2
x$Ci_transition

x$df
plot(x)

# All ACi outputs into a table
#Creating a loop for Vcmax, Jmax, Rd, Ci.transition dataframe from ACI Curves



ID_list <- matrix(unique((aci.sub$Curve))) 
datalist = list()

for (i in 1:length(ID_list)) {
  dat <- data.frame(ID_list[[i]], aci.fit[[i]]$pars[1,-2],aci.fit[[i]]$pars[2,-2],aci.fit[[i]]$pars[3,-2])#add here for column
  dat$Ci.transition <- aci.fit[[i]]$Ci_transition
  dat$i <- i
  datalist[[i]] <- dat # add it to your list
}

aci_output = do.call(rbind, datalist)
names(aci_output)= c("ID","Vcmax","Jmax","Rd", "Ci.transition", "i")#add names here before "i" 
aci_output$i <- NULL

# replicate ID column
aci_output$ID2 <- aci_output$ID
## new columns with Site, Year, Treatment
aci_output <- separate(aci_output, ID2, into = c("Site", "ID3", "Year", "Treatment"), sep="_")

aci_output

#Add in climate and weather covariates
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
aci_output$Year <- as.numeric(aci_output$Year)
aci_output$Site <- as.factor(aci_output$Site)
mydata <- left_join(aci_output, wna_all, by=c("Site", "Year"))

# Scale variables before running models
mydata <- mydata %>% 
  mutate(#Year.scaled = scale(Year), #treat year as category instead?
    Latitude.scaled = as.vector(scale(Latitude)),
    MAT.clim.scaled = as.vector(scale(MAT.clim)),
    MAP.clim.scaled = as.vector(scale(MAP.clim)),
    CMD.clim.scaled = as.vector(scale(CMD.clim)),    
    MAT.weath.scaled = as.vector(scale(MAT.weath)),
    MAP.weath.scaled = as.vector(scale(MAP.weath)),
    CMD.weath.scaled = as.vector(scale(CMD.weath)),    
    MAT.anom.scaled = as.vector(scale(MAT.anom)),
    MAP.anom.scaled = as.vector(scale(MAP.anom)),
    CMD.anom.scaled = as.vector(scale(CMD.anom)))

plot(CMD.anom.scaled~CMD.clim.scaled, data=mydata)
plot(CMD.weath.scaled~CMD.clim.scaled, data=mydata)

# make new column grouping pre-drought (2010 and 2011) and post drought (2015-2016)
mydata <- mydata %>%
  mutate(PrePost = ifelse(Year == "2010", "Pre", 
                      ifelse(Year == "2011", "Pre", 
                         ifelse(Year == "2015", "Post",
                            ifelse(Year == "2016", "Post", NA)))))

# Make sure factors are set correctly
mydata$Year <- as.factor(mydata$Year)
mydata$Site <- as.factor(mydata$Site)
mydata$Block <- as.factor(mydata$Block)
mydata$Plant.ID <- as.factor(mydata$Plant.ID)
mydata$PrePost <- as.factor(mydata$PrePost)


######## are there any group differences, regardless of climate/weather?
Vcmax1.group= lmer(Vcmax ~ Treatment*Site*PrePost + (1|Year), mydata)
summary(Vcmax1.group)
anova(Vcmax1.group) #3-way significant at P<0.05

# drop 3-way
Vcmax2.group= lmer(Vcmax ~ Treatment*Site + Site*PrePost + Treatment*PrePost + (1|Year), mydata)
summary(Vcmax2.group)
anova(Vcmax2.group) 

model.sel(Vcmax1.group, Vcmax2.group)
#wow, model with 3-way interaction is highly favored over reduced model
visreg(Vcmax1.group, xvar="PrePost", by="Site", cond=list(Treatment="D"))
visreg(Vcmax1.group, xvar="PrePost", by="Site", cond=list(Treatment="W"))

# Can we sub Latitude for Site to make interpretation easier?
Vcmax1.lat= lmer(Vcmax ~ Treatment*Latitude*PrePost + (1|Year), mydata)
summary(Vcmax1.lat)
anova(Vcmax1.lat) #3-way not significant 
visreg(Vcmax1.lat, xvar="Latitude", by="PrePost", overlay=T, cond=list(Treatment="D"))
visreg(Vcmax1.lat, xvar="Latitude", by="PrePost", overlay=T, cond=list(Treatment="W"))

# drop 3-way
Vcmax2.lat= lmer(Vcmax ~ Treatment*Latitude + Latitude*PrePost + Treatment*PrePost + (1|Year), mydata)
summary(Vcmax2.lat)
anova(Vcmax2.lat) 

model.sel(Vcmax1.lat, Vcmax2.lat) #about 50:50 model weight, AIC almost tied 

# Can we sub CMD for Site to make interpretation easier?
Vcmax1.cmd= lmer(Vcmax ~ Treatment*CMD.clim.scaled*PrePost + (1|Year), mydata)
summary(Vcmax1.cmd)
anova(Vcmax1.cmd) #3-way not significant 
visreg(Vcmax1.cmd, xvar="CMD.clim.scaled", by="PrePost", overlay=T, cond=list(Treatment="D"))
visreg(Vcmax1.lat, xvar="CMD.clim.scaled", by="PrePost", overlay=T, cond=list(Treatment="W"))

# drop 3-way
Vcmax2.cmd= lmer(Vcmax ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*PrePost + Treatment*PrePost + (1|Year), mydata)
summary(Vcmax2.cmd)
anova(Vcmax2.cmd) 

model.sel(Vcmax1.cmd, Vcmax2.cmd) #strong support for 3-way interaction

##### CMD for Vcmax with anomaly
Vcmax1.cmdanom <- lmer(Vcmax ~ Treatment*CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
summary(Vcmax1.cmdanom)
anova(Vcmax1.cmdanom)

# drop 3-way
Vcmax2.cmdanom <- lmer(Vcmax ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
summary(Vcmax2.cmdanom)
lrtest(Vcmax1.cmdanom, Vcmax2.cmdanom) # this model is better than more complicated Vcmax1; simplify to this one

# Drop 2-ways singly
##Drop Trt*anomaly
Vcmax3.cmdanom <- lmer(Vcmax ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax2.cmdanom,Vcmax3.cmdanom) #this model is not better or worse than more complicated Vcmax2, no support for retaining Trt*anom
## Drop Trt*climate
Vcmax4.cmdanom <- lmer(Vcmax ~ Treatment*CMD.anom.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax2.cmdanom,Vcmax4.cmdanom) # Vcmax4 is not better or worse than Vcmax2, drop Trt*clim
Vcmax4b.cmdanom <- lmer(Vcmax ~ CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax3.cmdanom, Vcmax4b.cmdanom) # Vcmax4b is not better or worse than Vcmax3, drop Trt*clim

## Drop anom*climate
Vcmax5.cmdanom <- lmer(Vcmax ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax2.cmdanom, Vcmax5.cmdanom) #Vcmax5 is no better or worse than Vcmax2, no support for retaining clim*anom

# Go down to main effects only
Vcmax6.cmdanom <- lmer(Vcmax ~ Treatment + CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), mydata)
# Drop Treatment
Vcmax7.cmdanom <- lmer(Vcmax ~  CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax7.cmdanom, Vcmax6.cmdanom) #Treatment not significant, keep it out

Vcmax8.cmdanom <- lmer(Vcmax ~ CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax8.cmdanom, Vcmax7.cmdanom)

Vcmax9.cmdanom <- lmer(Vcmax ~ CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax9.cmdanom, Vcmax7.cmdanom)

model.sel(Vcmax1.cmdanom, Vcmax2.cmdanom, Vcmax3.cmdanom, Vcmax4.cmdanom, Vcmax5.cmdanom, Vcmax6.cmdanom, Vcmax7.cmdanom, Vcmax8.cmdanom, Vcmax9.cmdanom)
# top 3 models have more support than any others and have similar weights to one another: Vcmax1 (3-way), Vcmax5 (2-ways except anom*clim), and Vcmax2 (all 2-ways)

visreg(Vcmax1.cmdanom, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Vcmax2.cmdanom, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Vcmax5.cmdanom, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Vcmax1.cmdanom, xvar="CMD.anom.scaled", by="Treatment", overlay=T)
visreg(Vcmax2.cmdanom, xvar="CMD.anom.scaled", by="Treatment", overlay=T)
visreg(Vcmax5.cmdanom, xvar="CMD.anom.scaled", by="Treatment", overlay=T)


########## Jmax 

# among groups
Jmax1.group= lmer(Jmax ~ Treatment*Site*PrePost + (1|Year), mydata)
summary(Jmax1.group)
anova(Jmax1.group)

Jmax2.group= lmer(Jmax ~ Treatment*Site + Treatment*PrePost + Site*PrePost + (1|Year), mydata)
summary(Jmax2.group)
anova(Jmax2.group)

model.sel(Jmax1.group, Jmax2.group) #strong support for 3-way interaction
visreg(Jmax1.group, xvar="PrePost", by="Site", cond=list(Treatment="D"))
visreg(Jmax1.group, xvar="PrePost", by="Site", cond=list(Treatment="W"))

# sub Latitude for Site
Jmax1.lat= lmer(Jmax ~ Treatment*Latitude*PrePost + (1|Year), mydata)
summary(Jmax1.lat)
anova(Jmax1.lat)

Jmax2.lat= lmer(Jmax ~ Treatment*Latitude + Treatment*PrePost + Latitude*PrePost + (1|Year), mydata)
summary(Jmax2.lat)
anova(Jmax2.lat)

model.sel(Jmax1.lat, Jmax2.lat) #weak-moderate support for 3-way interaction

# sub CMD for Latitude
Jmax1.cmd= lmer(Jmax ~ Treatment*CMD.clim.scaled*PrePost + (1|Year), mydata)
summary(Jmax1.cmd)
anova(Jmax1.cmd)

Jmax2.cmd= lmer(Jmax ~ Treatment*CMD.clim.scaled + Treatment*PrePost + CMD.clim.scaled*PrePost + (1|Year), mydata)
summary(Jmax2.cmd)
anova(Jmax2.cmd)

model.sel(Jmax1.cmd, Jmax2.cmd) #strong support for 3-way interaction

# CMD for Jmax with anomaly
Jmax1.cmdanom <- lmer(Jmax ~ Treatment*CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
summary(Jmax1.cmdanom)
anova(Jmax1.cmdanom)

# drop 3-way
Jmax2.cmdanom <- lmer(Jmax ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
summary(Jmax2.cmdanom)
lrtest(Jmax1.cmdanom, Jmax2.cmdanom) # this model is better than more complicated Jmax1; simplify to this one

# Drop 2-ways singly
##Drop Trt*anomaly
Jmax3.cmdanom <- lmer(Jmax ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Jmax2.cmdanom,Jmax3.cmdanom) #this model is better than more complicated Jmax2, drop Trt*anom

## Drop Trt*climate
Jmax4.cmdanom <- lmer(Jmax ~ Treatment*CMD.anom.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Jmax2.cmdanom, Jmax4.cmdanom) # Jmax4 is  better or worse than Jmax2, drop Trt*clim

## Drop anom*climate
Jmax5.cmdanom <- lmer(Jmax ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Jmax2.cmdanom,Jmax5.cmdanom) #Jmax5 is no better or worse than Jmax2, no support for retaining clim*anom

# Go down to main effects only
Jmax6.cmdanom <- lmer(Jmax ~ Treatment + CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), mydata)
# Drop Treatment
Jmax7.cmdanom <- lmer(Jmax ~  CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Jmax7.cmdanom, Jmax6.cmdanom) #Treatment is significant, keep it

# Drop Climate
Jmax8.cmdanom <- lmer(Jmax ~ Treatment + CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Jmax8.cmdanom, Jmax6.cmdanom) #Climate is not better or worse than simpler model, so drop it

# Drop Anomaly
Jmax9.cmdanom <- lmer(Jmax ~ Treatment + (1|Year) + (1|Site), mydata)
lrtest(Jmax9.cmdanom, Jmax8.cmdanom) #Anomaly is better than simpler model, so keep it

model.sel(Jmax1.cmdanom, Jmax2.cmdanom, Jmax3.cmdanom, Jmax4.cmdanom, Jmax5.cmdanom, Jmax6.cmdanom, Jmax7.cmdanom, Jmax8.cmdanom, Jmax9.cmdanom)

# top 3 models have more support than any others and have similar weights to one another: Jmax1 (3-way), Jmax1 (2-ways except anom*clim), and Jmax1 (all 2-ways)

visreg(Jmax1.cmd, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Jmax2.cmd, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Jmax5.cmd, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Jmax1.cmd, xvar="CMD.anom.scaled", by="Treatment", overlay=T)
visreg(Jmax2.cmd, xvar="CMD.anom.scaled", by="Treatment", overlay=T)
visreg(Jmax5.cmd, xvar="CMD.anom.scaled", by="Treatment", overlay=T)


### Amy stopped here - can follow similar workflow as above for Rd and Ci transition



### Rd
Rd1.cmd= lmer(Rd ~ Treatment*CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
summary(Rd1.cmd)
anova(Rd1.cmd)

# drop 3-way
Rd2.cmd <- lmer(Rd ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
summary(Rd2.cmd)
lrtest(Rd1.cmd, Rd2.cmd) # this model is better than more complicated Rd1; simplify to this one

# Drop 2-ways singly
##Drop Trt*anomaly
Rd3.cmd <- lmer(Rd ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Rd2.cmd,Rd3.cmd) #this model is better than more complicated Rd2, drop Trt*anom
## Drop Trt*climate
Rd4.cmd <- lmer(Rd ~ Treatment*CMD.anom.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Rd2.cmd,Rd4.cmd) # Rd4 is no better or worse than Rd2, drop Trt*clim

## Drop anom*climate
Rd5.cmd <- lmer(Rd ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Rd2.cmd,Rd5.cmd) #Rd5 is better than Rd2, drop clim*anom

# Go down to main effects only
Rd6.cmd <- lmer(Rd ~ Treatment + CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), mydata)
# Drop Treatment
Rd7.cmd <- lmer(Rd ~  CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Rd7.cmd, Rd6.cmd) #Treatment is not significant, remove

# Drop Climate
Rd8.cmd <- lmer(Rd ~ CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Rd8.cmd, Rd7.cmd) #Climate is better than simpler model
# Drop Anomaly
Rd9.cmd <- lmer(Rd ~ CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Rd9.cmd, Rd7.cmd) #Anomaly is better than simpler model, so keep it

summary(Rd7.cmd) # not significant

### Ci transition 
Ci1.cmd= lmer(Ci.transition ~ Treatment*CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
summary(Ci1.cmd)
anova(Ci1.cmd)

# drop 3-way
Ci2.cmd <- lmer(Ci.transition ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
summary(Ci2.cmd)
lrtest(Ci1.cmd, Ci2.cmd) # this model is better than more complicated Ci1; simplify to this one

# Drop 2-ways singly
##Drop Trt*anomaly
Ci3.cmd <- lmer(Ci.transition ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Ci2.cmd,Ci3.cmd) #this model is better than more complicated Ci2, drop Trt*anom
## Drop Trt*climate
Ci4.cmd <- lmer(Ci.transition ~ Treatment*CMD.anom.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Ci2.cmd,Ci4.cmd) # Ci4 is better than Ci2, drop Trt*clim

## Drop anom*climate
Ci5.cmd <- lmer(Ci.transition ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Ci2.cmd,Ci5.cmd) #Ci5 is better than Ci2, drop clim*anom

# Go down to main effects only
Ci6.cmd <- lmer(Ci.transition ~ Treatment + CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), mydata)
# Drop Treatment
Ci7.cmd <- lmer(Ci.transition ~  CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Ci7.cmd, Ci6.cmd) #Treatment is significant, keep it in

# Drop Climate
Ci8.cmd <- lmer(Ci.transition ~ Treatment + CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Ci8.cmd, Ci6.cmd) #Climate is better than simpler model, keep it
# Drop Anomaly
Ci9.cmd <- lmer(Ci.transition ~ Treatment + (1|Year) + (1|Site), mydata)
lrtest(Ci9.cmd, Ci8.cmd) #Anomaly is better than simpler model, so keep it

library(lmerTest)
summary(Ci8.cmd)
visreg(Ci8.cmd, xvar="Treatment", ylab = "Ci.transition") #No significant difference



### trying to get mesophyll conductance
Photosyn(VPD = 1.5, Ca = 400, PPFD = 1500, Tleaf = 25,
         Patm = 100, RH = NULL, gsmodel = c("BBOpti", "BBLeuning",
                                            "BallBerry", "BBdefine"), g1 = 4, g0 = 0, gk = 0.5, vpdmin = 0.5,
         D0 = 5, GS = NULL, BBmult = NULL, alpha = 0.24, theta = 0.85,
         Jmax = 100, Vcmax = 50, gmeso = TRUE, TPU = 1000, alphag = 0,
         Rd0 = 0.92, Q10 = 1.92, Rd = NULL, TrefR = 25, Rdayfrac = 1,
         EaV = 58550, EdVC = 2e+05, delsC = 629.26, EaJ = 29680,
         EdVJ = 2e+05, delsJ = 631.88, GammaStar = NULL, Km = NULL,
         Ci = NULL, Tcorrect = TRUE, returnParsOnly = FALSE,
         whichA = c("Ah", "Amin", "Ac", "Aj"))

gmeso[["S11_119_2010_D"]]
