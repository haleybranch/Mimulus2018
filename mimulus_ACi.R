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

for (i in 1:120) {
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

# make new column grouping pre-drought (2010 and 2011) and post drought (2015-2016)
mydata <- mydata %>%
  mutate(prepost = ifelse(Year = "2010", "pre", 
                      ifelse(Year = "2011", "pre", 
                         ifelse(Year = "2015", "post",
                            ifelse(Year = "2016", "post", NA)))))

# Make sure factors are set correctly
mydata$Year <- as.factor(mydata$Year)
mydata$Site <- as.factor(mydata$Site)
mydata$Block <- as.factor(mydata$Block)
mydata$Plant.ID <- as.factor(mydata$Plant.ID)




##### work in progress
aci_est$Treatment <- as.character(aci_est$Treatment)
aci_est$Site <- as.character(aci_est$Site)

ggplot(data=aci_est, aes(x=Year, y=Vcmax)) +
  geom_point(aes(color="Treatment")) +
  facet_wrap(aes(group="Site"))

str(aci_est)

aci_est$Treatment <- as.factor(aci_est$Treatment)
aci_est$Site <- as.factor(aci_est$Site)

###################

# CMD for Vcmax
Vcmax1.cmd= lmer(Vcmax ~ Treatment*CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
summary(Vcmax1.cmd)
anova(Vcmax1.cmd)

# drop 3-way
Vcmax2.cmd <- lmer(Vcmax ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
summary(Vcmax2.cmd)
lrtest(Vcmax1.cmd, Vcmax2.cmd) # this model is better than more complicated Vcmax1; simplify to this one

# Drop 2-ways singly
##Drop Trt*anomaly
Vcmax3.cmd <- lmer(Vcmax ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax2.cmd,Vcmax3.cmd) #this model is not better or worse than more complicated Vcmax2, no support for retaining Trt*anom
## Drop Trt*climate
Vcmax4.cmd <- lmer(Vcmax ~ Treatment*CMD.anom.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax2.cmd,Vcmax4.cmd) # Vcmax4 is not better or worse than Vcmax2, drop Trt*clim
Vcmax4b.cmd <- lmer(Vcmax ~ CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax3.cmd,Vcmax4b.cmd) # Vcmax4b is not better or worse than Vcmax3, drop Trt*clim

## Drop anom*climate
Vcmax5.cmd <- lmer(Vcmax ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax2.cmd,Vcmax5.cmd) #Vcmax5 is no better or worse than Vcmax2, no support for retaining clim*anom

# Go down to main effects only
Vcmax6.cmd <- lmer(Vcmax ~ Treatment + CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), mydata)
# Drop Treatment
Vcmax7.cmd <- lmer(Vcmax ~  CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax7.cmd, Vcmax6.cmd) #Treatment not significant, keep it out

Vcmax8.cmd <- lmer(Vcmax ~ CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax8.cmd, Vcmax7.cmd)

Vcmax9.cmd <- lmer(Vcmax ~ CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Vcmax9.cmd, Vcmax7.cmd)

library(lmerTest)
summary(Vcmax7.cmd) # no significant difference
visreg(Vcmax6.cmd, xvar="Treatment", ylab = "Vcmax")

ggplot(data = mydata, aes(x=Latitude, y=Vcmax, color = Year)) + geom_point() 

### Jmax 
Jmax1.cmd= lmer(Jmax ~ Treatment*CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
summary(Jmax1.cmd)
anova(Jmax1.cmd)

# drop 3-way
Jmax2.cmd <- lmer(Jmax ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
summary(Jmax2.cmd)
lrtest(Jmax1.cmd, Jmax2.cmd) # this model is better than more complicated Jmax1; simplify to this one

# Drop 2-ways singly
##Drop Trt*anomaly
Jmax3.cmd <- lmer(Jmax ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Jmax2.cmd,Jmax3.cmd) #this model is better than more complicated Jmax2, drop Trt*anom
## Drop Trt*climate
Jmax4.cmd <- lmer(Jmax ~ Treatment*CMD.anom.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Jmax2.cmd,Jmax4.cmd) # Jmax4 is  better or worse than Jmax2, drop Trt*clim

## Drop anom*climate
Jmax5.cmd <- lmer(Jmax ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Jmax2.cmd,Jmax5.cmd) #Jmax5 is no better or worse than Jmax2, no support for retaining clim*anom

# Go down to main effects only
Jmax6.cmd <- lmer(Jmax ~ Treatment + CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), mydata)
# Drop Treatment
Jmax7.cmd <- lmer(Jmax ~  CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), mydata)
lrtest(Jmax7.cmd, Jmax6.cmd) #Treatment is significant, keep it

# Drop Climate
Jmax8.cmd <- lmer(Jmax ~ Treatment + CMD.anom.scaled + (1|Year) + (1|Site), mydata)
lrtest(Jmax8.cmd, Jmax6.cmd) #Climate is not better or worse than simpler model, so drop it
# Drop Anomaly
Jmax9.cmd <- lmer(Jmax ~ Treatment + (1|Year) + (1|Site), mydata)
lrtest(Jmax9.cmd, Jmax8.cmd) #Anomaly is better than simpler model, so keep it

library(lmerTest)
summary(Jmax8.cmd)
visreg(Jmax8.cmd, xvar="Treatment", ylab = "Jmax") #No significant difference


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
