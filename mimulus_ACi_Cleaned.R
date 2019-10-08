#mimulus_ACi cleaned up version 
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

# combining columns into single string and replaces
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
acidata <- left_join(aci_output, wna_all, by=c("Site", "Year"))

# Scale variables before running models
acidata <- acidata %>% 
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

plot(CMD.anom.scaled~CMD.clim.scaled, data=acidata)
plot(CMD.weath.scaled~CMD.clim.scaled, data=acidata)

# make new column grouping pre-drought (2010 and 2011) and post drought (2015-2016)
acidata <- acidata %>%
  mutate(PrePost = ifelse(Year == "2010", "Pre", 
                          ifelse(Year == "2011", "Pre", 
                                 ifelse(Year == "2015", "Post",
                                        ifelse(Year == "2016", "Post", NA)))))

# Make sure factors are set correctly
acidata$Year <- as.factor(acidata$Year)
acidata$Site <- as.factor(acidata$Site)
acidata$Block <- as.factor(acidata$Block)
acidata$Plant.ID <- as.factor(acidata$Plant.ID)
acidata$PrePost <- as.factor(acidata$PrePost)

######## are there any group differences, regardless of climate/weather?
Vcmax1.group= lmer(Vcmax ~ Treatment*Site*PrePost + (1|Year), acidata)
summary(Vcmax1.group)
anova(Vcmax1.group) #3-way significant at P<0.05

# drop 3-way
Vcmax2.group= lmer(Vcmax ~ Treatment*Site + Site*PrePost + Treatment*PrePost + (1|Year), acidata)
summary(Vcmax2.group)
anova(Vcmax2.group) 

model.sel(Vcmax1.group, Vcmax2.group) #wow, model with 3-way interaction is highly favored over reduced model
visreg(Vcmax1.group, xvar="PrePost", by="Site", cond=list(Treatment="D"))
visreg(Vcmax1.group, xvar="PrePost", by="Site", cond=list(Treatment="W"))
visreg(Vcmax1.group, xvar="PrePost", by="Site", cond=list(Treatment=c("W", "D")))

###CMD for Vcmax with anomaly
# substituted site with CMD.clim.scaled and PrePost with CMD.anom.scaled
Vcmax1.cmdanom <- lmer(Vcmax ~ Treatment*CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
summary(Vcmax1.cmdanom)
anova(Vcmax1.cmdanom) #no significance 

# drop 3-way
Vcmax2.cmdanom <- lmer(Vcmax ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
summary(Vcmax2.cmdanom)
lrtest(Vcmax1.cmdanom, Vcmax2.cmdanom) # this model is not better or worse than more complicated Vcmax1; simplify to this one

# Drop 2-ways singly
##Drop Trt*anomaly
Vcmax3.cmdanom <- lmer(Vcmax ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
lrtest(Vcmax2.cmdanom,Vcmax3.cmdanom) #this model is better than more complicated Vcmax2, drop Trt*anom
## Drop Trt*climate
Vcmax4.cmdanom <- lmer(Vcmax ~ Treatment*CMD.anom.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
lrtest(Vcmax2.cmdanom,Vcmax4.cmdanom) # Vcmax4 is not better or worse than Vcmax2, drop Trt*clim
Vcmax4b.cmdanom <- lmer(Vcmax ~ CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
lrtest(Vcmax3.cmdanom, Vcmax4b.cmdanom) # Vcmax4b is not better or worse than Vcmax3, drop Trt*clim

## Drop anom*climate
Vcmax5.cmdanom <- lmer(Vcmax ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + (1|Year) + (1|Site), acidata)
lrtest(Vcmax2.cmdanom, Vcmax5.cmdanom) #Vcmax5 is no better or worse than Vcmax2, no support for retaining clim*anom

# Go down to main effects only
Vcmax6.cmdanom <- lmer(Vcmax ~ Treatment + CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), acidata)
# Drop Treatment
Vcmax7.cmdanom <- lmer(Vcmax ~  CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), acidata)
lrtest(Vcmax7.cmdanom, Vcmax6.cmdanom) #Treatment not significant, keep it out

Vcmax8.cmdanom <- lmer(Vcmax ~ CMD.anom.scaled + (1|Year) + (1|Site), acidata)
lrtest(Vcmax8.cmdanom, Vcmax7.cmdanom)

Vcmax9.cmdanom <- lmer(Vcmax ~ CMD.clim.scaled + (1|Year) + (1|Site), acidata)
lrtest(Vcmax9.cmdanom, Vcmax7.cmdanom)

model.sel(Vcmax1.cmdanom, Vcmax2.cmdanom, Vcmax3.cmdanom, Vcmax4.cmdanom, Vcmax5.cmdanom, Vcmax6.cmdanom, Vcmax7.cmdanom, Vcmax8.cmdanom, Vcmax9.cmdanom)
# top 3 models have more support than any others and have similar weights to one another: Vcmax1 (3-way), Vcmax5 (2-ways except anom*clim), and Vcmax2 (all 2-ways)

visreg(Vcmax1.cmdanom, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Vcmax2.cmdanom, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Vcmax5.cmdanom, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Vcmax1.cmdanom, xvar="CMD.anom.scaled", by="Treatment", overlay=T)
visreg(Vcmax2.cmdanom, xvar="CMD.anom.scaled", by="Treatment", overlay=T)
visreg(Vcmax5.cmdanom, xvar="CMD.anom.scaled", by="Treatment", overlay=T)

# Overall: There is a significant interaction between treatment:site:prepost, but no main effects 

# Jmax

# among groups
Jmax1.group= lmer(Jmax ~ Treatment*Site*PrePost + (1|Year), acidata)
summary(Jmax1.group)
anova(Jmax1.group) #no significance 

Jmax2.group= lmer(Jmax ~ Treatment*Site + Treatment*PrePost + Site*PrePost + (1|Year), acidata)
summary(Jmax2.group)
anova(Jmax2.group) #no significance 

model.sel(Jmax1.group, Jmax2.group) #strong support for 3-way interaction
visreg(Jmax1.group, xvar="PrePost", by="Site", cond=list(Treatment="D"))
visreg(Jmax1.group, xvar="PrePost", by="Site", cond=list(Treatment="W"))

# CMD for Jmax with anomaly
Jmax1.cmdanom <- lmer(Jmax ~ Treatment*CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
summary(Jmax1.cmdanom)
anova(Jmax1.cmdanom)

# drop 3-way
Jmax2.cmdanom <- lmer(Jmax ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
summary(Jmax2.cmdanom)
lrtest(Jmax1.cmdanom, Jmax2.cmdanom) # this model is better than more complicated Jmax1; simplify to this one

# Drop 2-ways singly
##Drop Trt*anomaly
Jmax3.cmdanom <- lmer(Jmax ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
lrtest(Jmax2.cmdanom,Jmax3.cmdanom) #this model is better than more complicated Jmax2, drop Trt*anom

## Drop Trt*climate
Jmax4.cmdanom <- lmer(Jmax ~ Treatment*CMD.anom.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
lrtest(Jmax2.cmdanom, Jmax4.cmdanom) # Jmax4 is  better than Jmax2, drop Trt*clim

## Drop anom*climate
Jmax5.cmdanom <- lmer(Jmax ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + (1|Year) + (1|Site), acidata)
lrtest(Jmax2.cmdanom,Jmax5.cmdanom) #Jmax5 is no better or worse than Jmax2, no support for retaining clim*anom

# Go down to main effects only
Jmax6.cmdanom <- lmer(Jmax ~ Treatment + CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), acidata)
# Drop Treatment
Jmax7.cmdanom <- lmer(Jmax ~  CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), acidata)
lrtest(Jmax7.cmdanom, Jmax6.cmdanom) #Treatment is significant, keep it

# Drop Climate
Jmax8.cmdanom <- lmer(Jmax ~ Treatment + CMD.anom.scaled + (1|Year) + (1|Site), acidata)
lrtest(Jmax8.cmdanom, Jmax6.cmdanom) #Climate is not better or worse than simpler model, so drop it

# Drop Anomaly
Jmax9.cmdanom <- lmer(Jmax ~ Treatment + (1|Year) + (1|Site), acidata)
lrtest(Jmax9.cmdanom, Jmax8.cmdanom) #Anomaly is better than simpler model, so keep it

model.sel(Jmax1.cmdanom, Jmax2.cmdanom, Jmax3.cmdanom, Jmax4.cmdanom, Jmax5.cmdanom, Jmax6.cmdanom, Jmax7.cmdanom, Jmax8.cmdanom, Jmax9.cmdanom)

# top 3 models have more support than any others and have similar weights to one another: Jmax1 (3-way), Jmax1 (2-ways except anom*clim), and Jmax1 (all 2-ways)
# Help: confused by this part. Jmax1 (3-way) is most supported with AIC, but through dropping we found Jmax8 had the most support 

visreg(Jmax1.cmd, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Jmax2.cmd, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Jmax5.cmd, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Jmax1.cmd, xvar="CMD.anom.scaled", by="Treatment", overlay=T)
visreg(Jmax2.cmd, xvar="CMD.anom.scaled", by="Treatment", overlay=T)
visreg(Jmax5.cmd, xvar="CMD.anom.scaled", by="Treatment", overlay=T)

#Overall: no significant affect 

### Rd
Rd1.group = lmer(Rd ~ Treatment*Site*PrePost + (1|Year), acidata)
summary(Rd1.group)
anova(Rd1.group)

Rd2.group= lmer(Rd ~ Treatment*Site + Treatment*PrePost + Site*PrePost + (1|Year), acidata)
summary(Rd2.group)
anova(Rd2.group)

model.sel(Rd1.group, Rd2.group) #strong support for 2-way interaction
visreg(Rd1.group, xvar="PrePost", by="Site", cond=list(Treatment="D"))
visreg(Rd1.group, xvar="PrePost", by="Site", cond=list(Treatment="W"))

# CMD for Jmax with anomaly
Rd1.cmdanom <- lmer(Rd ~ Treatment*CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
summary(Rd1.cmdanom)
anova(Rd1.cmdanom)

# drop 3-way
Rd2.cmdanom <- lmer(Rd ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
summary(Rd2.cmdanom)
lrtest(Rd1.cmdanom, Rd2.cmdanom) # this model is better than more complicated Rd1; simplify to this one

# Drop 2-ways singly
##Drop Trt*anomaly
Rd3.cmdanom <- lmer(Rd ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
lrtest(Rd2.cmdanom,Rd3.cmdanom) #this model is not better or worse than more complicated Rd2, drop Trt*anom
## Drop Trt*climate
Rd4.cmdanom <- lmer(Rd ~ Treatment*CMD.anom.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
lrtest(Rd3.cmdanom,Rd4.cmdanom) # Rd4 is better than Rd3, drop Trt*clim

## Drop anom*climate
Rd5.cmdanom <- lmer(Rd ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + (1|Year) + (1|Site), acidata)
lrtest(Rd2.cmdanom,Rd5.cmdanom) #Rd5 is better than Rd2, drop clim*anom

# Go down to main effects only
Rd6.cmdanom <- lmer(Rd ~ Treatment + CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), acidata)
# Drop Treatment
Rd7.cmdanom <- lmer(Rd ~  CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), acidata)
lrtest(Rd7.cmdanom, Rd6.cmdanom) #Treatment is not significant, remove

# Drop Climate
Rd8.cmdanom <- lmer(Rd ~ CMD.anom.scaled + (1|Year) + (1|Site), acidata)
lrtest(Rd8.cmdanom, Rd7.cmdanom) #Climate is better than simpler model
# Drop Anomaly
Rd9.cmdanom <- lmer(Rd ~ CMD.clim.scaled + (1|Year) + (1|Site), acidata)
lrtest(Rd9.cmdanom, Rd7.cmdanom) #Anomaly is better than simpler model, so keep it

model.sel(Rd1.cmdanom, Rd2.cmdanom, Rd3.cmdanom, Rd4.cmdanom, Rd5.cmdanom, Rd6.cmdanom, Rd7.cmdanom, Rd8.cmdanom, Rd9.cmdanom)
# best models: Rd9, Rd8, Rd7
# don't understand when we narrowed it down to Rd7 
anova(Rd7.cmdanom) #not significant 

visreg(Rd1.cmd, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Rd2.cmd, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Rd5.cmd, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Rd1.cmd, xvar="CMD.anom.scaled", by="Treatment", overlay=T)
visreg(Rd2.cmd, xvar="CMD.anom.scaled", by="Treatment", overlay=T)
visreg(Rd5.cmd, xvar="CMD.anom.scaled", by="Treatment", overlay=T)


### Ci
Ci1.group = lmer(Ci.transition ~ Treatment*Site*PrePost + (1|Year), acidata)
summary(Ci1.group)
anova(Ci1.group)

Ci2.group= lmer(Ci.transition ~ Treatment*Site + Treatment*PrePost + Site*PrePost + (1|Year), acidata)
summary(Ci2.group)
anova(Ci2.group)

model.sel(Ci1.group, Ci2.group) #strong support for 3-way interaction
visreg(Ci1.group, xvar="PrePost", by="Site", cond=list(Treatment="D"))
visreg(Ci1.group, xvar="PrePost", by="Site", cond=list(Treatment="W"))
visreg(Ci1.group, xvar="PrePost", by="Treatment", cond=list(Site=))

# CMD for Jmax with anomaly
Ci1.cmdanom <- lmer(Ci.transition ~ Treatment*CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
summary(Ci1.cmdanom)
anova(Ci1.cmdanom) #CMD.clim.scaled significant 

# drop 3-way
Ci2.cmdanom <- lmer(Ci.transition ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
summary(Ci2.cmdanom)
lrtest(Ci1.cmdanom, Ci2.cmdanom) # this model is better than more complicated Ci1; simplify to this one

# Drop 2-ways singly
##Drop Trt*anomaly
Ci3.cmdanom <- lmer(Ci.transition ~ Treatment*CMD.clim.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
lrtest(Ci2.cmdanom,Ci3.cmdanom) #this model is better than more complicated Rd2, drop Trt*anom
## Drop Trt*climate
Ci4.cmdanom <- lmer(Ci.transition ~ Treatment*CMD.anom.scaled + CMD.clim.scaled*CMD.anom.scaled + (1|Year) + (1|Site), acidata)
lrtest(Ci2.cmdanom,Ci4.cmdanom) # Ci4 is better, drop Trt*clim

## Drop anom*climate
Ci5.cmdanom <- lmer(Ci.transition ~ Treatment*CMD.anom.scaled + Treatment*CMD.clim.scaled + (1|Year) + (1|Site), acidata)
lrtest(Ci2.cmdanom,Ci5.cmdanom) #Ci5 is better, drop clim*anom

# Go down to main effects only
Ci6.cmdanom <- lmer(Ci.transition ~ Treatment + CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), acidata)
# Drop Treatment
Ci7.cmdanom <- lmer(Ci.transition ~  CMD.anom.scaled + CMD.clim.scaled + (1|Year) + (1|Site), acidata)
lrtest(Ci7.cmdanom, Ci6.cmdanom) #Treatment is significant, keep it

# Drop Climate
Ci8.cmdanom <- lmer(Ci.transition ~ Treatment + CMD.anom.scaled + (1|Year) + (1|Site), acidata)
lrtest(Ci8.cmdanom, Ci7.cmdanom) #Climate is better than simpler model
# Drop Anomaly
Ci9.cmdanom <- lmer(Ci.transition ~ CMD.clim.scaled + (1|Year) + (1|Site), acidata)
lrtest(Ci9.cmdanom, Ci7.cmdanom) #Anomaly is better than simpler model, so keep it

#model Ci6 is best 
anova(Ci6.cmdanom) # significant interaction between Ci and clim scaled (corresponds to site)
visreg(Ci6.cmdanom, xvar="CMD.clim.scaled", by="Treatment", overlay=T) 
# Drier sites have higher transition point than wetter sites 
visreg(Ci6.cmdanom, xvar="CMD.anom.scaled", by="Treatment", overlay=T)


model.sel(Ci1.cmdanom, Ci2.cmdanom, Ci3.cmdanom, Ci4.cmdanom, Ci5.cmdanom, Ci6.cmdanom, Ci7.cmdanom, Ci8.cmdanom, Ci9.cmdanom)
#that's not what we got here though for model select 

visreg(Ci1.cmdanom, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Ci2.cmdanom, xvar="CMD.clim.scaled", by="Treatment", overlay=T)
visreg(Ci5.cmdanom, xvar="CMD.clim.scaled", by="Treatment", overlay=T)


####### Adding a graph showing the ACi curves 
ggplot(data = before, aes(x = Ci, y = Photo)) +
  geom_point(aes(colour = factor(Site))) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(title = NULL))+
  ylab("Assimilation")+
  xlab("Ci")

# new column for aci.sub 
before <- aci.sub

before <- separate(aci.sub, Curve, into = c("Site", "ID3", "Year", "Treatment"), sep="_")

before

ggplot(subset(before,Site %in% c("S07")), aes(x = Ci, y = Photo)) + 
  geom_point() +
  geom_smooth() +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(title = NULL))+
  ylab("Assimilation")+
  xlab("Ci")

# Plotting each aci curve
plot((aci.fit[["S18_278_2016_D"]]), main="S18_278_2016_D")


