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

