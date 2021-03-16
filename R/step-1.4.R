# Step 1.4 of assignment

load("trial_dataset.RData");

library(fitdistrplus)

male_ages <- vector()
female_ages <- vector()

#iterate over ages and conditionally put the numbers in one of the two vectors
for (i in 1:length(data$Age)) {
  #if male
  if (data$Male[i] == 1) {
    male_ages <- c(male_ages, data$Age[i])
  }
  else {
    female_ages <- c(female_ages, data$Age[i])
  }
}

#Start with log normal as that seems like the best candidate. Age cannot be negative.

hist(male_ages)
dist_lognorm_male <- fitdist(male_ages, distr = "lnorm")
hist(female_ages)
dist_lognorm_female <- fitdist(female_ages, distr= "lnorm")

#Check
gofstat(dist_lognorm_male, fitnames=c("Male lognormal"))
plot(dist_lognorm_male)

gofstat(dist_lognorm_female, fitnames = "Female lognormal")
plot(dist_lognorm_female)

#Conclusion: lognormal is a good fit with the male age. 
#It matches the shape on the histogram and has a low Akaike’s score, indicating good fit.
#The P-P plot and CDF graph looks good as well.
# The female ages are skewed and look like they would fit better with a gamma distribution.

dist_gamma_female <- fitdist(female_ages, distr = "gamma")
plot(dist_gamma_female)
gofstat(list(dist_lognorm_female, dist_gamma_female), fitnames=c("Lognormal female", "Gamma female"))

#Gamma distribution still not a good fit

dist_weibull_female <- fitdist(female_ages, distr="weibull")
plot(dist_weibull_female)
gofstat(list(dist_weibull_female, dist_gamma_female, dist_lognorm_female),  
        fitnames=c("Weibull", "Gamma", "Lognormal"))

#The weibull distribution looks like the best fit graphically based on the histogram and Q-Q plot
# The statistics are pretty similar between the distributions (none of them are an amazing fit).
# However, the von Mises and Anderson Darling statistics are much better for the Weibull distribution. 
# These tests are appropriate because if the data is continuous (and age is continuous) then
# the tests are independent of the distribution being tested.
# The Anderson Darling test gives extra weight to the tails and this is where Weibull performs notably better.
# Small A-S statistic = better
# We will go with Weibull.

Tx1_average_response_rate <- sum(data$Tx1.C1.Dx.Pet==1) / length(data$Tx1.C1.Dx.Pet)

#Average Tx1 response rate is 48%, about 50/50

Tx1_poor_response <- vector()
for (i in 1:length(data$Poor)) {
  if (data$Poor[i] == 1) {
    Tx1_poor_response <- c(Tx1_poor_response, data$Tx1.C1.Dx.Pet[i])
  }
}
Tx1_poor_response_rate <- sum(Tx1_poor_response==1) / length(Tx1_poor_response)

Tx1_good_condition_response <- vector()
for (i in 1:length(data$Poor)) {
  if (data$Poor[i] == 0) {
    Tx1_good_condition_response <- c(Tx1_good_condition_response, data$Tx1.C1.Dx.Pet[i])
  }
}

Tx1_good_condition_response_rate <- sum(Tx1_good_condition_response==1) / length(Tx1_good_condition_response)

#44.8% vs 50%. Doesn’t make a huge amount of difference if the patient is in poor condition or not. 

#Let's see whether response in Tx2 is predicated on response in Tx1, or if they are independent.

#There should be NAs in Tx2.C1 only if the patient died sometime in the previous cycles.

#if you responded to treatment 1, how do you respond to treatment 2?
#For this, I compare Tx1.C1.Dx.Pet (first cycle in first treatment) with first cycle in second treatment
Tx2_responders <- vector()
Tx2_nonresponders <- vector()
for (i in 1:length(data$Tx1.C2.Dx.Pet)) {
  if (data$Tx1.C1.Dx.Pet[i] == 1 && !is.na(data$Tx2.C1.Dx.Pet[i])) {
    Tx2_responders <- c(Tx2_responders, data$Tx2.C1.Dx.Pet[i])
  }
  else if (!is.na(data$Tx2.C1.Dx.Pet[i])) { #we are ignoring NAs which signify death or major comps
    Tx2_nonresponders <- c(Tx2_nonresponders, data$Tx2.C1.Dx.Pet[i])
  }
}
Tx2_responders_rate <- sum(Tx2_responders==1) / length(Tx2_responders)
Tx2_nonresponders_rate <- sum(Tx2_nonresponders==1) / length(Tx2_nonresponders)

# The nonresponders had a response rate of 36% in treatment 2 cycle 1, and the responders had a rate of 59%
# I don't think a statistical test of significance is necessary with a difference this big and an n close to 200

prop.test(x=c(sum(Tx2_responders==1), sum(Tx2_nonresponders==1)), 
          n=c(length(Tx2_responders), length(Tx2_nonresponders)), p = NULL, alternative = "greater",
          correct = FALSE)
#The Z test of proportion shows us that the p value is tiny

# Logistic regression to predict death (0 or 1) based on the independent variables sex and age.
# Using data from C1 and C2 of Tx1 only.

# We will use the glm package, but first we have to define the variables

male <- data$Male #male is already 0 1 coded
age <- data$Age #age is numerical and continuous

# Let's look at the people who are still alive at C2 to figure out the total death rate during those 2 cycles.
C2_events <- data$Tx1.C2.Event
#We need to transform this data into alive or dead coded 0 or 1.
for (i in 1:length(C2_events)) {
  #0 = no comps, 2 = major, 3 = minor, but they are alive
  if (!is.na(C2_events[i]) && (C2_events[i] == 0 || C2_events[i] == 2 || C2_events[i] == 3)) {
    C2_events[i] <- 1
  }
  else {
    C2_events[i] <- 0
  }
}

logreg_death <- glm(C2_events ~ male + age, family = binomial)
summary(logreg_death)

#Age is not significant but gender is. Men are more likely to die from colon cancer.

#We want to find time-to-death and time-to-major-complication for Tx1 and Tx2 based on C1 + C2

# Time to death
C1_events <- data$Tx1.C1.Event
C1_time_to_death <- vector()
for (i in 1:length(C1_events)) {
  if (C1_events[i] == 1) { #if death
    C1_time_to_death <- c(C1_time_to_death, data$Tx1.C1.Time[i])
  }
}
C2_time_to_death <- vector()
for (i in 1:length(C2_events)) {
  if (C2_events[i] == 1) {
    C2_time_to_death <- c(C2_time_to_death, data$Tx1.C2.Time[i])
  }
}
Tx1_times <- c(C1_time_to_death, C2_time_to_death)
hist(Tx1_times)
# I am not convinced there is going to be a good distribution that fits Tx1_times. 
# The exponential is worth a try, it is often used for time-to-event fitting.

Tx1_time_exponential <- fitdist(Tx1_times, distr="exp")
plot(Tx1_time_exponential)
gofstat(Tx1_time_exponential, fitnames="Exponential distribution for time to death")

#Nope, it's miles off.

Tx1_time_gamma <- fitdist(Tx1_times, distr="gamma")
plot(Tx1_time_gamma)
Tx1_time_weibull <- fitdist(Tx1_times, distr="weibull")
plot(Tx1_time_weibull)
gofstat(list(Tx1_time_exponential, Tx1_time_gamma, Tx1_time_weibull), fitnames=c(
  "Exponential", "Gamma", "Weibull"
))

#Weibull is the best fit, but I would personally just use the mean.

#If the patient has a major complication they will discontinue treatment. So we will only look at C1.
Tx1_time_to_major <- vector()

for (i in 1:length(C1_events)) {
  if (C1_events[i] == 2) {
    Tx1_time_to_major <- c(Tx1_time_to_major, data$Tx1.C1.Time[i])
  }
}

Tx1_major_gamma <- fitdist(Tx1_time_to_major, distr="gamma")
Tx1_major_weibull <- fitdist(Tx1_time_to_major, distr="weibull")
plot(Tx1_major_gamma)
plot(Tx1_major_weibull)
gofstat(list(Tx1_major_gamma, Tx1_major_weibull), fitnames=c("gamma", "weibull"))

#For time to major, the gamma distribution seems to be a slightly better fit.

# I am going to assume that Tx2 is symmetrical with Tx1 and we can use the same values for both. 


#startified porbabilities
major_poor_response=0
major_good_response=0
major_poor_nonresponse=0
major_good_nonresponse=0
minor_poor_response=0
minor_good_response=0
minor_poor_nonresponse=0
minor_good_nonresponse=0

i=0
for(i in 1:length(data$Poor)){
  if(data$Poor[i] == 1) {
    if (data$Tx1.C1.Dx.Pet[i]==1){
      if (data$Tx1.C1.Event[i] == 2) {
        major_poor_response=major_poor_response+1
      }
      else if (data$Tx1.C1.Event[i] == 3) {
        minor_poor_response=minor_poor_response+1
      }
    }
    else if (data$Tx1.C1.Dx.Pet[i]==0){
      if (data$Tx1.C1.Event[i] == 2) {
        major_poor_nonresponse=major_poor_nonresponse+1
      }
      else if (data$Tx1.C1.Event[i] == 3) {
        minor_poor_nonresponse=minor_poor_nonresponse+1
      }
    }
  }
  else if (data$Poor[i] == 0) {
    if (data$Tx1.C1.Dx.Pet[i]==1){
      if (data$Tx1.C1.Event[i] == 2) {
        major_good_response=major_good_response+1
      }
      else if (data$Tx1.C1.Event[i] == 3) {
        minor_good_response=minor_good_response+1
      }
    }
    else if (data$Tx1.C1.Dx.Pet[i]==0){
      if (data$Tx1.C1.Event[i] == 2) {
        major_good_nonresponse=major_good_nonresponse+1
      }
      else if (data$Tx1.C1.Event[i] == 3) {
        minor_good_nonresponse=minor_good_nonresponse+1
      }
    }
  }
}

prob_major_good_response=(major_good_response/sum(data$Poor==0 & data$Tx1.C1.Dx.Pet==1))
prob_minor_good_response=(minor_good_response/sum(data$Poor==0 & data$Tx1.C1.Dx.Pet==1))

prob_major_good_nonresponse=(major_good_nonresponse/sum(data$Poor==0 & data$Tx1.C1.Dx.Pet==0))
prob_minor_good_nonresponse=(minor_good_nonresponse/sum(data$Poor==0 & data$Tx1.C1.Dx.Pet==0))

prob_major_poor_response=(major_poor_response/sum(data$Poor==1 & data$Tx1.C1.Dx.Pet==1))
prob_minor_poor_response=(minor_poor_response/sum(data$Poor==1 & data$Tx1.C1.Dx.Pet==1))

prob_major_poor_nonresponse=(major_poor_nonresponse/sum(data$Poor==1 & data$Tx1.C1.Dx.Pet==0))
prob_minor_poor_nonresponse=(minor_poor_nonresponse/sum(data$Poor==1 & data$Tx1.C1.Dx.Pet==0))