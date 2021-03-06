# Step 1.4 of assignment and some additional distribution fitting

load("Data/trial_dataset.RData");

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
mean_male_ages <- mean(male_ages)
sd_male_ages <- sd(male_ages)
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
C2_transformed_events <- data$Tx1.C2.Event
#We need to transform this data into alive or dead coded 0 or 1.
for (i in 1:length(C2_transformed_events)) {
  #If the value is 1, the patient died died in Tx1.C2
  #I have to test for NA values to make R evaluate the if condition properly
  if (!is.na(C2_transformed_events[i]) && (C2_transformed_events[i] == 1)) {
    C2_transformed_events[i] <- 0
  }
  #If the value is NA, the patient either died in cycle 1 OR had a major complication, and stopped treatment.
  #Check if they died in C1 and assign death if that is the case
  else if (is.na(C2_transformed_events[i]) && data$Tx1.C1.Event[i]==1) {
    C2_transformed_events[i] <- 0
  }
  else {
    C2_transformed_events[i] <- 1
  }
}

logreg_death <- glm(C2_transformed_events ~ male + age, family = binomial)
summary(logreg_death)

#Age is not significant. Male gender seems to have stronger significance, 
#but still it does not reach significance at the alpha = 0.10 level.
#We had 538 observations BUT only 31 people died. 
#It's hard to get statistical power for an unlikely event, 
#e.g. vaccine complications are hard to prove statistically: look at Guillian-Barre for the bird flu vaccine.

#We want to find time-to-death and time-to-major-complication for Tx1 and Tx2 based on C1 + C2

# Time to death
C1_events <- data$Tx1.C1.Event
C1_time_to_death <- vector()
for (i in 1:length(C1_events)) {
  if (C1_events[i] == 1) { #if death
    C1_time_to_death <- c(C1_time_to_death, data$Tx1.C1.Time[i])
  }
}
C2_events <- data$Tx1.C2.Event
C2_time_to_death <- vector()
for (i in 1:length(C2_events)) {
  if (!is.na(C2_events[i]) && C2_events[i] == 1) {
    #We are only calculating the time they spend in this individual subcycle!
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
Tx1_time_uniform <- fitdist(Tx1_times, distr="unif")
plot(Tx1_time_uniform)
gofstat(list(Tx1_time_exponential, Tx1_time_gamma, Tx1_time_weibull, Tx1_time_uniform), fitnames=c(
  "Exponential", "Gamma", "Weibull", "Uniform"
))

#Turns out uniform is the closest-fitting distribution.

#If the patient has a major complication they will discontinue treatment. So we will only look at C1.
Tx1_time_to_major <- vector()

for (i in 1:length(C1_events)) {
  if (C1_events[i] == 2) {
    Tx1_time_to_major <- c(Tx1_time_to_major, data$Tx1.C1.Time[i])
  }
}

Tx1_major_gamma <- fitdist(Tx1_time_to_major, distr="gamma")
Tx1_major_weibull <- fitdist(Tx1_time_to_major, distr="weibull")
Tx1_major_exp <- fitdist(Tx1_time_to_major, distr="exp")
plot(Tx1_major_gamma)
plot(Tx1_major_weibull)
plot(Tx1_major_exp) #exponential is miles off
gofstat(list(Tx1_major_gamma, Tx1_major_weibull), fitnames=c("gamma", "weibull"))

#For time to major, the gamma distribution seems to be a slightly better fit.

# I am going to assume that Tx2 is symmetrical with Tx1 and we can use the same values for both. 


#Calculating stratified probabilities of major and minor complications
#intialisations
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
  if(data$Poor[i] == 1) {  #dividing the data as per clinical condition
    if (data$Tx1.C1.Dx.Pet[i]==1){ #dividing the data as per response
      if (data$Tx1.C1.Event[i] == 2) {
        major_poor_response=major_poor_response+1 #Number of major complications in this subset
      }
      else if (data$Tx1.C1.Event[i] == 3) {
        minor_poor_response=minor_poor_response+1 #Number of minor complications in this subset
      }
    }
    else if (data$Tx1.C1.Dx.Pet[i]==0){
      if (data$Tx1.C1.Event[i] == 2) {
        major_poor_nonresponse=major_poor_nonresponse+1 #Number of major complications in this subset
      }
      else if (data$Tx1.C1.Event[i] == 3) {
        minor_poor_nonresponse=minor_poor_nonresponse+1 #Number of minor complications in this subset
      }
    }
  }
  else if (data$Poor[i] == 0) {
    if (data$Tx1.C1.Dx.Pet[i]==1){
      if (data$Tx1.C1.Event[i] == 2) {
        major_good_response=major_good_response+1 #Number of major complications in this subset
      }
      else if (data$Tx1.C1.Event[i] == 3) {
        minor_good_response=minor_good_response+1 #Number of minor complications in this subset
      }
    }
    else if (data$Tx1.C1.Dx.Pet[i]==0){
      if (data$Tx1.C1.Event[i] == 2) {
        major_good_nonresponse=major_good_nonresponse+1 #Number of major complications in this subset
      }
      else if (data$Tx1.C1.Event[i] == 3) {
        minor_good_nonresponse=minor_good_nonresponse+1 #Number of minor complications in this subset
      }
    }
  }
}

#Final probabilities
prob_major_good_response=(major_good_response/sum(data$Poor==0 & data$Tx1.C1.Dx.Pet==1))
prob_minor_good_response=(minor_good_response/sum(data$Poor==0 & data$Tx1.C1.Dx.Pet==1))

prob_major_good_nonresponse=(major_good_nonresponse/sum(data$Poor==0 & data$Tx1.C1.Dx.Pet==0))
prob_minor_good_nonresponse=(minor_good_nonresponse/sum(data$Poor==0 & data$Tx1.C1.Dx.Pet==0))

prob_major_poor_response=(major_poor_response/sum(data$Poor==1 & data$Tx1.C1.Dx.Pet==1))
prob_minor_poor_response=(minor_poor_response/sum(data$Poor==1 & data$Tx1.C1.Dx.Pet==1))

prob_major_poor_nonresponse=(major_poor_nonresponse/sum(data$Poor==1 & data$Tx1.C1.Dx.Pet==0))
prob_minor_poor_nonresponse=(minor_poor_nonresponse/sum(data$Poor==1 & data$Tx1.C1.Dx.Pet==0))

# Beta distributions for minor complications
n_good_response <- sum(data$Poor==0 & data$Tx1.C1.Dx.Pet==1)
n_good_nonresponse <- sum(data$Poor==0 & data$Tx1.C1.Dx.Pet==0)
n_poor_response <- sum(data$Poor==1 & data$Tx1.C1.Dx.Pet==1)
n_poor_nonresponse <- sum(data$Poor==1 & data$Tx1.C1.Dx.Pet==0)
beta_minor_good_response <- rbeta(1, minor_good_response, n_good_response - minor_good_response)
beta_minor_good_nonresponse <- rbeta(1, minor_good_nonresponse, n_good_nonresponse - minor_good_nonresponse)
beta_minor_poor_response <- rbeta(1, minor_poor_response, n_poor_response - minor_poor_response)
beta_minor_poor_nonresponse <- rbeta(1, minor_poor_nonresponse, n_poor_nonresponse - minor_poor_nonresponse)

#Stratify the patients’ utility during Tx1 treatment according to their treatment response, based on C1 of Tx1
Tx1_utility_responders <- vector()
Tx1_utility_nonresponders <- vector()
i=0
for (i in 1:length(data$Tx1.C1.Dx.Pet)) {
  if (data$Tx1.C1.Dx.Pet[i] == 1) {
    Tx1_utility_responders <- c(Tx1_utility_responders,data$Tx1.C1.QoL[i])
  }
  else if (data$Tx1.C1.Dx.Pet[i]==0) {
    Tx1_utility_nonresponders <- c(Tx1_utility_nonresponders, data$Tx1.C1.QoL[i])
  }
}

#Stratify the patients’ utility during Tx2 treatment according to their treatment response, based on C1 of Tx2 (do not stratify the disutility)
Tx2_utility_responders <- vector()
Tx2_utility_nonresponders <- vector()
i=0
for (i in 1:length(data$Tx2.C1.Dx.Pet)) {
  if (!is.na(data$Tx2.C1.Dx.Pet[i]) && data$Tx2.C1.Dx.Pet[i] == 1) {
    Tx2_utility_responders <-c(Tx2_utility_responders,data$Tx2.C1.QoL[i])
  }
  else if (!is.na(data$Tx2.C1.Dx.Pet[i]) && data$Tx2.C1.Dx.Pet[i]==0) {
    Tx2_utility_nonresponders <- c(Tx2_utility_nonresponders,data$Tx2.C1.QoL[i])
  }
}

#Step 2.7.3
hist(Tx1_utility_responders)
dist_lognorm_Tx1_u_r <- fitdist(Tx1_utility_responders, distr = "lnorm")
#Check
gofstat(dist_lognorm_Tx1_u_r, fitnames=c("Tx1 utility responders lognormal"))
plot(dist_lognorm_Tx1_u_r)
#function to determine utlity of patients depending on response
func.tx1_u_r <- function() {
  mean_tx1_u_r <- mean(Tx1_utility_responders); # mean utility of patients responding to treatment in tx1
  sd_tx1_u_r <- sd(Tx1_utility_responders); # SD of utility of patients responding to treatment in tx1
  location <- log(mean_tx1_u_r^2 / sqrt(sd_tx1_u_r^2 + mean_tx1_u_r^2))
  shape <- sqrt(log(1 + (sd_tx1_u_r^2 / mean_tx1_u_r^2)))
  tx1_u_r <- rlnorm(n=1, meanlog=location,sdlog=shape)
  return(tx1_u_r)
}

hist(Tx1_utility_nonresponders)
dist_lognorm_Tx1_u_nr <- fitdist(Tx1_utility_nonresponders, distr = "lnorm")
#Check
gofstat(dist_lognorm_Tx1_u_nr, fitnames=c("Tx1 utility nonresponders lognormal"))
plot(dist_lognorm_Tx1_u_nr)
#function to determine utlity of patients depending on response
func.tx1_u_nr <- function() {
  mean_tx1_u_nr <- mean(Tx1_utility_nonresponders); # mean utility of patients not responding to treatment in tx1
  sd_tx1_u_nr <- sd(Tx1_utility_nonresponders); # SD of utility of patients not responding to treatment in tx1
  location <- log(mean_tx1_u_nr^2 / sqrt(sd_tx1_u_nr^2 + mean_tx1_u_nr^2))
  shape <- sqrt(log(1 + (sd_tx1_u_nr^2 / mean_tx1_u_nr^2)))
  tx1_u_nr <- rlnorm(n=1, meanlog=location,sdlog=shape)
  return(tx1_u_nr)
}

hist(Tx2_utility_responders)
dist_lognorm_Tx2_u_r <- fitdist(Tx1_utility_responders, distr = "lnorm")
#Check
gofstat(dist_lognorm_Tx2_u_r, fitnames=c("Tx2 utility responders lognormal"))
plot(dist_lognorm_Tx2_u_r)
#function to determine utlity of patients depending on response
func.tx2_u_r <- function() {
  mean_tx2_u_r <- mean(Tx2_utility_responders); # mean utility of patients responding to treatment in tx2
  sd_tx2_u_r <- sd(Tx2_utility_responders); # SD of utility of patients responding to treatment in tx2
  location <- log(mean_tx2_u_r^2 / sqrt(sd_tx2_u_r^2 + mean_tx2_u_r^2))
  shape <- sqrt(log(1 + (sd_tx2_u_r^2 / mean_tx2_u_r^2)))
  tx2_u_r <- rlnorm(n=1, meanlog=location,sdlog=shape)
  return(tx2_u_r)
}

hist(Tx2_utility_nonresponders)
dist_lognorm_Tx2_u_nr <- fitdist(Tx2_utility_nonresponders, distr = "lnorm")
#Check
gofstat(dist_lognorm_Tx2_u_nr, fitnames=c("Tx2 utility nonresponders lognormal"))
plot(dist_lognorm_Tx2_u_nr)
#function to determine utlity of patients depending on response
func.tx2_u_nr <- function() {
  mean_tx2_u_nr <- mean(Tx2_utility_nonresponders); # mean utility of patients not responding to treatment in tx2
  sd_tx2_u_nr <- sd(Tx2_utility_nonresponders); # SD of utility of patients not responding to treatment in tx2
  location <- log(mean_tx2_u_nr^2 / sqrt(sd_tx2_u_nr^2 + mean_tx2_u_nr^2))
  shape <- sqrt(log(1 + (sd_tx2_u_nr^2 / mean_tx2_u_nr^2)))
  tx2_u_nr <- rlnorm(n=1, meanlog=location,sdlog=shape)
  return(tx2_u_nr)
}


#2.7.6
# Distribution for prob. of death in first follow up
fu <- na.omit(data$FU1.Event)
n_fu <- length(na.omit(data$FU1.Event)) #n persons 
r_fu <- sum(fu==1) #1 is death
beta_death_followup <- rbeta(1, r_fu, n_fu - r_fu)

