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