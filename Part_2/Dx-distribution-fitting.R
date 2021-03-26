#Script to determine distribution of Dx1, Dx2 and Dx3 from the dataset

library("fitdistrplus")
library("SimDesign")

load("Data/trial_dataset.RData")

#I am only going to fit the data on C1 of Tx1, the first cycle
#This is what it says on the PDF
#I will distinguish between responders and non-responders
#I am going to combine the data and fit a multivariate normal distribution
#This is because the test results are correlated

Test1 <- data$Tx1.C1.Dx.Test1
Test2 <- data$Tx1.C1.Dx.Test2
Test3 <- data$Tx1.C1.Dx.Test3

#Histograms for all 3 tests seem approximately normal
hist(Test1)
hist(Test2)
hist(Test3)

#Histogram for combined vector
hist(c(Test1, Test2, Test3))

#Separate by responders and non-responders

Test1_responders <- vector()
Test1_nonresponders <- vector()

Test2_responders <- vector()
Test2_nonresponders <- vector()

Test3_responders <- vector()
Test3_nonresponders <- vector()

for (i in 1:length(data$Tx1.C1.Dx.Pet)) {
  if (data$Tx1.C1.Dx.Pet[i]==1) { #is responder, assign values to tests
    Test1_responders <- c(Test1_responders, data$Tx1.C1.Dx.Test1[i])
    Test2_responders <- c(Test2_responders, data$Tx1.C1.Dx.Test2[i])
    Test3_responders <- c(Test3_responders, data$Tx1.C1.Dx.Test3[i])
  }
  else {
    Test1_nonresponders <- c(Test1_nonresponders, data$Tx1.C1.Dx.Test1[i])
    Test2_nonresponders <- c(Test2_nonresponders, data$Tx1.C1.Dx.Test2[i])
    Test3_nonresponders <- c(Test3_nonresponders, data$Tx1.C1.Dx.Test3[i])
  }
}

Responders <- c(Test1_responders, Test2_responders, Test3_responders)
Nonresponders <- c(Test1_nonresponders, Test2_nonresponders, Test3_nonresponders)

hist(Responders)
hist(Nonresponders)

#Those are some wacky histograms. 

Test1_responders_norm <- fitdist(Test1_responders, distr="norm")
summary(Test1_responders_norm)
gofstat(Test1_responders_norm)
plot(Test1_responders_norm)

Test2_responders_norm <- fitdist(Test2_responders, distr="norm")
plot(Test2_responders_norm) #OK fit

Test3_responders_norm <- fitdist(Test3_responders, distr="norm")
plot(Test3_responders_norm)

#Indeed a normal distribution is a very good fit for the individual tests

m_vect <- cbind(Test1_responders,
                Test2_responders,
                Test3_responders)

m_cov <- cov(m_vect) # create covariance matrix
v_means <-  c(mean(Test1_responders), mean(Test2_responders), mean(Test3_responders))
print(v_means)
rmvnorm(n = 1, mean = v_means, sigma = m_cov) #Test results for responders

plot(density(Test1_responders), col = 'red', xlim = c(-5, 10), ylim = c(0, 2), main = 'Density functions') # plot the different density functions
lines(density(Test2_responders), col = 'darkgray')
lines(density(Test3_responders), col = 'green')
legend("topright",                    # Add legend to plot
       legend = c("Test1_responders", "Test2_responders", "Test3_responders"),
       col = c('red', 'darkgray', 'green'),
       lty = 1)

non_m_vect <- cbind(Test1_nonresponders,
                    Test2_nonresponders,
                    Test3_nonresponders)
non_m_cov <- cov(non_m_vect)
non_v_means <-  c(mean(Test1_nonresponders), mean(Test2_nonresponders), mean(Test3_nonresponders))
print(non_v_means)

rmvnorm(n = 1, mean = non_v_means, sigma = non_m_cov)

#We are doing a multivariate distribution because the tests are correlated
#MGF theory tells us that if you add 3 normal distributions, you end up with a normal distribution 
#whose mean is the sum of the other 3 means

#Determine threshold values for the 3 tests

quantile(Test1_responders, 0.8)
quantile(Test1_nonresponders, 0.2)
quantile (data$Tx1.C1.Dx.Test1[data$Tx1.C1.Dx.Pet==0], probs=0.80, na.rm=FALSE)
test1_boundary <- quantile (data$Tx1.C1.Dx.Test1[data$Tx1.C1.Dx.Pet==0], probs=0.80, na.rm=FALSE)

quantile (data$Tx1.C1.Dx.Test2[data$Tx1.C1.Dx.Pet==0], probs=0.80, na.rm=FALSE)
quantile(Test2_responders, 0.8)
quantile(Test2_nonresponders, 0.2)
#0 for the 2nd test
test2_boundary <- quantile (data$Tx1.C1.Dx.Test2[data$Tx1.C1.Dx.Pet==0], probs=0.80, na.rm=FALSE)

quantile (data$Tx1.C1.Dx.Test3[data$Tx1.C1.Dx.Pet==0], probs=0.80, na.rm=FALSE)
quantile(Test3_responders, 0.8)
quantile(Test3_nonresponders, 0.2)
#They overlap more
test3_boundary <- quantile (data$Tx1.C1.Dx.Test3[data$Tx1.C1.Dx.Pet==0], probs=0.80, na.rm=FALSE)

validation <- vector()

for (i in 1:length(data$Tx1.C1.Dx.Pet)) {
  #Let's classify them as non-responder if ONE condition is higher than the threshold
  if (data$Tx1.C1.Dx.Test1[i] > test1_boundary || data$Tx1.C1.Dx.Test2[i] > test2_boundary 
      || data$Tx1.C1.Dx.Test3[i] > test3_boundary) {
    validation <- c(validation, 0)
  }
  else {
    validation <- c(validation, 1)
  }
}

hitrate <- vector()
for (i in 1:length(data$Tx1.C1.Dx.Pet)) {
  if (data$Tx1.C1.Dx.Pet[i] == validation[i]) {
    hitrate <- c(hitrate, TRUE)
  }
  else {
    hitrate <- c(hitrate, FALSE)
  }
}


print(sum(hitrate==TRUE) / length(hitrate))

# These boundaries seem to work!

apply( m_vect , 2 , quantile , probs = 0.8 , na.rm = TRUE )
apply( non_m_vect , 2 , quantile , probs = 0.2 , na.rm = TRUE ) #produces same numbers as doing it individually