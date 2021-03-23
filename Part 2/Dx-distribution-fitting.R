#Script to determine distribution of Dx1, Dx2 and Dx3 from the dataset

library("fitdistrplus")
library("SimDesign")

load("trial_dataset.RData")

#I am only going to fit the data on C1 of Tx1, the first cycle
#This is what it says on the PDF
#I will distinguish between responders and non-responders
#I am going to combine the data and fit a multivariate normal distribution
#I would rather just fit separate distributions for each of the 3 tests, but whatever

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

quantile(c(Test1_responders, Test2_responders, Test3_responders), c(0.8, 0.8, 0.8))
quantile(Test1_nonresponders, 0.2)
#0.71 is the exact value by which they diverge
test1_boundary <- 0.71

quantile(Test2_responders, 0.8)
quantile(Test2_nonresponders, 0.2)
#0 for the 2nd test
test2_boundary <- 0

quantile(Test3_responders, 0.8)
quantile(Test3_nonresponders, 0.2)
#They overlap more
test3_boundary <- mean(c(quantile(Test3_responders, 0.8),
     quantile(Test3_nonresponders, 0.2)))