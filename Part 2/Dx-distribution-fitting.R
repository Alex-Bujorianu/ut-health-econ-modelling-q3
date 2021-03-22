#Script to determine distribution of Dx1, Dx2 and Dx3 from the dataset

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
