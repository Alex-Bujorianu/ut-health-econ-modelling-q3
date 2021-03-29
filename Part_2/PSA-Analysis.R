# R Script to Analyse PSA Results Data

library("ggplot2")
library("rlist")

load("Data/psa_results.RData")

# Calculate ICER. Mean difference in costs / mean difference in effects

mean_effects_bsc <- mean(psa.results$e.bsc)
mean_effects_exp <- mean(psa.results$e.exp)
mean_costs_bsc <- mean(psa.results$c.bsc)
mean_costs_exp <- mean(psa.results$c.exp)

ICER <- (mean_costs_exp - mean_costs_bsc) / (mean_effects_exp - mean_effects_bsc)

# Plot the incremental cost effectiveness plane

# For each patient, find the difference in costs and effects

e.exp <- psa.results$e.exp
e.bsc <- psa.results$e.bsc
c.exp <- psa.results$c.exp
c.bsc <- psa.results$c.bsc
effects_vector <- vector()
costs_vector <- vector()

#Time for some real programming.

list_of_costs_effects <- list()

for (i in 1:length(psa.results$c.bsc)) {
  dif_effects <- e.exp[i] - e.bsc[i]
  dif_costs <- c.exp[i] - c.bsc[i]
  A <- c(dif_costs, dif_effects)
  list_of_costs_effects <- list.append(list_of_costs_effects, A)
  effects_vector <- c(effects_vector, dif_effects)
  costs_vector <- c(costs_vector, dif_costs)
}

plane <- data.frame(effects_vector, costs_vector)
plot(plane, xlab="Incremental effects", ylab="Incremental costs")

##ggplot
WTP_plot <- ggplot(plane, aes(x=effects_vector, y=costs_vector)) +
  scale_x_continuous(limits = c(-0.1, 0.15)) + scale_y_continuous(limits = c(-4000, 4000)) +
  geom_point(data=plane) +
  geom_abline(intercept = 0, slope = 20000) +
  geom_abline(intercept = 0, slope = 80000) +
  xlab("Incremental effects") +
  ylab("Incremental costs") # for the y axis label

WTP_plot
        
##We have to determine how many fall below the threshold. We do this by calculating the individual ICERs.

# I think calculating individuals ICERs is incorrect.
# Instead, let's compare the y value (cost) against x*WTP to see if it's lower.
# They are indeed much more plausible.
WTP_20K_vector <- vector()
WTP_80K_vector <- vector()
WTP_20K <- 20000
WTP_80K <- 80000
for (i in 1:length(list_of_costs_effects)) {
  WTP_20K_vector <- c(WTP_20K_vector, list_of_costs_effects[[i]][1] < list_of_costs_effects[[i]][2]*WTP_20K)
  WTP_80K_vector <- c(WTP_80K_vector, list_of_costs_effects[[i]][1] < list_of_costs_effects[[i]][2]*WTP_80K)
}
percentage_below_20000 <- sum(WTP_20K_vector==TRUE) / length(WTP_20K_vector)
percentage_below_80000 <- sum(WTP_80K_vector==TRUE) / length(WTP_80K_vector)

# This seems wrong. 
# The south-east quadrant should have negative ICERs (you are saving money while gaining utility)
# However the north-west quadrant should have large ICERs because the new treatment is costing you money 
# while also resulting in harm.


# Let's calculate the fraction for 200 WTPs: from 0 to 200,000€
percentages_vector <- vector()
for (i in seq(0, 200000, length.out = 201)) {
  temp_vector <- vector()
  for (j in 1:length(list_of_costs_effects)) {
    temp_vector <- c(temp_vector, list_of_costs_effects[[j]][1] < list_of_costs_effects[[j]][2]*i)
  }
  percentage <- sum(temp_vector==TRUE) / length(temp_vector)
  percentages_vector <- c(percentages_vector, percentage)
}

plot(x=seq(0, 200000, length.out = 201), y=percentages_vector, xlab="Willingness to Pay Threshold",
  ylab="Points below WTP curve", main="CEAC curve")
