# R Script to Analyse PSA Results Data

library("ggplot2")

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

ICERs_vector <- vector()
comparison_vector <- vector()
comparison_vector_80000 <- vector()

for (i in 1:length(psa.results$c.bsc)) {
  dif_effects <- e.exp[i] - e.bsc[i]
  dif_costs <- c.exp[i] - c.bsc[i]
  comparison <- dif_costs < (dif_effects * 20000)
  comparison_80000 <- dif_costs < (dif_effects * 80000)
  ICERs_vector <- c(ICERs_vector, (dif_costs / dif_effects))
  comparison_vector <- c(comparison_vector, comparison)
  comparison_vector_80000 <- c(comparison_vector_80000, comparison_80000)
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
length <- length(ICERs_vector)
percentage_below_20000 <- sum(ICERs_vector<20000) / length
percentage_below_80000 <- sum(ICERs_vector<80000) / length

# I think calculating individuals ICERs is incorrect.
# Instead, let's compare the y value (cost) against x*WTP to see if it's lower.

percentage_below_20000 <- sum(comparison_vector==TRUE) / length(comparison_vector)
percentage_below_80000 <- sum(comparison_vector_80000==TRUE) / length(comparison_vector_80000)

# This seems wrong. 
# The south-east quadrant should have negative ICERs (you are saving money while gaining utility)
# However the north-west quadrant should have large ICERs because the new treatment is costing you money 
# while also resulting in harm.


# Let's calculate the fraction for 200 WTPs: from 0 to 200,000€
percentages_vector <- vector()
for (i in seq(0, 200000, length.out = 201)) {
  percentage <- sum(ICERs_vector<i) / length
  percentages_vector <- c(percentages_vector, percentage)
}

plot(x=seq(0, 200000, length.out = 201), y=percentages_vector, xlab="Willingness to Pay Threshold",
  ylab="Points below WTP curve", main="CEAC curve")
