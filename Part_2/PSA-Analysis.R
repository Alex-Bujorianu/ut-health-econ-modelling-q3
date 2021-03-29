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

for (i in 1:length(psa.results$c.bsc)) {
  dif_effects <- e.exp[i] - e.bsc[i]
  dif_costs <- c.exp[i] - c.bsc[i]
  effects_vector <- c(effects_vector, dif_effects)
  costs_vector <- c(costs_vector, dif_costs)
}

plane <- data.frame(effects_vector, costs_vector)
plot(plane, xlab="Incremental effects", ylab="Incremental costs")

##ggplot
ggplot(plane, aes("Inc costs", "Inc effects")) +
        xlim(-0.1, 0.15) +
        ylim(-4000, 4000) +
         geom_point(data=plane) +
        geom_abline(intercept = 0, slope = 20000)