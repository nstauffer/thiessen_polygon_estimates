library(ggplot2)
#### Get the functions loaded ####
source("C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/functions.R")

# The path to the place where all the simulations are living
simulations_path <- "C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/simulations"

# The names of the folders in simulations_path with the relevant simulations
# simulations <- c("intensification_continuous",
#                  "intensification_continuous_midthiessen",
#                  "intensification_continuous_increasedthiessen")
simulations <- c("intensification_categorical",
                 "intensification_categorical_increasedthiessen")

# Read in the results from all the simulations!
results <- read_results(simulations = simulations,
                        simulations_path = simulations_path)


# Get the tolerance info
summary_5percent <- tolerance_summary(data = results,
                                      variable = "proportion",
                                      comparison_variable = "proportion_true",
                                      grouping_variables = c("sim_name", "weighting"),
                                      percent_tolerance = 5)
summary_10percent <- tolerance_summary(data = results,
                                       variable = "proportion",
                                       comparison_variable = "proportion_true",
                                       grouping_variables = c("sim_name", "weighting"),
                                       percent_tolerance = 10)
summary <- rbind(summary_5percent,
                 summary_10percent)

# PLOT! And then stare in confusion
ggplot(summary) +
  geom_point(aes(x = sim_name,
                 y = proportion_within_tolerance,
                 color = weighting),
             size = 4,
             shape = "diamond") +
  ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45)) +
  facet_wrap(~percent_tolerance,
             labeller = label_context)
