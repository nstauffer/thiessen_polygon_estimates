##### SIMULATE ####
n_sims <- 100
sim_seed_offset <- 0

# Find the sim files
sim_files <- list.files(path = "C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/simulations/intensification_categorical_decreasedthiessen/code/",
                        pattern = "\\.[rR]",
                        full.names = TRUE)

# Run the sims
for (simulation_seed in (1 + sim_seed_offset):(n_sims + sim_seed_offset)) {
  source(sim_files)
}

beepr::beep(5)

#### READ IN RESULTS ####
# Find the results files
results_files <- list.files(path = "C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/simulations/intensification_categorical_decreasedthiessen/output/results",
                            pattern = "\\.(csv|CSV)$",
                            full.names = TRUE)

# Read in the results
results_list <- lapply(X = results_files,
                       FUN = function(X){
                         # Get a unique ID from the filename's suffix (output_id_suffix in the sim script)
                         uid <- stringr::str_extract(X,
                                                     pattern = "_\\d+\\.(csv|CSV)$")
                         uid <- stringr::str_extract(uid,
                                                     pattern = "\\d+")
                         
                         # Read in the results
                         current_results <- read.csv(X,
                                                     stringsAsFactors = FALSE)
                         
                         # Add the unique ID to the results
                         current_results$sim_id <- uid
                         
                         current_results
                       })

# Combine all the results!
full_results <- do.call(rbind,
                        results_list)

# Why is this not a string???
full_results$category <- paste(full_results$category)


raster_summary <- unique(full_results[, c("sim_id", "raster_id", "aoi_id", "category", "proportion_true")])

#### PLOT ####
library(ggplot2)
# Boxplot of means
ggplot() +
  geom_boxplot(data = full_results,
               aes(x = category,
                   y = proportion)) +
  geom_point(data = raster_summary,
             aes(x = category,
                 y = proportion_true),
             color = "red") +
  facet_wrap(sim_id~weighted,
             ncol = 3)


#### TEST ####
# This is testing within each simulation run if the calculated proportions are distinguishable from the true proportions
# Start off by splitting by sim run and weighting approach
results_list <- split(full_results,
                      list(full_results$sim_id, full_results$weighted))
# For each weighting approach in each simulation run, run the wilcoxon ranked sign test for each category
wilcoxon_results_list <- lapply(X = results_list,
                                FUN = function(X){
                                  results_by_category <- split(X,
                                                               X$category)
                                  wilcoxon_results_list <- lapply(X = results_by_category,
                                                                  FUN = function(X){
                                                                    wilcoxon_results <- wilcox.test(x = X$proportion,
                                                                                                    y = X$proportion_true,
                                                                                                    paired = TRUE)
                                                                    data.frame(category = X$category[1],
                                                                               p_value = wilcoxon_results$p.value)
                                                                  })
                                  wilcoxon_results <- do.call(rbind,
                                                              wilcoxon_results_list)
                                  wilcoxon_results$weighting <- X$weighted[1]
                                  wilcoxon_results$sim_id <- X$sim_id[1]
                                  wilcoxon_results
                                })
# Combine results
wilcoxon_results <- do.call(rbind,
                            wilcoxon_results_list)
# Identify the ones that were "distinguishable", that is, have a p value > our alpha (0.05)
wilcoxon_results$indistinguishable <- wilcoxon_results$p_value > 0.05

# Summarize the counts of "accurate" versus "inaccurate" sims
wilcoxon_results_summary <- as.data.frame(cbind(table(wilcoxon_results$weighting[wilcoxon_results$indistinguishable]),
                                                table(wilcoxon_results$weighting[!wilcoxon_results$indistinguishable])))
names(wilcoxon_results_summary) <- c("n_indistinguishable_prediction", "n_distinguishable_prediction")
wilcoxon_results_summary$weighting <- names(table(wilcoxon_results$weighting[wilcoxon_results$accurate]))
# Add in a proportion that were indistinguishable from the true proportion
wilcoxon_results_summary$proportion_indistinguishable <- apply(X = wilcoxon_results_summary,
                                                               MARGIN = 1,
                                                               FUN = function(X){
                                                                 total <- X[["n_indistinguishable_prediction"]] +  X[["n_distinguishable_prediction"]]
                                                                 X[["n_indistinguishable_prediction"]] /  total
                                                               })
