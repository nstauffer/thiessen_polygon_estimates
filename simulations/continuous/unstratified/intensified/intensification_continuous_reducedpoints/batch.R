#### Get the functions loaded ------------------------------------------------
source("C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/functions.R")

#### CONFIGURATION ------------------------------------------------
n_sims <- 50
sim_seed_offset <- 99

sim_file <- "C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/code/workflow_continuous_intensification.R"

# Simulation
projection <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
output_path <- paste0("C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/simulations/",
                      "intensification_continuous_reducedpoints",
                      "/output/")

# Raster
raster_type <- "continuous"
raster_ncol <- 1000
raster_nrow <- 1000
raster_resolution = 1
raster_autocorr_range = 50
raster_mag_var = 10
raster_nug = 0.2
raster_mean = 1
raster_rescale = TRUE


# AOI
aoi_n_vertices <- 6
aoi_convex_hull <- TRUE


# Sample
frame_n_vertices <- 6
frame_convex_hull <- TRUE
sample_type <- "simple"
n_sample_points <- 10
sample_seeds <- 1:99

# Thiessen polygons
thiessen_distribution <- "simple"
thiessen_n_polygons <- 5
thiessen_minimum_sample <- 2

# Analysis
analysis_alpha <- 0.05
percent_tolerance <- 5

#### SIMULATE ------------------------------------------------
# Run the sims
for (simulation_seed in (1 + sim_seed_offset):(n_sims + sim_seed_offset)) {
  raster_seed <- 420 * simulation_seed
  aoi_seed <- 1123 * simulation_seed
  frame_seed <- 111 * simulation_seed
  source(sim_file)
}

beepr::beep(5)

#### READ IN RESULTS ------------------------------------------------
# Find the results files
results_files <- list.files(path = paste0(output_path, "/results"),
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
                         
                         # Just a quick correction if necessary
                         if ("weighted" %in% names(current_results)) {
                           names(current_results)[names(current_results) == "weighted"] <- "weighting"
                         }
                         
                         current_results
                       })

# Combine all the results!
full_results <- do.call(rbind,
                        results_list)

raster_summary <- unique(full_results[, c("sim_id", "raster_id", "aoi_id", "mean_true", "sd_true", "variance_true")])

#### PLOT ------------------------------------------------
library(ggplot2)
# Boxplot of means
ggplot() +
  geom_boxplot(data = full_results,
               aes(x = weighting,
                   y = mean)) +
  geom_hline(data = raster_summary,
             aes(yintercept = mean_true),
             color = "red") +
  facet_wrap(~sim_id)

# Boxplot of variance
ggplot() +
  geom_boxplot(data = full_results,
               aes(x = weighting,
                   y = variance)) +
  geom_hline(data = raster_summary,
             aes(yintercept = variance_true),
             color = "red") +
  facet_wrap(~sim_id)


#### TEST ------------------------------------------------
# This is testing within each simulation run if the calculated proportions are distinguishable from the true proportions
# Start off by splitting by sim run and weighting approach
results_list <- split(full_results,
                      list(full_results$sim_id, full_results$weighted))
# For each weighting approach in each simulation run, run the wilcoxon ranked sign test for each category
wilcoxon_results_list <- lapply(X = results_list,
                                FUN = function(X){
                                  
                                  wilcoxon_results <- wilcox.test(x = X$mean,
                                                                  y = X$mean_true,
                                                                  paired = TRUE)
                                  wilcoxon_results <- data.frame(weighting = X$weighted[1],
                                                                 sim_id = X$sim_id[1],
                                                                 p_value = wilcoxon_results$p.value)
                                  wilcoxon_results
                                })
# Combine results
wilcoxon_results <- do.call(rbind,
                            wilcoxon_results_list)
# Identify the ones that were "distinguishable", that is, have a p value > our alpha (0.05)
wilcoxon_results$indistinguishable <- wilcoxon_results$p_value > analysis_alpha

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

# Okay, but what about setting some kind of difference threshold????
# We'll test to see how the measured mean compares to the true mean by seeing if
# it falls within ±x% of the true value.
# That isn't true value in percent ± x% but true value in percent ± x% of the true value
threshold_results <- tolerance_summary(data = full_results,
                                       variable = "mean",
                                       comparison_variable = "mean_true",
                                       grouping_variables = c("weighting"),
                                       percent_tolerance = percent_tolerance)


# Plot the results to look for weird distributions
ggplot(full_results) +
  geom_point(aes(x = mean,
                 y = mean_true,
                 color = sim_id),
             alpha = 0.1) +
  geom_abline(slope = 1,
              intercept = 0) +
  facet_wrap(~weighting)
