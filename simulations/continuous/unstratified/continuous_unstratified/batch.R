#### Get the functions loaded ####
source("C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/functions.R")

#### CONFIGURATION ####
n_sims <- 50
sim_seed_offset <- 99

sim_file <- "C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/code/workflow_continuous_SIMPLIFIED.R"

# Simulation
projection <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
output_path <- paste0("C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/simulations/",
                      "continuous/",
                      "unstratified/",
                      "continuous_unstratified",
                      "/output/")

# Raster
raster_type <- "continuous"
raster_ncol <- 100
raster_nrow <- 100
raster_resolution = 1
raster_autocorr_range = 25
raster_mag_var = 10
raster_nug = 0.2
raster_mean = 1
raster_rescale = TRUE


# AOI
aoi_n_vertices <- 6
aoi_convex_hull <- TRUE


# Sample
sample_type <- "simple"
n_aoi_sample_points <- 10
sample_seeds <- 1:99

# Thiessen polygons
thiessen_distribution <- "simple"
thiessen_n_polygons <- 3
thiessen_minimum_sample <- 2

# Analysis
analysis_alpha <- 0.05
percent_tolerance <- 5

##### SIMULATE ####
# Run the sims
for (simulation_seed in (1 + sim_seed_offset):(n_sims + sim_seed_offset)) {
  raster_seed <- 420 * simulation_seed
  aoi_seed <- 1123 * simulation_seed
  source(sim_file)
}

beepr::beep(5)

#### READ IN RESULTS ####
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
                         
                         current_results
                       })

# Combine all the results!
full_results <- do.call(rbind,
                        results_list)

raster_summary <- unique(full_results[, c("sim_id", "raster_id", "aoi_id", "mean_true", "sd_true", "variance_true")])

#### PLOT ####
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


#### TEST ####
# This is testing within each simulation run if the calculated proportions are distinguishable from the true proportions
# Start off by splitting by sim run and weighting approach
results_list <- split(full_results,
                      list(full_results$sim_id, full_results$weighting))
# For each weighting approach in each simulation run, run the wilcoxon ranked sign test for each category
wilcoxon_results_list <- lapply(X = results_list,
                                FUN = function(X){
                                  
                                  wilcoxon_results <- wilcox.test(x = X$mean,
                                                                  y = X$mean_true,
                                                                  paired = TRUE)
                                  wilcoxon_results <- data.frame(weighting = X$weighting[1],
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
within_tolerance <- tolerance_test(data = full_results,
                       variable = "mean",
                       comparison_variable = "mean_true",
                       percent_tolerance = percent_tolerance)

threshold_var_name <- paste0("within_", percent_tolerance, "_percent")
full_results[[threshold_var_name]] <- within_tolerance

threshold_list <- lapply(X = split(full_results, full_results$weighting),
                         percent_tolerance = percent_tolerance,
                         FUN = function(X, percent_tolerance){
                           # How many sims are we looking at?
                           n_observations <- nrow(X)
                           # Which variable has the info about whether the threshold was met or not?
                           var_name <- paste0("within_", percent_tolerance, "_percent")
                           # How many sims were within the threshold?
                           within_tolerance_count <- sum(X[[var_name]])
                           # What's the proportion within the threshold?
                           proportion_within_tolerance = within_tolerance_count / n_observations
                           # Gimme those results!
                           data.frame(weighting = X[["weighting"]][1],
                                      n_sims = n_observations,
                                      percent_tolerance = percent_tolerance,
                                      n_within_tolerance = within_tolerance_count,
                                      proportion_within_tolerance = proportion_within_tolerance)
                         })

threshold_results <- do.call(rbind,
                             threshold_list)
# Plot the results to look for weird distributions
ggplot(full_results) +
  geom_point(aes(x = mean,
                 y = mean_true,
                 color = sim_id),
             alpha = 0.1) +
  geom_abline(slope = 1,
              intercept = 0) +
  facet_wrap(~weighting)
