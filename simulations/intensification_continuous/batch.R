##### SIMULATE ####
n_sims <- 4
sim_seed_offset <- 69

# Find the sim files
sim_files <- list.files(path = "C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/simulations/intensification_continuous/code/",
                        pattern = "\\.[rR]",
                        full.names = TRUE)

# Run the sims
for (simulation_seed in (1 + sim_seed_offset):(n_sims + sim_seed_offset)) {
  source(sim_files)
}

beepr::beep(5)

#### READ IN RESULTS ####
# Find the results files
results_files <- list.files(path = "C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/simulations/intensification_continuous/output/results",
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
               aes(x = weighted,
                   y = mean)) +
  geom_hline(data = raster_summary,
             aes(yintercept = mean_true),
             color = "red") +
  facet_wrap(~sim_id)

# Boxplot of variance
ggplot() +
  geom_boxplot(data = full_results,
               aes(x = weighted,
                   y = variance)) +
  geom_hline(data = raster_summary,
             aes(yintercept = variance_true),
             color = "red") +
  facet_wrap(~sim_id)


#### TEST ####
