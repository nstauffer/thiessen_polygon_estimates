#### Get the functions loaded ####
source("C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/functions.R")

##### CONFIGURATION ####
n_sims <- 25
sim_seed_offset <- 57

# Find the sim files
sim_file <- paste0("C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/code/",
                   "workflow_multisample_oneframe_intensification_onethiessen_categorical.R")

# Simluation
projection <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
output_path <- paste0("C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/simulations/",
                      "intensification_categorical_increasedthiessen",
                      "/output/")

# Raster
raster_type <- "categorical"
raster_ncol <- 1000
raster_nrow <- 1000
raster_resolution = 1
raster_autocorr_range = 50
raster_mag_var = 10
raster_nug = 0.2
raster_mean = 1
raster_rescale = TRUE
raster_n_categories <- 3

# AOI
aoi_n_vertices <- 6
aoi_convex_hull <- TRUE

# Sample
frame_n_vertices <- 6
frame_convex_hull <- TRUE
sample_type <- "simple"
n_sample_points <- 25
sample_seeds <- 1:99

# Thiessen polygons
thiessen_distribution <- "simple"
thiessen_n_polygons <- 13
thiessen_minimum_sample <- 2

# Analysis
analysis_alpha <- 0.05

#### SIMULATE ####
# Run the sims
for (simulation_seed in (1 + sim_seed_offset):(n_sims + sim_seed_offset)) {
  raster_seed <- 420 * simulation_seed
  aoi_seed <- 1123 * simulation_seed
  frame_seed <- 111 * simulation_seed
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
