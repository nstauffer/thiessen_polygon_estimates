#### Get the packages attached ####
library(ggplot2)

#### Get the functions loaded ####
source("C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/functions.R")

#### Set the config ####
# Simluation
n_sims <- 10
projection <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Raster
raster_type <- "continuous"
raster_ncol <- 100
raster_nrow <- 100
# raster_values <- c(1, 2, 3)
# raster_distribution <- "normal"
raster_resolution = 1
raster_autocorr_range = 10
raster_mag_var = 10
raster_nug = 0.2
raster_mean = 1
raster_user_seed = 420
raster_rescale = TRUE
raster_seed <- 420

# AOI
aoi_n_vertices <- 6
aoi_convex_hull <- TRUE
aoi_seed <- 690

# Sample
frame_seed <- 42
frame_n_vertices <- 6
frame_convex_hull <- TRUE
sample_type <- "simple"
n_sample_points <- 50
sample_seeds <- 1:99

# Thiessen polygons
thiessen_distribution <- "simple"
thiessen_proportion <- 0.1
thiessen_seed <- 96
thiessen_minimum_sample <- 2

# Analysis
analysis_alpha <- 0.2

#### Generate a raster ####
current_raster <- switch(raster_type,
                         "categorical" = {
                           landscape_gen_categorical(categories =  raster_values,
                                                     ncol = raster_ncol,
                                                     nrow = raster_nrow,
                                                     seed_number = raster_seed,
                                                     projection = NULL)
                         },
                         "continuous" = {
                           raster <- NLMR::nlm_gaussianfield(ncol = raster_ncol,
                                                             nrow = raster_nrow,
                                                             resolution = raster_resolution,
                                                             autocorr_range = raster_autocorr_range,
                                                             mag_var = raster_mag_var,
                                                             nug = raster_nug,
                                                             mean = raster_mean,
                                                             user_seed = raster_seed,
                                                             rescale = raster_rescale)
                           raster::projection(raster) <- projection
                           raster
                         })

raster_metadata <- data.frame(raster_id = paste0("raster_", raster_seed),
                              raster_seed = raster_seed,
                              raster_type = raster_type,
                              raster_ncol = raster_ncol,
                              raster_nrow = raster_nrow)

#### Generate AOI ####
aoi <- aoi_gen(xmax = raster_ncol,
               xmin = 0,
               ymax = raster_nrow,
               ymin = 0,
               n_vertices = aoi_n_vertices,
               convex_hull = aoi_convex_hull,
               seed_number = aoi_seed)

aoi$raster_id <- raster_metadata[["raster_id"]]

aoi$aoi_id <- paste0(raster_metadata[["raster_id"]],
                     "-",
                     "aoi_",
                     aoi_seed)
aoi$frame_id <- paste0("frame_",
                       aoi_seed)
aoi$frame_seed <- aoi_seed

#### Generate sampling points ####
sample_points_list_1 <- lapply(X = sample_seeds,
                               frame = aoi,
                               sample_type = sample_type,
                               n_points = n_sample_points,
                               projection = projection,
                               raster = current_raster,
                               FUN = function(X,
                                              frame,
                                              sample_type,
                                              n_points,
                                              projection,
                                              raster){
                                 sample_points <- points_gen(frame = frame,
                                                             sample_type = sample_type,
                                                             n_points = n_points,
                                                             seed_number = X,
                                                             projection = projection)
                                 
                                 # Attribute them with raster values
                                 sample_points_spdf <- methods::as(sample_points,
                                                                   "Spatial")
                                 
                                 raster_values <- raster::extract(x = raster,
                                                                  y = sample_points_spdf)
                                 
                                 sample_points$value <- raster_values
                                 
                                 sample_points$frame_id <- unique(frame$frame_id)
                                 
                                 sample_points
                               })

# Create the sampling frame that isn't the AOI
frame_2 <- aoi_gen(xmax = raster_ncol,
                   xmin = 0,
                   ymax = raster_nrow,
                   ymin = 0,
                   n_vertices = frame_n_vertices,
                   convex_hull = frame_convex_hull,
                   seed_number = frame_seed)
frame_2$raster_id <- raster_metadata$raster_id
frame_2$frame_seed <- frame_seed
frame_2$frame_id <- paste0(raster_metadata$raster_id,
                           "-",
                           "frame_",
                           frame_2$frame_seed)

# Generate sampling points for that frame, restricted to those that fall in the AOI
sample_point_list_2 <- lapply(X = sample_seeds,
                              frame = frame_2,
                              aoi = aoi,
                              sample_type = sample_type,
                              n_points = n_sample_points,
                              projection = projection,
                              raster = current_raster,
                              FUN = function(X,
                                             frame,
                                             aoi,
                                             sample_type,
                                             n_points,
                                             projection,
                                             raster){
                                
                                sample_points <- points_gen(frame = frame,
                                                            sample_type = sample_type,
                                                            n_points = n_points,
                                                            seed_number = X,
                                                            projection = projection)
                                
                                sample_point_indices_in_aoi <- as.vector(sf::st_intersects(x = aoi,
                                                                                           y = sample_points,
                                                                                           sparse = FALSE))
                                
                                sample_points <- sample_points[sample_point_indices_in_aoi, ]
                                
                                # Attribute them with raster values
                                sample_points_spdf <- methods::as(sample_points,
                                                                  "Spatial")
                                
                                raster_values <- raster::extract(x = raster,
                                                                 y = sample_points_spdf)
                                
                                sample_points$value <- raster_values
                                
                                sample_points$frame_id <- unique(frame$frame_id)
                                
                                sample_points
                              })

# Create a third sampling frame that isn't the AOI
frame_3 <- aoi_gen(xmax = raster_ncol,
                   xmin = 0,
                   ymax = raster_nrow,
                   ymin = 0,
                   n_vertices = frame_n_vertices,
                   convex_hull = frame_convex_hull,
                   seed_number = frame_seed * 2)
frame_3$raster_id <- raster_metadata$raster_id
frame_3$frame_seed <- frame_seed * 2
frame_3$frame_id <- paste0(raster_metadata$raster_id,
                           "-",
                           "frame_",
                           frame_3$frame_seed)

# Generate sampling points for that frame, restricted to those that fall in the AOI
sample_point_list_3 <- lapply(X = sample_seeds,
                              frame = frame_3,
                              aoi = aoi,
                              sample_type = sample_type,
                              n_points = n_sample_points,
                              projection = projection,
                              raster = current_raster,
                              FUN = function(X,
                                             frame,
                                             aoi,
                                             sample_type,
                                             n_points,
                                             projection,
                                             raster){
                                
                                sample_points <- points_gen(frame = frame,
                                                            sample_type = sample_type,
                                                            n_points = n_points,
                                                            seed_number = X,
                                                            projection = projection)
                                
                                sample_point_indices_in_aoi <- as.vector(sf::st_intersects(x = aoi,
                                                                                           y = sample_points,
                                                                                           sparse = FALSE))
                                
                                sample_points <- sample_points[sample_point_indices_in_aoi, ]
                                
                                # Attribute them with raster values
                                sample_points_spdf <- methods::as(sample_points,
                                                                  "Spatial")
                                
                                raster_values <- raster::extract(x = raster,
                                                                 y = sample_points_spdf)
                                
                                sample_points$value <- raster_values
                                
                                sample_points$frame_id <- unique(frame$frame_id)
                                
                                sample_points
                              })

# Combine the two lists into a list where all draws that used the same seed are in the same sf object
# So, for the first seed there was a draw for the AOI and each frame and those are put into the same
# sf points object as the first object in this list
sample_points_list <- mapply(X = sample_points_list_1,
                             Y = sample_point_list_2,
                             Z = sample_point_list_3,
                             SIMPLIFY = FALSE,
                             FUN = function(X, Y, Z){
                               rbind(X, Y, Z)
                             })

#### Generate Thiessen polygons ####
sample_points_attributed_thiessen_list_list <- lapply(X = sample_points_list,
                                                      frame = aoi,
                                                      n_polygons = round(n_sample_points * thiessen_proportion),
                                                      points = sample_points,
                                                      points_min = thiessen_minimum_sample,
                                                      seed_number = thiessen_seed,
                                                      seed_increment = 100000,
                                                      use_albers = TRUE,
                                                      verbose = TRUE,
                                                      FUN = function(X,
                                                                     frame,
                                                                     n_polygons,
                                                                     points,
                                                                     points_min,
                                                                     seed_number,
                                                                     seed_increment,
                                                                     use_albers,
                                                                     verbose){
                                                        sample_points <- X
                                                        current_sample_seed <- unique(sample_points$sample_seed)
                                                        
                                                        thiessen_polygons <- thiessen_polygons_gen_random(frame = frame,
                                                                                                          n_polygons = n_polygons,
                                                                                                          points = sample_points,
                                                                                                          points_min = points_min,
                                                                                                          # seed_number = seed_number,
                                                                                                          seed_number = unique(sample_points$sample_seed),
                                                                                                          seed_increment = seed_increment,
                                                                                                          use_albers = use_albers,
                                                                                                          verbose = verbose)
                                                        
                                                        thiessen_polygons$tpoly_id <- paste0(thiessen_polygons$aoi_id,
                                                                                             "-",
                                                                                             thiessen_polygons$tpoly_id)
                                                        thiessen_polygons$sample_seed <- current_sample_seed
                                                        
                                                        ## Attribute the points with Thiessen weights
                                                        sample_points <- sf::st_join(x = sample_points,
                                                                                     y = thiessen_polygons[, c("tpoly_id", "weight")])
                                                        
                                                        sample_points <- sample_points[!is.na(sample_points$tpoly_id), ]
                                                        
                                                        # We're returning both the points and the polygons so that we can preserve the polygons
                                                        # Without this, they'd vanish as soon as the lapply() is done
                                                        output <- list(sample_points = sample_points,
                                                                       thiessen_polygons = thiessen_polygons)
                                                        
                                                        return(output)
                                                      })

# Make a list of just the sample points with the Thiessen info attached
sample_points_attributed_thiessen_list <- lapply(X = sample_points_attributed_thiessen_list_list,
                                                 FUN = function(X){
                                                   X[["sample_points"]]
                                                 })

# Make an sf polygon object of the Thiessen polygons
thiessen_list <- lapply(X = sample_points_attributed_thiessen_list_list,
                        FUN = function(X){
                          X[["thiessen_polygons"]]
                        })
thiessen_polygons <- do.call(rbind,
                             thiessen_list)

#### Generate wgtcat polygons ####
# Give both sets of polygons unique IDs in the same variable
aoi$uid <- "aoi"
frame_2$uid <- "frame_2"
frame_3$uid <- "frame_3"

# Combine the polygons
polygons <- rbind(aoi[, "uid"],
                  frame_2[, "uid"],
                  frame_3[, "uid"])

# Intersect them to create weight category polygons
wgtcat_polygons <- wgtcat_gen(polygons,
                              aoi_index = 1)
wgtcat_polygons$raster_id <- raster_metadata$raster_id
wgtcat_polygons$aoi_id <- aoi$aoi_id

# Attribute the sample point with wgtcat info
sample_points_attributed_wgtcat_list <- lapply(X = sample_points_list,
                                               wgtcat_polygons = wgtcat_polygons,
                                               FUN = function(X, wgtcat_polygons){
                                                 sample_points <- X
                                                 
                                                 
                                                 # Join with the points to get wgtcat_id added
                                                 sample_points <- sf::st_join(x = sample_points,
                                                                              y = wgtcat_polygons[, c("wgtcat_id")])
                                                 
                                                 # Restrict to the AOI (just to be safe)
                                                 sample_points <- sample_points[!is.na(sample_points$wgtcat_id), ]
                                                 
                                                 # Summarize the number of points in each wgtcat_id
                                                 wgtcat_summary <- table(sample_points$wgtcat_id)
                                                 wgtcat_summary <- data.frame(wgtcat_id = names(wgtcat_summary),
                                                                              n_points = as.vector(wgtcat_summary),
                                                                              stringsAsFactors = FALSE)
                                                 
                                                 # Add the point counts to wgtcat_polygons
                                                 wgtcat_polygons <- merge(x = wgtcat_polygons,
                                                                          y = wgtcat_summary,
                                                                          by = "wgtcat_id")
                                                 
                                                 # Calculate the weight for each wgtcat_id
                                                 wgtcat_polygons$weight <- wgtcat_polygons$area_m2 / wgtcat_polygons$n_points
                                                 
                                                 # Add that weight to the points
                                                 sample_points <- sf::st_join(x = sample_points,
                                                                              y = wgtcat_polygons[, c("weight")])
                                                 
                                                 sample_points
                                               })

#### Run weighted analysis ####
sample_point_summary_thiessen_list <- lapply(X = sample_points_attributed_thiessen_list,
                                             possible_values = raster_values,
                                             alpha = analysis_alpha,
                                             FUN = function(X,
                                                            possible_values,
                                                            alpha){
                                               sample_points <- X
                                               sample_point_summary <- continuous_analysis(data = sample_points,
                                                                                           alpha = alpha)
                                               sample_point_summary$sample_seed <- sample_points[["sample_seed"]][1]
                                               sample_point_summary
                                             })

sample_point_summary_thiessen <- do.call(rbind,
                                         sample_point_summary_thiessen_list)


sample_point_summary_wgtcat_list <- lapply(X = sample_points_attributed_wgtcat_list,
                                           possible_values = raster_values,
                                           alpha = analysis_alpha,
                                           FUN = function(X,
                                                          possible_values,
                                                          alpha){
                                             sample_points <- X
                                             sample_point_summary <- continuous_analysis(data = sample_points,
                                                                                         alpha = alpha)
                                             sample_point_summary$sample_seed <- sample_points[["sample_seed"]][1]
                                             sample_point_summary
                                           })

sample_point_summary_wgtcat <- do.call(rbind,
                                       sample_point_summary_wgtcat_list)

#### Summarize AOI ####
raster_values_in_aoi <- unlist(raster::extract(x = current_raster,
                                               y = aoi))

raster_summary <- switch(raster_type,
                         "categorical" = {
                           # Count the number of cells in each category
                           category_counts <- table(raster_values_in_aoi)
                           # Make that into a data frame
                           raster_summary <- data.frame(category = names(category_counts),
                                                        n = as.vector(category_counts),
                                                        stringsAsFactors = FALSE)
                           # Calculate the proportions
                           raster_summary$proportion <- raster_summary$n / sum(raster_summary$n)
                           raster_summary
                         },
                         "continuous" = {
                           raster_mean <- mean(raster_values_in_aoi)
                           raster_sd <- sd(raster_values_in_aoi)
                           raster_summary <- data.frame(mean = raster_mean,
                                                        sd = raster_sd)
                           raster_summary
                         })

#### Plot the results ####
unweighted_results <- sample_point_summary_thiessen[, c("sample_seed", "n", "mean", "sd")]
unweighted_results$weighted <- "Unweighted"
weighted_thiessen_results <- sample_point_summary_thiessen[, c("sample_seed", "n", "mean_weighted", "sd_weighted")]
weighted_thiessen_results$weighted <- "Weighted (Thiessen)"
names(weighted_thiessen_results) <- names(unweighted_results)
weighted_wgtcat_results <- sample_point_summary_wgtcat[, c("sample_seed", "n", "mean_weighted", "sd_weighted")]
weighted_wgtcat_results$weighted <- "Weighted (WgtCat)"
names(weighted_wgtcat_results) <- names(unweighted_results)

results <- rbind(unweighted_results,
                 weighted_thiessen_results,
                 weighted_wgtcat_results)

# Boxplot of means
ggplot() +
  geom_boxplot(data = results,
               aes(x = weighted,
                   y = mean)) +
  geom_hline(data = raster_summary,
             aes(yintercept = mean),
             color = "red")

# Bounds
ggplot() +
  geom_segment(data = results,
               aes(x = value,
                   xend = value,
                   y = lower_bound,
                   yend = upper_bound),
               size = 2,
               alpha = 0.01) +
  geom_point(data = raster_summary,
             aes(x = category,
                 y = proportion),
             color = "red") +
  facet_wrap(~weighted)

#### Write out results ####
# AOI
aoi_output_variables <- c("raster_id", "aoi_id", "aoi_seed")
aoi_output_path <- paste0(output_path,
                          "/",
                          "spatial",
                          "/",
                          "aoi.shp")
sf::st_write(aoi[, aoi_output_variables],
             dsn = aoi_output_path,
             driver = "ESRI Shapefile")

# Sample frames
frame_output_variables <- c("raster_id", "frame_id", "frame_seed")
frames <- rbind(frame_2,
                frame_3)[, frame_output_variables]
frame_output_path <- paste0(output_path,
                            "/",
                            "spatial",
                            "/",
                            "frames.shp")
sf::st_write(frames[, frame_output_variables],
             dsn = frame_output_path,
             driver = "ESRI Shapefile")

# Thiessen polygons
thiessen_polygons_output_variables <- c("raster_id", "aoi_id", "tpoly_id", "sample_seed", "tpoly_seed", "area_m2", "n_points", "weight")
thiessen_polygons_output_path <- paste0(output_path,
                                        "/",
                                        "spatial",
                                        "/",
                                        "thiessen_polygons.shp")
sf::st_write(thiessen_polygons[, thiessen_polygons_output_variables],
             dsn = thiessen_polygons_output_path,
             driver = "ESRI Shapefile")

# WgtCat polygons
wgtcat_output_variables <- c("raster_id", "aoi_id", "wgtcat_id", "area_m2")
wgtcat_output_path <- paste0(output_path,
                             "/",
                             "spatial",
                             "/",
                             "wgtcat.shp")
sf::st_write(wgtcat_polygons[, wgtcat_output_variables],
             dsn = wgtcat_output_path,
             driver = "ESRI Shapefile")

# Points
sample_points_attributed_thiessen <- do.call(rbind,
                                             sample_points_attributed_thiessen_list)
sample_points_attributed_thiessen$wgtcat_id <- NA
sample_points_attributed_thiessen$raster_id <- raster_metadata$raster_id
sample_points_attributed_thiessen$aoi_id <- aoi$aoi_id
sample_points_attributed_wgtcat <- do.call(rbind,
                                           sample_points_attributed_wgtcat_list)
sample_points_attributed_wgtcat$tpoly_id <- NA
sample_points_attributed_wgtcat$raster_id <- raster_metadata$raster_id
sample_points_attributed_wgtcat$aoi_id <- aoi$aoi_id

sample_points_output_variables <- c("raster_id", "aoi_id", "frame_id", "tpoly_id", "wgtcat_id", "sample_id", "sample_seed", "value", "weight")
sample_points_attributed_combined <- rbind(sample_points_attributed_thiessen[, sample_points_output_variables],
                                           sample_points_attributed_wgtcat[, sample_points_output_variables])
sample_points_output_path <- paste0(output_path,
                                    "/",
                                    "spatial",
                                    "/",
                                    "points.shp")
sf::st_write(sample_points_attributed_combined,
             dsn = sample_points_output_path,
             driver = "ESRI Shapefile")

# Results
results_output_variables <- c("raster_id", "aoi_id", "sample_seed", "weighted", "n", "mean", "sd", "mean_true", "sd_true")
results_output_path <- paste0(output_path,
                              "/",
                              "results",
                              "/",
                              "results.csv")
write.csv(x = results[, results_output_variables],
