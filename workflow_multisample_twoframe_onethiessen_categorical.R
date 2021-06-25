#### Get the packages attached ####
library(ggplot2)

#### Get the functions loaded ####
source("C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/functions.R")

#### Set the config ####
# Simluation
n_sims <- 10
projection <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Raster
raster_type <- "categorical"
raster_ncol <- 100
raster_nrow <- 100
raster_values <- c(1, 2, 3)
raster_distribution <- "normal"
raster_mean <- NULL
raster_sd <- NULL
raster_seed <- 420

# AOI
aoi_n_vertices <- 6
aoi_convex_hull <- TRUE
aoi_seed <- 69

# Sample
frame_seed <- 42
frame_n_vertices <- 6
frame_convex_hull <- TRUE
sample_type <- "simple"
n_sample_points <- 50
sample_seeds <- 1:100

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
                           landscape_gen_continuous(max = max(raster_values),
                                                    min = min(raster_values),
                                                    distribution = raster_distribution,
                                                    ncol = raster_ncol,
                                                    nrow = raster_nrow,
                                                    mean = raster_mean,
                                                    sd = raster_sd,
                                                    seed_number = NULL,
                                                    projection = NULL)
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
                                 
                                 sample_points
                               })

# Create the sampling frame that isn't the AOI
frame <- aoi_gen(xmax = raster_ncol,
                 xmin = 0,
                 ymax = raster_nrow,
                 ymin = 0,
                 n_vertices = frame_n_vertices,
                 convex_hull = frame_convex_hull,
                 seed_number = frame_seed)

# Generate sampling points for that frame, restricted to those that fall in the AOI
sample_point_list_2 <- lapply(X = sample_seeds,
                              frame = frame,
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
                                
                                sample_points
                              })

# Combine the two lists into a list where both draws that used the same seed are in the same sf object
sample_points_list <- mapply(X = sample_points_list_1,
                             Y = sample_point_list_2,
                             SIMPLIFY = FALSE,
                             FUN = function(X, Y){
                               rbind(X, Y)
                             })

#### Generate Thiessen polygons ####
sample_points_attributed_thiessen_list <- lapply(X = sample_points_list,
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
                                          
                                          thiessen_polygons <- thiessen_polygons_gen_random(frame = frame,
                                                                                            n_polygons = n_polygons,
                                                                                            points = sample_points,
                                                                                            points_min = points_min,
                                                                                            # seed_number = seed_number,
                                                                                            seed_number = unique(sample_points$sample_seed),
                                                                                            seed_increment = seed_increment,
                                                                                            use_albers = use_albers,
                                                                                            verbose = verbose)
                                          
                                          ## Attribute the points with Thiessen weights
                                          sample_points <- sf::st_join(x = sample_points,
                                                                       y = thiessen_polygons[, c("tpoly_id", "weight")])
                                          
                                          sample_points <- sample_points[!is.na(sample_points$tpoly_id), ]
                                          
                                          sample_points
                                        })

#### Generate wgtcat polygons ####
# Give both sets of polygons unique IDs in the same variable
aoi$uid <- "aoi"
frame$uid <- "frame"

# Combine the polygons
polygons <- rbind(aoi[, "uid"],
                  frame[, "uid"])

# Intersect them to create weight category polygons
wgtcat_polygons <- sf::st_intersection(polygons)
# Add in a unique ID for them
wgtcat_polygons$wgtcat_id <- sapply(X = wgtcat_polygons$origins,
                                    FUN = function(X){
                                      paste(X, collapse = "-")
                                    })
# Drop any that are just from the frame outside the AOI
wgtcat_polygons <- wgtcat_polygons[!(wgtcat_polygons$wgtcat_id %in% c("2")), ]

# Add the areas
wgtcat_polygons$area_m2 <- as.vector(sf::st_area(wgtcat_polygons))


sample_points_attributed_wgtcat_list <- lapply(X = sample_points_list,
                                               wgtcat_polygons,
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
                                      sample_point_summary <- categorical_analysis(data = sample_points,
                                                                                   possible_values = possible_values,
                                                                                   alpha = alpha)
                                      sample_point_summary$sample_seed <- sample_points[["sample_seed"]][1]
                                      sample_point_summary$value <- paste(sample_point_summary$value)
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
                                               sample_point_summary <- categorical_analysis(data = sample_points,
                                                                                            possible_values = possible_values,
                                                                                            alpha = alpha)
                                               sample_point_summary$sample_seed <- sample_points[["sample_seed"]][1]
                                               sample_point_summary$value <- paste(sample_point_summary$value)
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
unweighted_results <- sample_point_summary_thiessen[, c("sample_seed", "value", "n", "proportion", "lower_bound", "upper_bound")]
unweighted_results$weighted <- "Unweighted"
weighted_thiessen_results <- sample_point_summary_thiessen[, c("sample_seed", "value", "n", "proportion_weighted", "lower_bound_weighted", "upper_bound_weighted")]
weighted_thiessen_results$weighted <- "Weighted (Thiessen)"
names(weighted_thiessen_results) <- names(unweighted_results)
weighted_wgtcat_results <- sample_point_summary_wgtcat[, c("sample_seed", "value", "n", "proportion_weighted", "lower_bound_weighted", "upper_bound_weighted")]
weighted_wgtcat_results$weighted <- "Weighted (WgtCat)"
names(weighted_wgtcat_results) <- names(unweighted_results)

results <- rbind(unweighted_results,
                 weighted_thiessen_results,
                 weighted_wgtcat_results)

# Boxplot of proportions
ggplot() +
  geom_boxplot(data = results,
               aes(x = value,
                   y = proportion)) +
  geom_point(data = raster_summary,
             aes(x = category,
                 y = proportion),
             color = "red") +
  facet_wrap(~weighted)

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


ggplot() +
  geom_boxplot(data = sample_point_summary,
               aes(x = value,
                   y = proportion_weighted)) +
  geom_point(data = raster_summary,
             aes(x = category,
                 y = proportion),
             color = "red")

ggplot() +
  geom_boxplot(data = sample_point_summary,
               aes(x = value,
                   y = proportion)) +
  geom_point(data = raster_summary,
             aes(x = category,
                 y = proportion),
             color = "red")

#### Write out results ####