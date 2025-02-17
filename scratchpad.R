source("C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/functions.R")
library(ggplot2)

#### CLUSTER-BASED THIESSEN POLYGONS ####
raster_ncol <- 100
raster_nrow <- 100

thiessen_n_polygons <- 3

projection <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Get a polygon and points
aoi <- aoi_gen(xmin = 0,
               xmax = raster_ncol,
               ymin = 0,
               ymax = raster_nrow,
               n_vertices = 10,
               convex_hull = TRUE,
               seed_number = 420)

sample_points <- points_gen(frame = aoi,
                            sample_type = "simple",
                            n_points = 16,
                            seed_number = 70)

test <- thiessen_polygons_gen_clustered(frame = aoi,
                                        points = sample_points,
                                        n_polygons = 5)

ggplot() + 
  geom_sf(data = aoi) +
  geom_sf(data = test) +
  geom_sf(data = sample_points)

# Get the point coordinates. We'll need them to calcualte distances
points_coords <- as.data.frame(sf::st_coordinates(points))
names(points_coords) <- c("x", "y")

# Get a distance matrix
points_distance_matrix <- geosphere::distm(x = points_coords)

# Do some hierarchical clustering based on the distances
hierarchical_clusters <- hclust(as.dist(m = points_distance_matrix),
                                method = "complete")

# Put them into a number of clusters matching the Thiessen polygon count
cluster_membership <- cutree(tree = hierarchical_clusters,
                             k = thiessen_n_polygons)

# Write that info into the points object
points$cluster <- cluster_membership
points_coords$cluster <- cluster_membership

ggplot() + 
  geom_sf(data = aoi) +
  geom_sf(data = points,
          aes(color = cluster_membership))

# For each cluster, make an sf object for the centroid
centroid_sf_list <- lapply(X = split(points_coords, points_coords$cluster),
                         projection = projection,
                         FUN = function(X,
                                        projection) {
                           coords <- X
                           current_cluster <- coords$cluster[1]
                           centroid_x <- mean(coords$x)
                           centroid_y <- mean(coords$y)
                           
                           centroid_df <- data.frame(cluster = current_cluster,
                                                     x = centroid_x,
                                                     y = centroid_y)
                           
                           coords_matrix <- as.matrix(centroid_df[, c("x", "y")])
                           
                           centroid_sfc <- sf::st_point(x = coords_matrix)
                           
                           # For some reason, I have to do this to feed into sf::st_sf()
                           # instead of just giving it the sfc object
                           centroid_sfc_geometry <- sf::st_geometry(centroid_sfc)
                           
                           centroid_sf <- sf::st_sf(centroid_sfc_geometry,
                                                    crs = projection)
                           
                           centroid_sf$cluster <- current_cluster
                           
                           centroid_sf
                         })

centroids_sf <- do.call(rbind,
                        centroid_sf_list)

ggplot() + 
  geom_sf(data = aoi) +
  geom_sf(data = points,
          aes(color = cluster_membership)) +
  geom_sf(data = centroids_sf,
          color = "red")

tpolys <- thiessen_polygons_gen_fixed(centroids = centroids_sf,
                                      frame = aoi)

ggplot() + 
  geom_sf(data = aoi) +
  geom_sf(data = tpolys,
          alpha = 0.25) +
  geom_sf(data = points,
          aes(color = cluster_membership)) +
  geom_sf(data = centroids_sf,
          color = "red")

#### WILCOXON INTERPRETATION ####
summarize_wilcoxon <- function(results,
                               thresholds) {
  if (class(results) != "data.frame") {
    stop("`results` must be a data frame.")
  }
  if (!all(c("weighting", "p_value") %in% names(results))) {
    stop("`results` must contain the variables 'weighting' and 'p_value'.")
  }
  if (!is.numeric(thresholds)) {
    stop("`thresholds` must be numeric.")
  }
  if (any(thresholds > 1 | thresholds <= 0)) {
    stop("All values in `thresholds` must be between 0 and 1.")
  }
  
  # Create single-column data frames for each threshold indicating
  # if the p value was above the threshold
  checked_list <- lapply(X = thresholds,
                         results = wilcoxon_results,
                         FUN = function(X, results) {
                           p_value_threshold <- X
                           var_name <- paste0("p_value_above_",
                                              p_value_threshold)
                           check_vector <- results$p_value > p_value_threshold
                           output <- data.frame("above_p_value" = check_vector)
                           names(output) <- var_name
                           output
                         })
  
  checked_df <- do.call(cbind,
                        checked_list)
  
  wilcoxon_results_checked <- cbind(wilcoxon_results,
                                    checked_df)
  
  # Summarize by weighting approach
  wilcoxon_checked_summary_list <- lapply(X = split(wilcoxon_results_checked,
                                                    wilcoxon_results_checked$weighting),
                                          p_value_thresholds = thresholds,
                                          FUN = function(X, p_value_thresholds) {
                                            current_results <- X
                                            
                                            n_observations <- nrow(current_results)
                                            
                                            output <- data.frame("weighting" = current_results$weighting[1])
                                            
                                            for (p_value_threshold in p_value_thresholds) {
                                              current_var <- paste0("p_value_above_",
                                                                    p_value_threshold)
                                              output[[paste0("n_", current_var)]] <- sum(current_results[[current_var]])
                                              output[[paste0("proportion_", current_var)]] <- sum(current_results[[current_var]]) / n_observations
                                            }
                                            
                                            output
                                          })
  
  wilcoxon_checked_summary <- do.call(rbind,
                                      wilcoxon_checked_summary_list)
}

test <- summarize_wilcoxon(results = wilcoxon_results,
                           thresholds = c(0.2, 0.1, 0.05, 0.01))

p_value_thresholds <- c(0.2, 0.1, 0.05, 0.01)

# Create single-column data frames for each threshold indicating
# if the p value was above the threshold
checked_list <- lapply(X = p_value_thresholds,
                       results = wilcoxon_results,
                       FUN = function(X, results) {
                         p_value_threshold <- X
                         var_name <- paste0("p_value_above_",
                                            p_value_threshold)
                         check_vector <- results$p_value > p_value_threshold
                         output <- data.frame("above_p_value" = check_vector)
                         names(output) <- var_name
                         output
                       })

checked_df <- do.call(cbind,
                      checked_list)

wilcoxon_results_checked <- cbind(wilcoxon_results,
                                  checked_df)

# Summarize by weighting approach
wilcoxon_checked_summary_list <- lapply(X = split(wilcoxon_results_checked,
                                                  wilcoxon_results_checked$weighting),
                                        p_value_thresholds = p_value_thresholds,
                                        FUN = function(X, p_value_thresholds) {
                                          current_results <- X
                                          
                                          n_observations <- nrow(current_results)
                                          
                                          output <- data.frame("weighting" = current_results$weighting[1])
                                          
                                          for (p_value_threshold in p_value_thresholds) {
                                            current_var <- paste0("p_value_above_",
                                                                  p_value_threshold)
                                            output[[paste0("n_", current_var)]] <- sum(current_results[[current_var]])
                                            output[[paste0("proportion_", current_var)]] <- sum(current_results[[current_var]]) / n_observations
                                          }
                                          
                                          output
                                        })

wilcoxon_checked_summary <- do.call(rbind,
                                    wilcoxon_checked_summary_list)

#### TESTING WEIGHTING ####
test_list <- lapply(X = sample_points_attributed_thiessen_list,
                    FUN = function(X) {
                      X$weight <- 1
                      continuous_analysis(data = X,
                                          alpha = 0.05)
                    })
test <- do.call(rbind,
                test_list)

identical(round(test$mean, digits = 6), round(test$mean_weighted, digits = 6))
identical(round(test$sd, digits = 6), round(test$sd_weighted, digits = 6))
identical(round(test$variance, digits = 6), round(test$variance_weighted, digits = 6))

#### RASTER TO POLYGOJN FOR STRATA ####
library(ggplot2)
projection <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
raster_ncol <- 100
raster_nrow <- 100
raster_resolution = 1
raster_autocorr_range = 10
raster_mag_var = 10
raster_nug = 0.2
raster_mean = 1
raster_rescale = TRUE
raster_seed <- 666

landscape_raster <- NLMR::nlm_gaussianfield(ncol = raster_ncol,
                                            nrow = raster_nrow,
                                            resolution = raster_resolution,
                                            autocorr_range = raster_autocorr_range,
                                            mag_var = raster_mag_var,
                                            nug = raster_nug,
                                            mean = raster_mean,
                                            user_seed = raster_seed,
                                            rescale = raster_rescale)
raster::projection(landscape_raster) <- projection
raster_plotting_df <- data.frame(raster::rasterToPoints(landscape_raster))
frame <- aoi_gen(xmax = raster_ncol,
                 xmin = 0,
                 ymax = raster_nrow,
                 ymin = 0,
                 n_vertices = 6,
                 convex_hull = TRUE,
                 seed_number = 112358)

n_strata <- 3
type <- "partitioned"
landscape_raster <- landscape_raster
seed_number <- NULL
seed_increment <- 10000

strata_polygons <- strata_gen(frame = frame,
                              n_strata = n_strata,
                              type = type,
                              landscape_raster = landscape_raster,
                              seed_number = seed_number,
                              seed_increment = seed_increment)

ggplot() +
  # geom_raster(data = raster_plotting_df,
  #             aes(x = x,
  #                 y = y,
  #                 fill = layer)) +
  # geom_sf(data = frame,
  #         alpha = 0.25) +
  geom_sf(data = strata_polygons,
          aes(fill = stratum_id),
          alpha = 0.1) +
  geom_sf(data = sample_points_list[[1]])

#### RASTER BOUNDARY SFC ####
projection <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
raster_boundary_coords <- data.frame(x = c(0, raster_ncol, raster_ncol, 0, 0),
                                     y = c(0, 0, raster_nrow, raster_nrow, 0))
raster_boundary_sfc <- sf::st_sfc(sf::st_polygon(x = list(outer = as.matrix(raster_boundary_coords))))
raster_boundary_sf <- sf::st_sf(raster_boundary_sfc,
                                crs = projection)
frame <- raster_boundary_sf

ggplot() +
  geom_sf(data = raster_boundary_sfc,
          alpha = 0.25)

#### THRESHOLD ANALYSIS ####
#' Check whether results are within a certain tolerance of the true value
#' @param data Data frame. Must contain the variables \code{variable} and \code{comparison_variable}.
#' @param varaible Character string. The name of the variable in \code{data} containing the values to check against the tolerance.
#' @param comparison_variable Character string. The name of the variable in \code{data} containing the values to to use to calculate the tolerance.
#' @param percent_tolerance Numeric. The percent difference from the value in \code{comparison_variable} that its paired value in \code{variable} is allowed to be.
#' @returns Logical vector. The value is \code{TRUE} for any index where the value in \code{variable} was within the permitted tolerance and \code{FALSE} for all other indices.
tolerance_test <- function(data,
                           variable,
                           comparison_variable,
                           percent_tolerance = 5){
  if (class(data) != "data.frame") {
    stop("data must be a data frame")
  }
  if (nrow(data) < 1) {
    stop("data must contain at least one row of values")
  }
  if (!(variable %in% names(data))) {
    stop(paste0("The variable ", variable, " is missing from data"))
  }
  if (!(comparison_variable %in% names(data))) {
    stop(paste0("The variable ", comparison_variable, " is missing from data"))
  }
  if (percent_tolerance < 0 | percent_tolerance > 100) {
    stop("percent_tolerance must be a value between 0 and 100")
  }
  
  proportion_tolerance <- percent_tolerance / 100
  magnitude_tolerance <- abs(proportion_tolerance * data[[comparison_variable]])
  magnitude_difference <- abs(data[[variable]] - data[[comparison_variable]])
  magnitude_difference < magnitude_tolerance
}

#### PLOTS ####
library(ggplot2)
raster_plotting_df <- data.frame(raster::rasterToPoints(current_raster))
raster_plotting_df$category <- paste(raster_plotting_df$layer)

# Categorical
ggplot() +
  geom_raster(data = raster_plotting_df,
              aes(x = x,
                  y = y,
                  fill = category)) +
  scale_fill_viridis_d() +
  geom_sf(data = aoi,
          alpha = 0.25) +
  geom_sf(data = frame,
          alpha = 0.25) +
  geom_sf(data = sample_points_list[[1]])

# Continuous
ggplot() +
  # geom_raster(data = raster_plotting_df,
  #             aes(x = x,
  #                 y = y,
  #                 fill = layer)) +
  scale_fill_viridis_c() +
  # geom_sf(data = raster_boundary_sf,
  #         alpha = 0.25) +
  geom_sf(data = aoi,
          alpha = 0.25) +
  # geom_sf(data = wgtcat_polygons,
  #         alpha = 0.25) +
  # geom_sf(data = frame,
  #         alpha = 0.25) +
  geom_sf(data = thiessen_list[[3]],
          alpha = 0.25) +
  # geom_sf(data = thiessen_polygons,
  #         alpha = 0.25) +
  # geom_sf(data = thiessen_polygons_clipped,
  #         alpha = 0.25) +
  # geom_sf(data = test,
  #         alpha = 0.25) +
  # geom_sf(data = centroids,
  #         color = "red") +
  geom_sf(data = sample_points_list[[3]],
          aes(color = frame_id))


#### WEIGHTED ANALYSIS OF CONTINUOUS VARIABLE ####
test <- weighted_variance(values = sample_points_attributed_thiessen$value,
                          weights = sample_points_attributed_thiessen$weight,
                          na_remove = FALSE)


data <- sample_points_attributed_list[[1]]
alpha <- 0.2

continuous_analysis <- function(data,
                                alpha) {
  data$weighted_value <- data$value * data$weight / sum(data$weight)
  n <- nrow(data)
  mean <- sum(data$weighted_value)
  sd <- sqrt(sum(data$weight * (data$value - mean)^2) / ((n - 1) / n * sum(data$weight)))
  
  output <- data.frame(n = n,
                       mean = mean,
                       sd = sd)
  
  return(output)
}

test <- continuous_analysis(data = sample_points_attributed_list[[1]],
                            alpha = analysis_alpha)

### WILCOXON SIGNED RANK TEST ####
# Categorical
results_test <- merge(x = sample_point_summary_thiessen,
                      y = raster_summary,
                      by.x = "value",
                      by.y = "category",
                      all.x = TRUE)

test <- wilcox.test(x = results_test$proportion_weighted,
                    y = results_test$proportion.y,
                    paired = TRUE)

# Continuous
results_test <- sample_point_summary_thiessen
results_test$mean_raster <- raster_summary$mean

test <- wilcox.test(x = results_test$mean_weighted,
                    y = results_test$mean_raster,
                    paired = TRUE)

#### LANDSCAPE RASTERS ####
projection <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
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
raster_n_categories <- 3

test_raster <- NLMR::nlm_gaussianfield(ncol = raster_ncol,
                                       nrow = raster_nrow,
                                       resolution = raster_resolution,
                                       autocorr_range = raster_autocorr_range,
                                       mag_var = raster_mag_var,
                                       nug = raster_nug,
                                       mean = raster_mean,
                                       user_seed = raster_seed,
                                       rescale = raster_rescale)

raster::projection(test_raster) <- projection

category_increment <- 1 / raster_n_categories

for (category in raster_n_categories:1) {
  test_raster[test_raster >= ((category - 1) * category_increment) & test_raster < (category * category_increment)] <- category
}

raster::plot(test_raster)

frame <- aoi_gen(xmin = 0,
                 xmax = 100,
                 ymin = 0,
                 ymax = 100,
                 n_vertices = 5,
                 convex_hull = TRUE,
                 seed_number = 666,
                 projection = projection)

raster_values_in_aoi <- unlist(raster::extract(x = test_raster,
                                               y = frame))

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


#### THIESSEN POLYGONS ####
aea_proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Load in sample frame
# sample_frame <- sf::st_read("C:/Users/Nelson/Documents/Projects/Sample_Designs_2021/california/walkerridge/data/ca_walkerridge_2021_data.gdb",
# layer = "sample_frame")
# sample_frame <- sf::st_read("C:/Users/Nelson/Documents/Projects/Sample_Designs_2021/colorado/whiteriver/data/co_whiteriver_2021_data.gdb",
#                             layer = "sample_frame")
# sample_frame <- sf::st_read("C:/Users/Nelson/Documents/Projects/Sample_Designs_2021/arizona/kingman/data/az_kingman_2021_data.gdb",
#                             layer = "sample_frame")
sample_frame <- sf::st_read("C:/Users/Nelson/Documents/Projects/Sample_Designs_2021/utah/moab/data/ut_moab_2021_data.gdb",
                            layer = "sample_frame")
# Convert the CRS of the sample_frame to match
sample_frame <- sf::st_transform(x = sf::st_zm(sample_frame,
                                               drop = TRUE),
                                 crs = aea_proj)

# Convert to Spatial Polygons Data Frame
# This is because I need an sfc object, but I can't convert directly from sf to sfc
# without getting a deeply uninformative error message with no clear resolution
sample_frame_spdf <- methods::as(sample_frame, "Spatial")

# Also repair this!
sample_frame_spdf <- aim.analysis::repair_geometry(sample_frame_spdf,
                                                   force = TRUE,
                                                   verbose = TRUE)

# Convert to sfc object because that's what's needed for the envelope argument of st_voronoi()
sample_frame_sfc <- sf::st_as_sfc(sample_frame_spdf)

# And an sf object again because we've repaired it. I know this is inefficient.
sample_frame <- sf::st_as_sf(sample_frame_spdf)

# Generate design
design <- sample.design::allocate_panels(spdf = sample_frame_spdf,
                                         stratum_field = "sf_id",
                                         panel_names = c("test"),
                                         panel_sample_size = 25,
                                         points_min = 3,
                                         oversample_proportion = 0,
                                         oversample_min = 0)

# Generate points
points_spdf <- sample.design::grts_aim(design_object = design,
                                       design_name = "Test",
                                       source_frame = "sp.object",
                                       frame = sample_frame_spdf,
                                       stratum_field = "sf_id",
                                       seed_number = 420)
# Convert the CRS to match
points_spdf <- sp::spTransform(x = points_spdf,
                               CRSobj = aea_proj)

# Convert to an sf object
points <- methods::as(points_spdf, "sf")

# Plot the situation
base_map <- ggplot() +
  # Using the sfc object instead of the sf just to prove to myself that it's valid
  geom_sf(data = sample_frame_sfc,
          fill = "gray95") +
  geom_sf(data = points,
          color = "blue")

base_map


# Draw Thiessen polygons
# Here's where it gets weird
# The points need to be a multipoint object, apparently
points_multipoint <- sf::st_combine(sf::st_geometry(points))

# Generate the Thiessen polygons, trying to use the sample frame as an envelope
# Note that it does jack shit with the envelope argument if the envelope polygons are smaller than the default boundaries
thiessen_polygons_raw <- sf::st_voronoi(points_multipoint,
                                        envelope = sample_frame_sfc)

# This is making sure that we only have polygon features
thiessen_polygons_raw <- sf::st_collection_extract(thiessen_polygons_raw,
                                                   type = "POLYGON")

# Convert the polygons to an sf object
thiessen_polygons <- sf::st_sf(thiessen_polygons_raw)

# Plot for reassurance
ggplot() +
  geom_sf(data = thiessen_polygons_raw,
          fill = "gray95")

# Clip to sample frame
thiessen_polygons_clipped <- sf::st_intersection(x = thiessen_polygons,
                                                 y = sample_frame)

# Plot once again
ggplot() +
  geom_sf(data = thiessen_polygons_clipped,
          fill = "gray95") +
  geom_sf(data =points,
          color = "blue")

# Add in the areas for the polygons
thiessen_polygons_clipped$area_m2 <- sf::st_area(x = thiessen_polygons_clipped)

thiessen_polygons_clipped$polygon_unique_id <- 1:nrow(thiessen_polygons_clipped)

# Plot once again
ggplot() +
  geom_sf(data = thiessen_polygons_clipped,
          aes(fill = polygon_unique_id)) +
  geom_sf(data =points,
          color = "orange")

# Attribute the points with the areas!
points_attributed <- sf::st_intersection(x = points,
                                         y = thiessen_polygons_clipped)





test <- thiessen_polygons(points = points,
                          frame = sample_frame,
                          use_albers = TRUE)
# Plot once again
ggplot() +
  geom_sf(data = test,
          aes(fill = polygon_unique_id)) +
  geom_sf(data =points,
          color = "orange")


########## CONCAVE POLYGONS ##########
projection <- sp::CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
seed_number <- 69

polygon <- aoi_gen(xmin = 0,
                   xmax = 100,
                   ymin = 0,
                   ymax = 100,
                   n_vertices = 5,
                   convex_hull = TRUE,
                   seed_number = 420,
                   projection = projection)



ggplot() +
  geom_sf(data = polygon)

polygon_bbox <- sf::st_bbox(polygon)

set.seed(seed_number)
x_coords <- runif(n = 10,
                  min = polygon_bbox["xmin"],
                  max = polygon_bbox["xmax"])
set.seed(seed_number + 1)
y_coords <- runif(n = 10,
                  min = polygon_bbox["ymin"],
                  max = polygon_bbox["ymax"])
coords <- data.frame(x_coord = x_coords,
                     y_coord = y_coords)

# Making the coordinates into points this way so that every observation in the data frame is a point
potential_vertices <- sf::st_as_sf(x = coords,
                                   coords = c("x_coord", "y_coord"),
                                   crs = projection)

# Which of those are inside the polygon?
concave_vertex_indices <- as.vector(sf::st_intersects(x = potential_vertices,
                                                      y = polygon,
                                                      sparse = FALSE))

# Keep only the inside-the-polygon points
valid_concave_vertices <- potential_vertices[concave_vertex_indices, ]

ggplot() +
  geom_sf(data = polygon) +
  geom_sf(data = valid_concave_vertices)

# Get the polygon coordinates
polygon_coordinates <- sf::st_coordinates(polygon)

# Keep only the unique ones
# polygon_coordinates <- polygon_coordinates[1:(nrow(polygon_coordinates) - 1), c("X", "Y")]

valid_concave_vertex_coordinates <- sf::st_coordinates(valid_concave_vertices)

# Get a random insertion point for the new vertex
set.seed(seed_number)
insertion_index <- sample(x = 1:(nrow(polygon_coordinates) - 1),
                          size = 1)

# Make new vectors for x and y coordinates inserting the first valid concave vertex into the vector
x_coords_new <- c(polygon_coordinates[1:insertion_index, "X"],
                  valid_concave_vertex_coordinates[1, "X"],
                  polygon_coordinates[(insertion_index + 1):nrow(polygon_coordinates), "X"])
y_coords_new <- c(polygon_coordinates[1:insertion_index, "Y"],
                  valid_concave_vertex_coordinates[1, "Y"],
                  polygon_coordinates[(insertion_index + 1):nrow(polygon_coordinates), "Y"])

# Combine into a data frame to convert into an sf object
coords_new <- data.frame(x_coord = x_coords_new,
                         y_coord = y_coords_new)

# Convert into a point sf object
new_vertices <- sf::st_as_sf(x = coords_new,
                             coords = c("x_coord", "y_coord"))

# Combine all the vertex points and cast them as a polygon
concave_polygon <- sf::st_sf(sf::st_cast(sf::st_combine(new_vertices), "POLYGON"))

# Write in the original properties
for (variable in names(polygon)) {
  concave_polygon[[variable]] <- polygon[[variable]]
}

# Make sure the output is only the original properties
output <- concave_polygon[, names(polygon)]

ggplot() +
  geom_sf(data = output)



#### SAMPLING POINTS ####
source("C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/functions.R")
projection <- sp::CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
seed_number <- 69

sample_type <- "simple"
# sample_type <- "balanced"
# sample_type <- "cluster"
frame <- aoi_gen(xmin = 0,
                 xmax = 100,
                 ymin = 0,
                 ymax = 100,
                 n_vertices = 5,
                 convex_hull = TRUE,
                 seed_number = 420,
                 projection = projection)
n_points <- 50
seed_number <- 666

# do the correct kind of sample draw
points <- switch(sample_type,
                 "simple" = {
                   set.seed(seed_number)
                   sp::spsample(x = methods::as(frame, "Spatial"),
                                n = n_points,
                                type = "random")
                 },
                 "balanced" = {
                   design <- list(None = list(panel = c("1" = n_points),
                                              seltype = "Equal",
                                              over = 0))
                   set.seed(seed_number)
                   points <- spsurvey::grts(design = design,
                                            DesignID = "",
                                            type.frame = "area",
                                            src.frame = "sf.object",
                                            sf.object = frame,
                                            shapefile = FALSE)
                 },
                 "cluster" = {
                   set.seed(seed_number)
                   
                 })

# Convert from SPDF to sf
points <- methods::as(points, "sf")

points$sample_id <- paste0("sample_",
                           seed_number,
                           "-",
                           1:nrow(points))

points$sample_seed <- seed_number

output <- points[, c("sample_id", "sample_seed")]

ggplot() +
  geom_sf(data = frame) +
  geom_sf(data = output)
