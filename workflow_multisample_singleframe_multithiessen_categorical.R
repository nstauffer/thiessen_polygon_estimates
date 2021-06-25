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
sample_type <- "simple"
n_sample_points <- 50
sample_seeds <- 666:676

# Thiessen polygons
thiessen_distribution <- "simple"
thiessen_proportion <- 0.1
thiessen_seeds <- 96:106
thiessen_minimum_sample <- 2


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
sample_points_list <- lapply(X = sample_seeds,
                             frame = aoi,
                             sample_type = sample_type,
                             n_points = n_sample_points,
                             raster = current_raster,
                             projection = projection,
                             FUN = function(X,
                                            frame,
                                            sample_type,
                                            n_points,
                                            raster,
                                            projection){
                               # Generate the points
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
                               
                               return(sample_points)
                             })


#### Generate Thiessen polygons ####
### Maybe do this inside a while() and name the output list with the seed number?
## Generate centroids
centroids <- points_gen(frame = aoi,
                        sample_type = thiessen_distribution,
                        n_points = round(n_sample_points * thiessen_proportion),
                        seed_number = thiessen_seed,
                        projection = projection)

# Get the variables right
centroids$raster_id <- raster_metadata$raster_id
centroids$aoi_id <- aoi$aoi_id
centroids$tpoly_id <- paste0(centroids$aoi_id,
                             "-",
                             "tpoly_",
                             thiessen_seed,
                             "-",
                             1:nrow(centroids))
centroids$tpoly_seed <- thiessen_seed
centroid_variables <- c("raster_id",
                        "aoi_id",
                        "tpoly_id",
                        "tpoly_seed")
centroids <- centroids[, centroid_variables]

## Generate polygons
thiessen_polygons <- thiessen_polygons_gen(points = centroids,
                                       frame = aoi,
                                       use_albers = TRUE)
thiessen_polygons$tpoly_id <- paste0(thiessen_polygons$aoi_id,
                                     "-",
                                     "tpoly_",
                                     thiessen_seed,
                                     "-",
                                     thiessen_polygons$polygon_unique_id)

## Check that polygons contain enough points
sample_points_attributed <- sf::st_join(x = sample_points,
                                        y = thiessen_polygons[, c("raster_id", "aoi_id", "tpoly_id")])

tpoly_summary <- data.frame(tpoly_id = names(table(sample_points_attributed$tpoly_id)),
                            n_points = as.vector(table(sample_points_attributed$tpoly_id)),
                            stringsAsFactors = FALSE)

# This will draw new polygons until all of them have at least two points
while (!all(tpoly_summary$n_points > 1)) {
  thiessen_seed<- thiessen_seed + 1000000
  ## Generate centroids
  centroids <- points_gen(frame = aoi,
                          sample_type = thiessen_distribution,
                          n_points = round(n_sample_points * thiessen_proportion),
                          seed_number = thiessen_seed,
                          projection = projection)
  
  # Get the variables right
  centroids$raster_id <- raster_metadata$raster_id
  centroids$aoi_id <- aoi$aoi_id
  centroids$tpoly_id <- paste0(centroids$aoi_id,
                               "-",
                               "tpoly_",
                               thiessen_seed,
                               "-",
                               1:nrow(centroids))
  centroids$tpoly_seed <- thiessen_seed
  centroid_variables <- c("raster_id",
                          "aoi_id",
                          "tpoly_id",
                          "tpoly_seed")
  centroids <- centroids[, centroid_variables]
  
  ## Generate polygons
  thiessen_polygons <- thiessen_polygons_gen(points = centroids,
                                             frame = aoi,
                                             use_albers = TRUE)
  thiessen_polygons$tpoly_id <- paste0(thiessen_polygons$aoi_id,
                                       "-",
                                       "tpoly_",
                                       thiessen_seed,
                                       "-",
                                       thiessen_polygons$polygon_unique_id)
  
  ## Check that polygons contain enough points
  sample_points_attributed <- sf::st_join(x = sample_points,
                                          y = thiessen_polygons[, c("raster_id", "aoi_id", "tpoly_id")])
  
  tpoly_summary <- data.frame(tpoly_id = names(table(sample_points_attributed$tpoly_id)),
                              n_points = as.vector(table(sample_points_attributed$tpoly_id)),
                              stringsAsFactors = FALSE)
}

# Add the point counts to the polygons
thiessen_polygons <- merge(x = thiessen_polygons,
                           y = tpoly_summary)

# Calculate the weights
thiessen_polygons$weight <- thiessen_polygons$area_m2 / thiessen_polygons$n_points

## Attribute the points with Thiessen weights
sample_points <- sf::st_join(x = sample_points,
                    y = thiessen_polygons[, c("tpoly_id", "weight")])

#### Run weighted analysis ####
# Add the point counts to the polygons
thiessen_polygons <- merge(x = thiessen_polygons,
                           y = tpoly_summary)

# Calculate the weights
thiessen_polygons$weight <- thiessen_polygons$area_m2 / thiessen_polygons$n_points

## Attribute the points with Thiessen weights
sample_points <- sf::st_join(x = sample_points,
                             y = thiessen_polygons[, c("tpoly_id", "weight")])

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
                           return(raster_summary)
                         },
                         "continuous" = {
                           raster_mean <- mean(raster_values_in_aoi)
                           raster_sd <- sd(raster_values_in_aoi)
                           raster_summary <- data.frame(mean = raster_mean,
                                                        sd = raster_sd)
                           return(raster_summary)
                         })

#### Write out results ####