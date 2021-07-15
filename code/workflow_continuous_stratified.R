# FOR SIMULATING A CONTINUOUS VARIABLE WITH A STRATIFIED SAMPLING DESIGN

#### Get the packages attached ------------------------------------------------
library(ggplot2)

#### Get the functions loaded ------------------------------------------------
source("C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/functions.R")

#### Set the config ------------------------------------------------
# THIS IS ALL INHERITED FROM THE BATCH.R RUN BUT IS LEFT HERE AS REFERENCE

# # Simulation
# projection <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
# output_path <- "C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/simulations/intensification_continuous/output/"
# # simulation_seed <- 1
# 
# # Raster
# raster_type <- "continuous"
# raster_ncol <- 1000
# raster_nrow <- 1000
# raster_resolution = 1
# raster_autocorr_range = 50
# raster_mag_var = 10
# raster_nug = 0.2
# raster_mean = 1
# raster_rescale = TRUE
# raster_seed <- 420 * simulation_seed
# 
# # AOI
# aoi_n_vertices <- 6
# aoi_convex_hull <- TRUE
# aoi_seed <- 1123 * simulation_seed
# 
# # Sample
# frame_seed <- 111 * simulation_seed
# frame_n_vertices <- 6
# frame_convex_hull <- TRUE
# sample_type <- "simple"
# n_sample_points <- 50
# sample_seeds <- 1:199
# 
# # Thiessen polygons
# thiessen_distribution <- "simple"
# thiessen_n_polygons <- 5
# thiessen_seed <- 96
# thiessen_minimum_sample <- 2
# 
# # Analysis
# analysis_alpha <- 0.05

#### Generate a raster ------------------------------------------------
##### Create the raster ####
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

##### Build the metadata dataframe ####
raster_metadata <- data.frame(raster_id = paste0("raster_", raster_seed),
                              raster_seed = raster_seed,
                              raster_type = raster_type,
                              raster_ncol = raster_ncol,
                              raster_nrow = raster_nrow)

##### Build a bounding box sfc object ####
# We're also making an sf object that represents the frame for the raster
# We'll need this for generating Thiessen polygons that fully cover the AOI
raster_boundary_coords <- data.frame(x = c(0, raster_ncol, raster_ncol, 0, 0),
                                     y = c(0, 0, raster_nrow, raster_nrow, 0))
raster_boundary_sfc <- sf::st_sfc(sf::st_polygon(x = list(outer = as.matrix(raster_boundary_coords))))
raster_boundary_sf <- sf::st_sf(raster_boundary_sfc,
                                crs = projection)

#### Generate AOI ------------------------------------------------
##### Create the AOI polygon ####
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

#### Generate strata ------------------------------------------------
strata_polygons <- strata_gen(frame = aoi,
                              n_strata = n_strata,
                              type = strata_type,
                              landscape_raster = current_raster,
                              seed_number = strata_seed,
                              seed_increment = 10000)

strata_polygons$area_m2 <- as.vector(sf::st_area(strata_polygons))
strata_polygons$aoi_id <- aoi$aoi_id

#### Generate sampling points ------------------------------------------------
##### Allocate the effort to the strata proportionally by area ####
strata_polygons$proportional_area <- strata_polygons$area_m2 / sum(strata_polygons$area_m2)

strata_polygons$n_points <- round(strata_polygons$proportional_area * n_sample_points)

# Handle too many (or too few) points due to rounding by adjusting the count in the stratum with the most points
n_excess_points <- sum(strata_polygons$n_points) - n_sample_points

most_points_index <- which(strata_polygons$n_points == max(strata_polygons$n_points))

strata_polygons$n_points[most_points_index] <- strata_polygons$n_points[most_points_index] - n_excess_points

# Create a design object for grts()
base_design <- lapply(X = strata_polygons$n_points,
                      FUN = function(X) {
                        list(panel = c("1" = X),
                             seltype = "Equal",
                             over = 0)
                      })

names(base_design) <- strata_polygons$stratum_id

##### Create points within the strata ####
sample_points_list <- lapply(X = sample_seeds,
                             strata = strata_polygons,
                             design = base_design,
                             sample_type = sample_type,
                             projection = projection,
                             raster = current_raster,
                             FUN = function(X,
                                            strata,
                                            design,
                                            sample_type,
                                            projection,
                                            raster){
                               sample_seed <- X
                               sample_points <- switch(sample_type,
                                                       "simple" = {
                                                         # To draw the appropriate number of points in each stratum
                                                         points_list <- lapply(X = names(design),
                                                                               design = design,
                                                                               strata = strata,
                                                                               sample_seed = sample_seed,
                                                                               FUN = function(X, design, strata, sample_seed) {
                                                                                 
                                                                                 points_gen(frame = strata[strata$stratum_id == X, ],
                                                                                            sample_type = "simple",
                                                                                            n_points = design[[X]][["panel"]],
                                                                                            seed_number = sample_seed,
                                                                                            projection = projection)
                                                                               })
                                                         do.call(rbind,
                                                                 points_list)
                                                       },
                                                       "balanced" = {
                                                         set.seed(sample_seed)
                                                         points <- spsurvey::grts(design = design,
                                                                                  DesignID = "",
                                                                                  type.frame = "area",
                                                                                  src.frame = "sf.object",
                                                                                  sf.object = strata,
                                                                                  stratum = "stratum_id",
                                                                                  shapefile = FALSE)
                                                         
                                                         if (!("sf" %in% class(points))) {
                                                           points <- methods::as(points, "sf")
                                                         }
                                                         
                                                         # Add our ID in the format *I* want
                                                         points$sample_id <- paste0("sample_",
                                                                                    seed_number,
                                                                                    "-",
                                                                                    1:nrow(points))
                                                         
                                                         points$sample_seed <- sample_seed
                                                         
                                                         points[, c("sample_id", "sample_seed")]
                                                       })
                               
                               
                               # Attribute them with raster values
                               # First we need them as an SPDF for raster::extract()
                               sample_points_spdf <- methods::as(sample_points,
                                                                 "Spatial")
                               
                               # Then we get a vector of the values
                               raster_values <- raster::extract(x = raster,
                                                                y = sample_points_spdf)
                               
                               # And write them into the sf object!
                               sample_points$value <- raster_values
                               
                               # Add the metadata about which frame these go to
                               sample_points$frame_id <- unique(frame$frame_id)
                               
                               # Return the points
                               sample_points
                             })


#### Generate Thiessen polygons ------------------------------------------------
##### Create Thiessen polygons and attribute points ####
# This gets kinda messy because we need to keep generating Thiessen polygons until they meet the criteria
# So we're going to go through the sample points list and for each draw in that list we'll generate
# Thiessen polygons until we get a set that works, then we attribute the points,
# and return BOTH the sample points with the weights AND the Thiessen polygons
sample_points_attributed_thiessen_list_list <- lapply(X = sample_points_list,
                                                      frame = aoi,
                                                      n_polygons = thiessen_n_polygons,
                                                      points_min = thiessen_minimum_sample,
                                                      # The envelope has to be an sfc, not sf, object for sf::st_voron
                                                      envelope = raster_boundary_sfc,
                                                      seed_increment = 100000,
                                                      use_albers = TRUE,
                                                      verbose = TRUE,
                                                      FUN = function(X,
                                                                     frame,
                                                                     n_polygons,
                                                                     points_min,
                                                                     envelope,
                                                                     seed_increment,
                                                                     use_albers,
                                                                     verbose){
                                                        # For ease of reading
                                                        sample_points <- X
                                                        
                                                        # The Thiessen polygon draws will start with the same seed as the sample point draw
                                                        # This may not be the final seed number used, depending on if the Thiessen polygons all
                                                        # Contain enough of the sample points
                                                        current_sample_seed <- unique(sample_points$sample_seed)
                                                        
                                                        # OKAY! So, we're generating the Thiessen polygons here
                                                        # This function uses centroids distributed in a simple, random way
                                                        # The clipping frame is the AOI (see the lapply() arguments)
                                                        thiessen_polygons <- thiessen_polygons_gen_random(frame = frame,
                                                                                                          # How many Thiessen polygons to draw
                                                                                                          n_polygons = n_polygons,
                                                                                                          # What points the polygons are being compared against
                                                                                                          points = sample_points,
                                                                                                          # The minimum number of points that need to be in each polygon
                                                                                                          points_min = points_min,
                                                                                                          # The maximum extent of the Thiessen polygons
                                                                                                          # We're using the boundary of the raster to make sure that
                                                                                                          # the Thiessen polygons stretch to cover the whole frame/AOI
                                                                                                          envelope = envelope,
                                                                                                          seed_number = current_sample_seed,
                                                                                                          # If a Thiessen draw doesn't have enough points in each polygon
                                                                                                          # then the seed number will be incremented by this much for the
                                                                                                          # next attempt. I keep this large so that no set of tpolys will
                                                                                                          # be the same for multiple sample point draws
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
                                                        
                                                        # This shouldn't be necessary because all points should fall within the
                                                        # boundaries of the Thiessen polygons
                                                        if (any(is.na(sample_points$tpoly_id))) {
                                                          warning("Some points did not overlap with the polygons and were dropped")
                                                          sample_points <- sample_points[!is.na(sample_points$tpoly_id), ]
                                                        }
                                                        
                                                        
                                                        # We're returning both the points and the polygons so that we can preserve the polygons
                                                        # Without this, they'd vanish as soon as the lapply() is done
                                                        output <- list(sample_points = sample_points,
                                                                       thiessen_polygons = thiessen_polygons)
                                                        
                                                        return(output)
                                                      })

##### Create lists of attributed points and thiessen polygons ####
# Make a list of just the sample points with the Thiessen info attached
sample_points_attributed_thiessen_list <- lapply(X = sample_points_attributed_thiessen_list_list,
                                                 FUN = function(X){
                                                   X[["sample_points"]]
                                                 })

# Make an sf polygon object of the Thiessen polygons
# We'll write this out later
thiessen_list <- lapply(X = sample_points_attributed_thiessen_list_list,
                        FUN = function(X){
                          X[["thiessen_polygons"]]
                        })
thiessen_polygons <- do.call(rbind,
                             thiessen_list)

#### Generate wgtcat polygons ------------------------------------------------
# There's no intensification, so these are just the strata
wgtcat_polygons <- strata_polygons

wgtcat_polygons$wgtcat_id <- wgtcat_polygons$stratum_id

# Add in some metadata
wgtcat_polygons$raster_id <- raster_metadata$raster_id
wgtcat_polygons$aoi_id <- aoi$aoi_id

##### Attribute the sample point with wgtcat info ####
sample_points_attributed_wgtcat_list <- lapply(X = sample_points_list,
                                               wgtcat_polygons = wgtcat_polygons,
                                               FUN = function(X, wgtcat_polygons){
                                                 sample_points <- X
                                                 
                                                 # Get the current sample seed
                                                 current_sample_seed <- unique(sample_points$sample_seed)
                                                 
                                                 # Add it to the wgtcat polygons
                                                 wgtcat_polygons$sample_seed <- current_sample_seed
                                                 
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
                                                 
                                                 # What if no points fall in a weight category????
                                                 missing_wgtcats <- wgtcat_polygons$wgtcat_id[!(wgtcat_polygons$wgtcat_id %in% wgtcat_summary$wgtcat_id)]
                                                 
                                                 if (length(missing_wgtcats) > 0) {
                                                   missing_wgtcats_summary <- data.frame(wgtcat_id = missing_wgtcats,
                                                                                         n_points = 0,
                                                                                         stringsAsFactors = FALSE)
                                                   wgtcat_summary <- rbind(wgtcat_summary,
                                                                           missing_wgtcats_summary)
                                                 }
                                                 
                                                 
                                                 # Add the point counts to wgtcat_polygons
                                                 wgtcat_polygons <- merge(x = wgtcat_polygons,
                                                                          y = wgtcat_summary,
                                                                          by = "wgtcat_id")
                                                 
                                                 # Calculate the weight for each wgtcat_id
                                                 wgtcat_polygons$weight <- wgtcat_polygons$area_m2 / wgtcat_polygons$n_points
                                                 
                                                 # Add that weight to the points
                                                 sample_points <- sf::st_join(x = sample_points,
                                                                              y = wgtcat_polygons[, c("wgtcat_id", "weight")])
                                                 
                                                 # Get the variables right
                                                 # I'm not really sure why wgtcat_id gets an x and y version, but the values
                                                 # are identical, so we'll just strip out the ".x" from "wgtcat_id.x"
                                                 names(sample_points) <- gsub(names(sample_points),
                                                                              pattern = "\\.x$",
                                                                              replacement = "")
                                                 
                                                 # Return only the desired variables
                                                 sample_points[, c("sample_id", "sample_seed", "value", "frame_id", "wgtcat_id", "weight", "geometry")]
                                               })

sample_points_attributed_wgtcat <- do.call(rbind,
                                           sample_points_attributed_wgtcat_list)

#### Run weighted analysis ####
##### Analyze using the Thiessen weights ####
# The list is still broken up by sampling design, so we'll use a lapply() to analyze each independently
# just as the statistics gods intended
sample_point_summary_thiessen_list <- lapply(X = sample_points_attributed_thiessen_list,
                                             alpha = analysis_alpha,
                                             FUN = function(X,
                                                            alpha){
                                               # Just to keep things simple to read
                                               sample_points <- X
                                               
                                               # Run the actual analysis
                                               sample_point_summary <- continuous_analysis(data = sample_points,
                                                                                           alpha = alpha)
                                               
                                               # Add in some metadata
                                               sample_point_summary$sample_seed <- sample_points[["sample_seed"]][1]
                                               sample_point_summary$method <- "Thiessen"
                                               
                                               # Return the data frame
                                               sample_point_summary
                                             })

# Combine the results into a single object
sample_point_summary_thiessen <- do.call(rbind,
                                         sample_point_summary_thiessen_list)


##### Analyze using the weight category weights ####
sample_point_summary_wgtcat_list <- lapply(X = sample_points_attributed_wgtcat_list,
                                           alpha = analysis_alpha,
                                           FUN = function(X,
                                                          alpha){
                                             sample_points <- X
                                             sample_point_summary <- continuous_analysis(data = sample_points,
                                                                                         alpha = alpha)
                                             sample_point_summary$sample_seed <- sample_points[["sample_seed"]][1]
                                             sample_point_summary$method <- "WgtCat"
                                             sample_point_summary
                                           })

sample_point_summary_wgtcat <- do.call(rbind,
                                       sample_point_summary_wgtcat_list)

#### Summarize AOI ------------------------------------------------
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
                           raster_variance <- var(x = raster_values_in_aoi)
                           raster_summary <- data.frame(mean = raster_mean,
                                                        sd = raster_sd,
                                                        variance = raster_variance)
                           raster_summary
                         })

#### Build the results outputs ------------------------------------------------
# So, both the sample_point_summary_* objects contain the unweighted results
# But we'll grab them from the Thiessen summary just because
unweighted_results <- sample_point_summary_thiessen[, c("sample_seed", "n", "mean", "sd", "variance")]
# Write in the weighting approach (in this case "Unweighted")
unweighted_results$weighting <- "Unweighted"

# Get just the Thiessen-weighted information
weighted_thiessen_results <- sample_point_summary_thiessen[, c("sample_seed", "n", "mean_weighted", "sd_weighted", "variance_weighted")]
# Write in that that's the weighting approach here
weighted_thiessen_results$weighting <- "Weighted (Thiessen)"
# Rename the variables to match the other data frames we'll combine this with
names(weighted_thiessen_results) <- names(unweighted_results)

# Now do the same for the weight category-weighted results
weighted_wgtcat_results <- sample_point_summary_wgtcat[, c("sample_seed", "n", "mean_weighted", "sd_weighted", "variance_weighted")]
weighted_wgtcat_results$weighting <- "Weighted (WgtCat)"
names(weighted_wgtcat_results) <- names(unweighted_results)

# Combine everything into a results data frame
results <- rbind(unweighted_results,
                 weighted_thiessen_results,
                 weighted_wgtcat_results)

# Write in the true values for each measure from the raster
results$mean_true <- raster_summary$mean
results$sd_true <- raster_summary$sd
results$variance_true <- raster_summary$variance
results$raster_id <- raster_metadata$raster_id
results$aoi_id <- aoi$aoi_id


#### Write out results and spatial objects ------------------------------------------------
##### Raster ####
raster_output_path <- paste0(output_path,
                             "/",
                             "spatial",
                             "/",
                             "raster",
                             "_",
                             simulation_seed,
                             ".tif")
raster::writeRaster(current_raster,
                    filename = raster_output_path,
                    overwrite = TRUE)

##### AOI ####
aoi_output_variables <- c("raster_id", "aoi_id", "aoi_seed")
aoi_output_path <- paste0(output_path,
                          "/",
                          "spatial",
                          "/",
                          "aoi",
                          "_",
                          simulation_seed,
                          ".shp")
sf::st_write(aoi[, aoi_output_variables],
             dsn = aoi_output_path,
             driver = "ESRI Shapefile",
             append = FALSE)


##### Thiessen polygons ####
thiessen_polygons_output_variables <- c("raster_id", "aoi_id", "tpoly_id", "sample_seed", "tpoly_seed", "area_m2", "n_points", "weight")
thiessen_polygons_output_path <- paste0(output_path,
                                        "/",
                                        "spatial",
                                        "/",
                                        "thiessen_polygons",
                                        "_",
                                        simulation_seed,
                                        ".shp")
sf::st_write(thiessen_polygons[, thiessen_polygons_output_variables],
             dsn = thiessen_polygons_output_path,
             driver = "ESRI Shapefile",
             append = FALSE)

##### Weight category polygons ####
wgtcat_output_variables <- c("raster_id", "aoi_id", "wgtcat_id", "area_m2")
wgtcat_output_path <- paste0(output_path,
                             "/",
                             "spatial",
                             "/",
                             "wgtcat",
                             "_",
                             simulation_seed,
                             ".shp")
sf::st_write(wgtcat_polygons[, wgtcat_output_variables],
             dsn = wgtcat_output_path,
             driver = "ESRI Shapefile",
             append = FALSE)

##### Points ####
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
                                    "points",
                                    "_",
                                    simulation_seed,
                                    ".shp")
sf::st_write(sample_points_attributed_combined,
             dsn = sample_points_output_path,
             driver = "ESRI Shapefile",
             append = FALSE)

##### Results ####
results_output_variables <- c("raster_id", "aoi_id", "sample_seed", "weighting", "n", "mean", "sd", "variance", "mean_true", "sd_true", "variance_true")
results_output_path <- paste0(output_path,
                              "/",
                              "results",
                              "/",
                              "results",
                              "_",
                              simulation_seed,
                              ".csv")
write.csv(x = results[, results_output_variables],
          file = results_output_path)
