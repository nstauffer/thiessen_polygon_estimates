library(ggplot2)
ggplot() +
  geom_sf(data = aoi) +
  geom_sf(data = frame,
          alpha = 0.5) +
  geom_sf(data = sample_point_list[[1]])

#### WEIGHTED ANALYSIS OF CONTINUOUS VARIABLE ####
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
