#### SETUP #####################################################################
# Making estimates of percentages and acres of ecoregional polygons

##### Packages and functions ---------------------------------------------------
# trex has functions for accessing the Landscape Data Commons API and can be
# used to retrieve the data for the areas of interest
# remotes::install_github(repo = "landscape-data-commons/trex")
library(trex)

# terradactyl has functions for calculating ecological indicators (e.g., cover)
# from the raw data.
# remotes::install_github(repo = "landscape-data-commons/terradactyl")
library(terradactyl)

# This is necessary because simply calling from the namespace doesn't access the
# support functions in the package (or maybe from support packages it wants to
# attach?) in order to find the density distribution of sampling.
library(spatstat)

# Thiessen polygon function filepath
# This is where the various functions for weighted areal estimates using
# Thiessen polygons and density partitions (TPs and DPs) live. Eventually these
# will be turned into a package.

tp_dp_function_filepath <- file.path("C:/Users/Nelson.Stauffer/OneDrive - USDA/Documents/Projects/thiessen_polygon_estimates/code",
                                     "functions.R")
source(tp_dp_function_filepath)


##### Filepaths ----------------------------------------------------------------
# Folder containing the locally-stored source data
data_path <- "C:/Users/Nelson.Stauffer/OneDrive - USDA/Documents/Projects/thiessen_polygon_estimates/analysis/ecoregions_pj/data"

# Folder to write outputs to
output_path <- "C:/Users/Nelson.Stauffer/OneDrive - USDA/Documents/Projects/thiessen_polygon_estimates/analysis/ecoregions_pj/output"

# Input polygon shapefile filename
polygons_filename <- "Regions_PJ.shp"

# Output geodatabase filename. Spatial outputs will be written here.
output_gdb_filename <- "ecoregional_pinyon_analysis_20250822.gdb"


##### Parameters ---------------------------------------------------------------
analysis_name <- "ecoregional_pinyon_cover"
file_date <- "20250822"

# The level of confidence for calculating the confidence intervals around the
# estimates during analysis.
# AIM defaults to 80%
conf <- 80
continuous_indicators <- "pinyon_cover"
categorical_indicators <- "cover_category"

###### Landscape Data Commons --------------------------------------------------
# If you want to retrieve the data from the LDC, use these.
# Otherwise, set use_ldc to FALSE and provide the filename of the geodatabase in
# the data_path that contains the AIM data. The code below assumes that the
# provided geodatabase follows the standard format and contains tables with the
# default names.
# Using a geodatabase might be slower because there will be a reformatting step
# that converts the raw data from the original wide format into a long format
# whereas the LDC already has the LPI data stored in the long format.
use_ldc <- FALSE
ldc_username <- ""
ldc_password <- ""

# These won't need to be separate in the near future, but right now AIM and LMF
# are stored in different geodatabases.
aim_geodatabase_filename <- "BLM_Natl_AIM_TerrADat_Public.gdb"
lmf_geodatabase_filename <- "BLM_Natl_AIM_LMF_Public.gdb"

###### Polygon variables -------------------------------------------------------
# Which variables in the polygons contain unique identifiers.
# Note that this will not be unique combinations across multiple variables if
# more than one variable is provided. Instead it will use unique values in each
# provided variable without considering other variables.
polygon_uid_variables <- c("US_L3NAME",
                           "US_L4NAME")

###### Pinyon-Juniper species --------------------------------------------------
# These are the species codes in the data which represent pinyon or juniper
# species. These were provided by Chris Witts.
# The juniper species are unused in this process.
pinyon_species <- c("PIED",
                    "PIMO",
                    "PIDI3")
# juniper_species <- c("JUCO6",
#                     "JUCOC3",
#                     "JUCOC4",
#                     "JUCOD",
#                     "JUCOM",
#                     "JUCOM2",
#                     "JUOS",
#                     "JUDE2",
#                     "JUMO",
#                     "JUMOG",
#                     "JUMOM",
#                     "JUSC2")


###### Cover categories --------------------------------------------------------
# What're the ranges to classify the cover into?
# These will be considered inclusive only at the minimum.
pinyon_cover_classes <- list("No pinyon" = c(min = 0,
                                         max = 0.00001),
                             "1-9%" = c(min = 0.00001,
                                        max = 10),
                             "Over 9%" = c(min = 10,
                                           max = Inf))


#### READING ###################################################################
# Because things like buffering won't work with degrees as units, we'll use this
# to transform the geographic inputs sometimes.
aea_proj <- "+proj=aea +lat_1=29.5 +lat_2=42.5"

##### Polygons -----------------------------------------------------------------
# Read in the polygons and make a list of the unique subsets according to the
# values in the polygon_uid_variables vector.
polygons <- sf::st_read(dsn = file.path(data_path,
                                        polygons_filename)) |>
  sf::st_transform(x = _,
                   crs = aea_proj)

polygons_list <- lapply(X = polygon_uid_variables,
                        polygons = polygons,
                        FUN = function(X, polygons){
                          # Get just the current UID variable with the geometry.
                          current_polygons <- dplyr::select(.data = polygons,
                                                            tidyselect::all_of(x = c("uid" = X))) |>
                            # We can't work with a third dimension, so this is
                            # to toss out those data if they're here.
                            sf::st_zm(x = _,
                                      drop = TRUE)
                          
                          # For each unique value in the current UID variable,
                          # we'll grab the associated polygons and dissolve them
                          # which'll give us a list of polygons with only one
                          # record in the data frame.
                          lapply(X = unique(current_polygons$uid),
                                 current_polygons = current_polygons,
                                 current_polygon_type = X,
                                 FUN = function(X, current_polygons, current_polygon_type){
                                   # Dissolve!
                                   current_polygon_subset <- dplyr::filter(.data = current_polygons,
                                                                           uid %in% X) |>
                                     sf::st_union(x = _) |>
                                     sf::st_as_sf(x = _) |>
                                     dplyr::mutate(.data = _,
                                                   uid = X,
                                                   polygon_type = current_polygon_type)
                                   current_polygon_subset
                                 })
                          
                        }) |>
  # This unlists one level so that we don't have nested lists specific to each
  # UID variable.
  unlist(x = _,
         recursive = FALSE)

# Grab LMF rangeland polygons??????????
# If we can find appropriate rangeland polygons, it'd be nice to make estimates
# to only rangeland and not other types.

# This is the overall frame, which we can use to query the Landscape Data
# Commons to get only the data that're relevant.
polygons_extent <- sf::st_union(x = polygons) |>
  sf::st_as_sf(x = _) |>
  dplyr::mutate(.data = _,
                frame_id = dplyr::row_number())

##### Make the output polygons -------------------------------------------------
output_polygons <- dplyr::bind_rows(polygons_list) |>
  dplyr::select(.data = _,
                tidyselect::all_of(x = c("uid",
                                         "polygon_type"))) |>
  dplyr::mutate(.data = _,
                polygon_type = stringr::str_extract(string = polygon_type,
                                                    pattern = "L\\d")) |>
  # Eventually we'll write these out to a geodatabase, but they have to all be
  # the same geometry type, so we're going to make even the single-part polygons
  # into multipolygons.
  sf::st_cast(x = _,
              to = "MULTIPOLYGON")
output_polygons$area_acres <- sf::st_area(output_polygons) |>
  as.numeric(x = _) / 4047

##### Raw data -----------------------------------------------------------------
# Snag the relevant AIM/LMF LPI data from the Landscape Data Commons.
# Alternatively, we could use a downloaded geodatabase.
if (use_ldc) {
  ###### Data from the LDC -----------------------------------------------------
  # Some sanitization to make sure that we properly pass in NULL if the script
  # has the default empty string above.
  if (nchar(ldc_username) < 1) {
    ldc_username <- NULL
  }
  if (nchar(ldc_password) < 1) {
    ldc_username <- NULL
  }
  headers <- trex::fetch_ldc_spatial(polygons = polygons_extent,
                                     data_type = "header",
                                     username = ldc_username,
                                     password = ldc_password)
  
  lpi_data <- trex::fetch_ldc_spatial(polygons = polygons_extent,
                                      data_type = "lpi",
                                      username = ldc_username,
                                      password = ldc_password)
} else {
  ###### Data from a geodatabase -----------------------------------------------
  available_layers <- dplyr::bind_rows(sf::st_layers(dsn = file.path(data_path,
                                                                     aim_geodatabase_filename)) |>
                                         dplyr::pull(.data = _,
                                                     name) |>
                                         as.data.frame() |>
                                         setNames(object = _,
                                                  nm = "name"),
                                       sf::st_layers(dsn = file.path(data_path,
                                                                     lmf_geodatabase_filename)) |>
                                         dplyr::pull(.data = _,
                                                     name) |>
                                         as.data.frame() |>
                                         setNames(object = _,
                                                  nm = "name"))
  # We'll make a lookup table with the actual layer names and the versions
  # without the prefixes so we can easily reference those.
  layers_lookup <- data.frame(internal_layer_name = available_layers$name) |>
    dplyr::mutate(.data = _,
                  # layer_name = stringr::str_extract(string = internal_layer_name,
                  #                                   pattern = "(?<=_[IFD]_)[A-z]+$") |>
                  #   unlist(),
                  layer_name = internal_layer_name,
                  source = dplyr::case_when(stringr::str_detect(string = internal_layer_name,
                                                                pattern = "^tbl") ~ "AIM",
                                            stringr::str_detect(string = internal_layer_name,
                                                                pattern = "^[A-Z]") ~ "LMF",
                                            .default = NA))
  # We're going to need to do this for AIM and LMF separately because they have
  # different formats in the geodatabase.
  input_data_list <- list(aim = list(),
                          lmf = list(),
                          lpi_tall = list())
  ####### AIM ------------------------------------------------------------------
  if ("tblPlots" %in% layers_lookup$layer_name) {
    input_data_list[["aim"]][["headers"]] <- sf::st_read(dsn = file.path(data_path,
                                                                         aim_geodatabase_filename),
                                                         layer = layers_lookup$internal_layer_name[layers_lookup$layer_name == "tblPlots"]) |>
      sf::st_transform(x = _,
                       crs = sf::st_crs(polygons_extent)) |>
      sf::st_intersection(x = _,
                          y = polygons_extent)
  } else {
    warning("Unable to find a tblLPIHeader table in the geodatabase. No AIM points will be used.")
    input_data_list[["aim"]][["headers"]] <- NULL
  }
  if ("tblLPIHeader" %in% layers_lookup$layer_name & !is.null(input_data_list[["aim"]][["headers"]])) {
    input_data_list[["aim"]][["lpi_headers"]] <- sf::st_read(dsn = file.path(data_path,
                                                                             aim_geodatabase_filename),
                                                             layer = layers_lookup$internal_layer_name[layers_lookup$layer_name == "tblLPIHeader"],
                                                             quiet = TRUE) |>
      dplyr::filter(.data = _,
                    PrimaryKey %in% input_data_list[["aim"]][["headers"]][["PrimaryKey"]])
  } else {
    warning("Unable to find either a tblLPIHeader table or tblPlots in the geodatabase. No AIM points will be used.")
    input_data_list[["aim"]][["lpi_headers"]] <- NULL
  }
  if ("tblLPIDetail" %in% layers_lookup$layer_name & !is.null(input_data_list[["aim"]][["headers"]])) {
    input_data_list[["aim"]][["lpi_detail"]] <- sf::st_read(dsn = file.path(data_path,
                                                                            aim_geodatabase_filename),
                                                            layer = layers_lookup$internal_layer_name[layers_lookup$layer_name == "tblLPIDetail"],
                                                            quiet = TRUE) |>
      dplyr::filter(.data = _,
                    PrimaryKey %in% input_data_list[["aim"]][["headers"]][["PrimaryKey"]])
  } else {
    warning("Unable to find either a tblLPIDetail table or tblPlots in the geodatabase. No AIM points will be used.")
    input_data_list[["aim"]][["lpi_detail"]] <- NULL
  }
  if ("tblSpecRichHeader" %in% layers_lookup$layer_name & !is.null(input_data_list[["aim"]][["headers"]])) {
    input_data_list[["aim"]][["species_headers"]] <- sf::st_read(dsn = file.path(data_path,
                                                                                 aim_geodatabase_filename),
                                                                 layer = layers_lookup$internal_layer_name[layers_lookup$layer_name == "tblSpecRichHeader"],
                                                                 quiet = TRUE) |>
      dplyr::filter(.data = _,
                    PrimaryKey %in% input_data_list[["aim"]][["headers"]][["PrimaryKey"]])
  } else {
    warning("Unable to find either a tblSpecRichHeader table or tblPlots in the geodatabase. No AIM points will be used.")
    input_data_list[["aim"]][["species_headers"]] <- NULL
  }
  if ("tblSpecRichDetail" %in% layers_lookup$layer_name & !is.null(input_data_list[["aim"]][["headers"]])) {
    input_data_list[["aim"]][["species_detail"]] <- sf::st_read(dsn = file.path(data_path,
                                                                                aim_geodatabase_filename),
                                                                layer = layers_lookup$internal_layer_name[layers_lookup$layer_name == "tblSpecRichDetail"],
                                                                quiet = TRUE) |>
      dplyr::filter(.data = _,
                    PrimaryKey %in% input_data_list[["aim"]][["headers"]][["PrimaryKey"]])
  } else {
    warning("Unable to find either a tblSpecRichDetail table or tblPlots in the geodatabase. No AIM points will be used.")
    input_data_list[["aim"]][["species_detail"]] <- NULL
  }
  
  if (!is.null(input_data_list[["aim"]][["lpi_headers"]]) & !is.null(input_data_list[["aim"]][["lpi_detail"]])) {
    # Now we can convert the wide format into a long one!
    input_data_list[["lpi_tall"]][["aim"]] <- terradactyl::gather_lpi_terradat(tblLPIDetail = input_data_list[["aim"]][["lpi_detail"]],
                                                                               tblLPIHeader = input_data_list[["aim"]][["lpi_headers"]],
                                                                               verbose = TRUE)
  }
  if (!is.null(input_data_list[["aim"]][["species_headers"]]) & !is.null(input_data_list[["aim"]][["species_detail"]])) {
    # Now we can convert the wide format into a long one!
    input_data_list[["species_tall"]][["aim"]] <- terradactyl::gather_species_inventory_terradat(tblSpecRichDetail = input_data_list[["aim"]][["species_detail"]],
                                                                                                 tblSpecRichHeader = input_data_list[["aim"]][["species_headers"]],
                                                                                                 verbose = TRUE)
  }
  
  ####### LMF ------------------------------------------------------------------
  if ("POINTCOORDINATES" %in% layers_lookup$layer_name) {
    input_data_list[["lmf"]][["headers"]] <- sf::st_read(dsn = file.path(data_path,
                                                                         lmf_geodatabase_filename),
                                                         layer = layers_lookup$internal_layer_name[layers_lookup$layer_name == "POINTCOORDINATES"]) |>
      sf::st_transform(x = _,
                       crs = sf::st_crs(polygons_extent)) |>
      sf::st_intersection(x = _,
                          y = polygons_extent)
  } else {
    warning("Unable to find a POINTS layer in the geodatabase. No LMF points will be used.")
    input_data_list[["lmf"]][["headers"]] <- NULL
  }
  if ("PINTERCEPT" %in% layers_lookup$layer_name & !is.null(input_data_list[["lmf"]][["headers"]])) {
    input_data_list[["lmf"]][["lpi_detail"]] <- sf::st_read(dsn = file.path(data_path,
                                                                            lmf_geodatabase_filename),
                                                            layer = layers_lookup$internal_layer_name[layers_lookup$layer_name == "PINTERCEPT"],
                                                            quiet = TRUE) |>
      dplyr::filter(.data = _,
                    PrimaryKey %in% input_data_list[["lmf"]][["headers"]][["PrimaryKey"]])
  } else {
    warning("Unable to find either a PINTERCEPT table or POINTS feature in the geodatabase. No LMF points will be used.")
    input_data_list[["aim"]][["lpi_detail"]] <- NULL
  }
  if ("PLANTCENSUS" %in% layers_lookup$layer_name & !is.null(input_data_list[["lmf"]][["headers"]])) {
    input_data_list[["lmf"]][["species_detail"]] <- sf::st_read(dsn = file.path(data_path,
                                                                            lmf_geodatabase_filename),
                                                            layer = layers_lookup$internal_layer_name[layers_lookup$layer_name == "PLANTCENSUS"],
                                                            quiet = TRUE) |>
      dplyr::filter(.data = _,
                    PrimaryKey %in% input_data_list[["lmf"]][["headers"]][["PrimaryKey"]])
  } else {
    warning("Unable to find either a PLANTCENSUS table or POINTS feature in the geodatabase. No LMF points will be used.")
    input_data_list[["aim"]][["lpi_detail"]] <- NULL
  }
  
  if (!is.null(input_data_list[["lmf"]][["lpi_detail"]])) {
    # Now we can convert the wide format into a long one!
    input_data_list[["lpi_tall"]][["lmf"]] <- terradactyl::gather_lpi_lmf(PINTERCEPT = input_data_list[["lmf"]][["lpi_detail"]],
                                                                          verbose = TRUE)
  }
  if (!is.null(input_data_list[["lmf"]][["species_detail"]])) {
    # Now we can convert the wide format into a long one!
    input_data_list[["species_tall"]][["lmf"]] <- terradactyl::gather_species_inventory_lmf(PLANTCENSUS = input_data_list[["lmf"]][["species_detail"]],
                                                                          verbose = TRUE)
  }
  
  lpi_data <- dplyr::bind_rows(input_data_list[["lpi_tall"]])
  species_data <- dplyr::bind_rows(input_data_list[["species_tall"]]) |>
    dplyr::select(.data = _,
                  tidyselect::all_of(x = c("PrimaryKey",
                                           "Species")))
  headers <- dplyr::bind_rows(dplyr::select(.data = input_data_list[["aim"]][["headers"]],
                                            PrimaryKey),
                              dplyr::select(.data = input_data_list[["lmf"]][["headers"]],
                                            PrimaryKey))
}


##### Species information ------------------------------------------------------
# This isn't actualy necessary! We can calculate just based on the raw codes.


#### INDICATORS ################################################################
# This was written extensibly-ish, but we're only calculating the percentage
# of pinyon-juniper cover here.

##### Prepping data ------------------------------------------------------------
# Add in a pinyon variable
lpi_data <- dplyr::mutate(.data = lpi_data,
                          pinyon = dplyr::case_when(code %in% pinyon_species ~ "pinyon_cover",
                                                .default = NA))

pinyon_presence <- dplyr::mutate(.data = species_data,
                                 pinyon = Species %in% pinyon_species) |>
  dplyr::summarize(.data = _,
                   .by = tidyselect::all_of(x = c("PrimaryKey")),
                   pinyon_present = any(pinyon))

##### Calculating indicators ---------------------------------------------------
# We'll calculate any hit cover (a pin drop will count towards cover if any
# layer record at that location is associated with a pinyon species.)

# Note that if a plot had a record for pinyon on the plot but a cover of 0, that
# gets treated as having 0.1% cover, but it's likely more.
pinyon_cover <- terradactyl::pct_cover(lpi_tall = lpi_data,
                                   tall = FALSE,
                                   hit = "any",
                                   by_line = FALSE,
                                   indicator_variables = c("pinyon"),
                                   digits = 1,
                                   verbose = TRUE) |>
  dplyr::left_join(x = _,
                   y = pinyon_presence,
                   by = "PrimaryKey",
                   relationship = "one-to-one") |>
  dplyr::mutate(.data = _,
                pinyon_cover = dplyr::case_when(pinyon_cover == 0 & pinyon_present ~ 0.1,
                                                .default = pinyon_cover))

##### Classifying by pinyon cover --------------------------------------------------
# We want to assign the classification based on the cover percentage using the
# list defined above in the parameters.
# Sloppy to do this as a for loop, but we'll live.

for (current_category in names(pinyon_cover_classes)) {
  pinyon_cover <- dplyr::mutate(.data = pinyon_cover,
                            # We stated above that the min is inclusive but the
                            # max is not!
                            meets_min = pinyon_cover >= pinyon_cover_classes[[current_category]]["min"],
                            meets_max = pinyon_cover < pinyon_cover_classes[[current_category]]["max"],
                            qualifies = meets_min & meets_max)
  pinyon_cover[["cover_category"]][pinyon_cover$qualifies] <- current_category
  pinyon_cover <- dplyr::select(.data = pinyon_cover,
                            -tidyselect::all_of(x = c("meets_min",
                                                      "meets_max",
                                                      "qualifies")))
}


#### WEIGHTING #################################################################
# This is the standard workflow.
# For each polygon (frame) that we want a weighted analysis output for we'll:
# 1) Get the subset of points that occur within the frame.
# 2) Calculate the density of sampling across the frame.
# 3) Break the frame into density partitions.
# 4) Generate Thiessen polygons within each partition.
# 5) Use the Thiessen polygons to weight each point.

# This uses a for loop instead of a lapply() because I want to make sure that we
# don't lose stuff that's already been calculated if there's a memory allocation
# issue or something similar with a polygon.

##### TP/DP parameters ---------------------------------------------------------
## Thiessen polygons and density partitions ---
# These parameters are used to control the use of Thiessen polygons and density
# partitioning to find estimates. The script will attempt to generate estimates
# using every requested configuration, so if you specify 2, 3, and 4 Thiessen
# polygons you'll get three different estimates per frame (one for 2 TPs, one
# for 3 TPs, and one for 4 TPs).

# The minimum number of data points that must fall within a Thiessen polygon in
# order for it to be valid. Higher numbers will result in longer processing
# times and may make some estimates impossible, e.g., if you try to use 3
# Thiessen polygons with a minimum of 5 points per polygon, the script will not
# attempt to make estimates for a frame that contains fewer than 15 points.
min_points_per_tpoly <- 1

# Thiessen polygons are randomized, but we would like to make these estimates
# repeatable. This is the set of numeric values to use as seed numbers which
# set the random state at the time of Thiessen polygon generation.
# This also controls how many Thiessen polygon solutions the script will attempt
# to find. If you want to find 100 Thiessen polygons solutions for each estimate
# you would need to have 100 seed numbers here.
# If the script can't find a valid solution within an iteration limit, you may
# have fewer solutions going into the final estimate than you requested. This is
# more common with very uneven spatial distributions for the data and low sample
# sizes.
tpoly_seeds <- c(1:25)

# We're going to partition the frames by point density. We'll try with a couple
# levels of partitioning, so this tells us how many levels to try.
# For each number in the vector there will be one estimate using only the
# density partitions and additional estimates where Thiessen polygons were drawn
# within them according to dp_tpoly_counts below.
density_partition_counts <- 3

# This is the same as tpoly_counts above except that this is the number of
# Thiessen polygons that will be drawn within the density partitions.
# tpoly_counts controls only estimates made without using density partitioning.
dp_tpoly_counts <- c(3)

# This is the sampling distance in meters for calculating density.
# The smaller this number, the finer-grained the density map but that can cause
# serious memory allocation issues when you have larger areas.
# Don't worry about changing this unless you're analyzing very large areas and
# are running into time or memory constraints.
density_sample_spacing <- 500

# In some cases, the frame polygons may be buffered for the sake of calculating
# density across them. This is used so that the density raster extends very
# slightly beyond the edge of the frame in order to completely cover the frame.
buffer_distance <- min(100,
                       density_sample_spacing / 2)


## Bootstrapping ---
# Bootstrapping can be used to combine multiple valid Thiessen polygon solutions
# into a single estimate. It will NOT be applied to situations with only one
# valid solution, e.g., equal weighting or density partitions without Thiessen
# polygons.

# bootstrap_count is the number of bootstrap replicates to use. Generally, more
# is "better" because it's sort of like bumping up the sample size, but that
# won't always hold true.
bootstrap_count <- 10000

# The bootstrapping is done using the package boot. bootstrap_types is passed to
# boot::boot_ci() to determine what kind(s) of bootstrapping should be used to
# calculate the confidence intervals. See the documentation for that function
# for details. If you request multiple types, the output will provide each
# separately.
bootstrap_types <- c("bca")

##### Finding weighting polygons -----------------------------------------------
# We refer to these internally as "weight categories" because "poststrata" has
# been used to describe using the design polygons for weighting in AIM.
point_weights_list <- list()
points_attributed_dp_tpoly_list <- list()
for (current_frame in polygons_list) {
  current_points <- sf::st_intersection(x = headers,
                                        y = current_frame)
  
  ###### Density partitions ####################################################
  # This next stretch is basically a vivisected density_polygon_gen_clustered()
  # and it's just finding the requested number of density partitions for this
  # set of analyses.
  
  density_partitions_list <- lapply(X = density_partition_counts,
                                    current_frame_id = current_frame$uid[1],
                                    frame = current_frame,
                                    points = current_points,
                                    buffer_distance = buffer_distance,
                                    density_sample_spacing = density_sample_spacing,
                                    FUN = function(X, current_frame_id, frame, points, buffer_distance, density_sample_spacing){
                                      intended_partition_count <- X
                                      current_partition_count <- X
                                      # Get the point coordinates to feed into spatstat.geom::ppp()
                                      point_coords <- unique(sf::st_coordinates(points))
                                      
                                      # Get the owin object for spatstat.geom::ppp()
                                      # We're going to get weird here. The polygons are often too complicated to
                                      # generate an owin object from, so we're going to remove holes and buffer a bit
                                      # to get simpler polygons that are slightly larger than the inference area so
                                      # we can calculate density then trim those results down to the inference area.
                                      density_frame <- sf::st_buffer(x = frame,
                                                                     dist = buffer_distance)
                                      # density_frame <- nngeo::st_remove_holes(density_frame,
                                      #                                         max_area = 20^2)
                                      frame_owin <- spatstat.geom::as.owin(density_frame)
                                      # frame_owin <- spatstat.geom::as.owin(inference_area)
                                      
                                      
                                      message("Creating point pattern object")
                                      # Make the ppp object (a point pattern)
                                      points_ppp <- spatstat.geom::ppp(x = point_coords[, 1],
                                                                       y = point_coords[, 2],
                                                                       window = frame_owin)
                                      
                                      message("Finding distribution density from point pattern")
                                      # Get the density info from the point pattern
                                      points_density <- density(points_ppp,
                                                                # For the CRS we're using, the units on this will be meters
                                                                # But it's CRS-dependent
                                                                eps = density_sample_spacing)
                                      
                                      # Make a data frame of coordinates with the density at each coordinate
                                      density_df <- expand.grid(y = points_density$yrow,
                                                                x = points_density$xcol)
                                      density_df$density <- as.vector(points_density$v)
                                      
                                      
                                      # # And now that we have the density data frame
                                      # current_density_partitions_list <- list()
                                      # for (current_partition_count in density_partition_counts) {
                                      message("Finding partition breaks. This can take a while.")
                                      # Figure out where the breaks are for the partitions
                                      partition_breaks <- BAMMtools::getJenksBreaks(var = density_df[["density"]],
                                                                                    # This'll a number of values equal to k
                                                                                    # and the terminal values will be the min and max
                                                                                    # So in order to get partition ranges, we need to add
                                                                                    # 1 to k so we get enough breakpoints
                                                                                    k = current_partition_count + 1)
                                      
                                      message("Classifying area by densities")
                                      for (partition_id in length(partition_breaks):2) {
                                        message(paste0("Identifying partition ", partition_id - 1))
                                        # Get the upper and lower cutoff values
                                        upper <- partition_breaks[partition_id]
                                        lower <- partition_breaks[partition_id - 1]
                                        
                                        # Determine if a value is below the current upper
                                        # bound and above the current lower bound for
                                        # the quantile
                                        below_upper <- sapply(X = density_df$density,
                                                              upper = upper,
                                                              FUN = function(X, upper){
                                                                if (is.na(X)) {
                                                                  FALSE
                                                                } else {
                                                                  X <= upper
                                                                }
                                                              })
                                        above_lower <- sapply(X = density_df$density,
                                                              lower = lower,
                                                              FUN = function(X, lower){
                                                                if (is.na(X)) {
                                                                  FALSE
                                                                } else {
                                                                  X >= lower
                                                                }
                                                              })
                                        
                                        applicable_indices <- mapply(X = below_upper,
                                                                     Y = above_lower,
                                                                     FUN = function(X, Y){
                                                                       X & Y
                                                                     })
                                        
                                        # Write in the current quantile ID to the relevant indices
                                        density_df$partition_id[applicable_indices] <- partition_id - 1
                                      }
                                      
                                      # Get a stars object, which is basically a raster
                                      density_stars <- stars::st_as_stars(.x = density_df,
                                                                          # This defaults to 1:2, but that transposes
                                                                          # the x and y axes, so we'll do 2:1
                                                                          coords = 2:1)
                                      
                                      # Convert the stars object to polygons
                                      # We have to use quantile_id instead of quantile because this'll only work with
                                      # numeric values. We can always get quantiles in there later with a join
                                      density_sf <- sf::st_as_sf(density_stars["partition_id"],
                                                                 as_points = FALSE,
                                                                 merge = TRUE)
                                      # Make sure that the CRS is assigned
                                      sf::st_crs(density_sf) <- sf::st_crs(frame)
                                      
                                      # OKAY!!!!!!!!!
                                      # So here's where we do a while() thing if the user wants to make sure there are
                                      # points in every polygon
                                      message("Carrying out point check. This involves a spatial join, so it can be slow.")
                                      partition_ids <- sf::st_drop_geometry(density_sf)$partition_id
                                      partitions_with_points <- unique(sf::st_drop_geometry(sf::st_intersection(x = points,
                                                                                                                y = density_sf[, "partition_id"]))[["partition_id"]])
                                      
                                      # And while not all partition IDs are represented, try again
                                      while (!all(partition_ids %in% partitions_with_points)) {
                                        warning("Not all density partitions contained points. Attempting to with an additional partition with empty partitions being merged up.")
                                        current_partition_count <- current_partition_count + 1
                                        
                                        if (current_partition_count > 1.5 * intended_partition_count) {
                                          warning("Hit the density partition limit without finding enough. Returning the frame as if it were a single density partition.")
                                          density_sf <- dplyr::mutate(.data = frame,
                                                                      partition_id = 1) |>
                                            dplyr::select(.data = _,
                                                          partition_id)
                                          break()
                                        }
                                        
                                        message(paste("Starting new attempt with", current_partition_count, "partitions. Don't worry: the final density partition count will be the actual requested one. These will be rolled up as necessary."))
                                        
                                        # Figure out where the breaks are for the partitions
                                        partition_breaks <-BAMMtools::getJenksBreaks(var = density_df[["density"]],
                                                                                     # This'll a number of values equal to k
                                                                                     # and the terminal values will be the min and max
                                                                                     # So in order to get partition ranges, we need to add
                                                                                     # 1 to k so we get enough breakpoints
                                                                                     k = current_partition_count + 1)
                                        
                                        
                                        message("Classifying area by densities")
                                        
                                        for (partition_id in length(partition_breaks):2) {
                                          # Get the upper and lower cutoff values
                                          upper <- partition_breaks[partition_id]
                                          lower <- partition_breaks[partition_id - 1]
                                          
                                          # Determine if a value is below the current upper
                                          # bound and above the current lower bound for
                                          # the quantile
                                          below_upper <- sapply(X = density_df$density,
                                                                upper = upper,
                                                                FUN = function(X, upper){
                                                                  if (is.na(X)) {
                                                                    FALSE
                                                                  } else {
                                                                    X <= upper
                                                                  }
                                                                })
                                          above_lower <- sapply(X = density_df$density,
                                                                lower = lower,
                                                                FUN = function(X, lower){
                                                                  if (is.na(X)) {
                                                                    FALSE
                                                                  } else {
                                                                    X >= lower
                                                                  }
                                                                })
                                          
                                          applicable_indices <- mapply(X = below_upper,
                                                                       Y = above_lower,
                                                                       FUN = function(X, Y){
                                                                         X & Y
                                                                       })
                                          
                                          # Write in the current quantile ID to the relevant indices
                                          density_df$partition_id[applicable_indices] <- partition_id - 1
                                        }
                                        
                                        
                                        # Get a stars object, which is basically a raster
                                        density_stars <- stars::st_as_stars(.x = density_df,
                                                                            # This defaults to 1:2, but that transposes
                                                                            # the x and y axes, so we'll do 2:1
                                                                            coords = 2:1)
                                        
                                        # Convert the stars object to polygons
                                        # We have to use quantile_id instead of quantile because this'll only work with
                                        # numeric values. We can always get quantiles in there later with a join
                                        density_sf <- sf::st_as_sf(density_stars["partition_id"],
                                                                   as_points = FALSE,
                                                                   merge = TRUE)
                                        # Make sure that the CRS is assigned
                                        sf::st_crs(density_sf) <- sf::st_crs(frame)
                                        
                                        # And now so we can combine the polygons
                                        partition_ids <- sf::st_drop_geometry(density_sf)[["partition_id"]]
                                        partitions_with_points <- unique(sf::st_drop_geometry(sf::st_intersection(x = points,
                                                                                                                  y = density_sf[, "partition_id"]))[["partition_id"]])
                                        partitions_without_points <- partition_ids[!(partition_ids %in% partitions_with_points)]
                                        
                                        # This needs to be ordered so that they roll up (e.g., the lowest empty
                                        # partition gets added to the next up THEN that partition which may've been
                                        # empty itself can be added up as well to avoid infinite loops)
                                        for (empty_partition in partitions_without_points[order(partitions_without_points, decreasing = FALSE)]) {
                                          density_sf$partition_id[density_sf$partition_id == empty_partition] <- density_sf$partition_id[density_sf$partition_id == empty_partition] + 1
                                        }
                                        
                                        partition_ids <- unique(sf::st_drop_geometry(density_sf)[["partition_id"]])
                                        # UPDATE THE PARTITION IDS TO START AT 1 AND INCREMENT BY 1
                                        # MERGE THE PARTITIONS BY ID
                                        current_density_partition_list <- lapply(X = 1:length(partition_ids),
                                                                                 current_ids = partition_ids[order(partition_ids)],
                                                                                 density_polygons = density_sf,
                                                                                 FUN = function(X, current_ids, density_polygons){
                                                                                   current_polygons <- density_polygons[density_polygons$partition_id == current_ids[X], ]
                                                                                   
                                                                                   current_polygons <- sf::st_as_sf(x = sf::st_union(x = current_polygons))
                                                                                   
                                                                                   current_polygons[["partition_id"]] <- X
                                                                                   
                                                                                   current_polygons
                                                                                 })
                                        
                                        density_sf <- do.call(rbind,
                                                              current_density_partition_list)
                                        
                                        # And now for the while() loop's benefit
                                        partition_ids <- sf::st_drop_geometry(density_sf)[["partition_id"]]
                                        partitions_with_points <- unique(sf::st_drop_geometry(sf::st_intersection(x = points,
                                                                                                                  y = density_sf[, "partition_id"]))[["partition_id"]])
                                        
                                      }
                                      density_sf[["analysis_frame_id"]] <- current_frame_id
                                      # current_density_partitions_list[[paste0(X, "_",
                                      #                                         length(unique(density_sf$partition_id)))]] <- density_sf
                                      # }
                                      # current_density_partitions_list
                                      
                                      # This should return NULL only if the density partitions are real messed up
                                      if (!(length(unique(density_sf$partition_id)) %in% c(1, intended_partition_count))) {
                                        NULL
                                      } else {
                                        density_sf
                                      }
                                    })
  
  names(density_partitions_list) <- paste0(current_frame$uid, "_",
                                           density_partition_counts)
  density_partitions_list <- density_partitions_list[!sapply(X = density_partitions_list,
                                                             FUN = is.null)]
  
  ###### Within-partition Thiessen polygons -------------------------------------
  # This is a loop instead of a lapply() so that we get to keep solutions even
  # if one of them fails.
  dp_tpoly_list <- list()
  
  # We'll go through this one level of partitioning at a time
  # Get ready for a gnarly nested list situation because it's the most obvious solution
  # to help keep things organized
  for (tpoly_count in dp_tpoly_counts) {
    message(paste0("Beginning to find tpoly solutions with ", tpoly_count, " polygons per density partition"))
    current_count_tpoly_list <- list()
    
    for (tpoly_seed_index in 1:length(tpoly_seeds)) {
      message(paste0("Starting with seed ", tpoly_seed_index, " of ", length(tpoly_seeds)))
      current_tpoly_seed <- tpoly_seeds[tpoly_seed_index]
      
      current_seed_tpoly_list <- list()
      
      for (density_partitions_name in names(density_partitions_list)) {
        message(paste0("Working with ", density_partitions_name, " density partitions"))
        density_partition_count <- as.numeric(stringr::str_extract(string = density_partitions_name,
                                                                   pattern = "\\d+$"))
        
        current_density_partitions <- density_partitions_list[[density_partitions_name]]
        
        frame_name <- current_density_partitions$analysis_frame_id[1]
        
        current_partition_tpoly_list <- list()
        for (partition in unique(current_density_partitions$partition_id)) {
          previous_partitions <- setdiff(x = unique(current_density_partitions$partition_id)[1:which(unique(current_density_partitions$partition_id) == partition)],
                                         partition)
          has_failed_tpoly_solution <- !all(previous_partitions %in% names(current_partition_tpoly_list))
          
          if (!has_failed_tpoly_solution) {
            message(paste0("Solving for partition ", which(unique(current_density_partitions$partition_id) == partition), " of ", length(unique(current_density_partitions$partition_id))))
            current_density_partition <- dplyr::filter(.data = current_density_partitions,
                                                       partition_id == partition)
            
            points_available <- sf::st_join(x = headers,
                                            y = current_density_partition) |>
              dplyr::filter(.data = _,
                            !is.na(partition_id))
            if (nrow(points_available) < min_points_per_tpoly * tpoly_count) {
              message("Inadequate number of points to meet minimum points-per-polygon requirement. Skipping.")
              current_tpolys <- NULL
            } else {
              current_tpolys <- thiessen_polygons_gen_random(frame = current_density_partition,
                                                             n_polygons = tpoly_count,
                                                             points = headers,
                                                             points_min = min_points_per_tpoly,
                                                             envelope = sf::st_union(current_density_partition),
                                                             seed_number = current_tpoly_seed,
                                                             seed_increment = 1000,
                                                             iteration_limit = 100,
                                                             use_albers = TRUE,
                                                             verbose = TRUE)
              
              if (is.null(current_tpolys)) {
                message("Unable to find a solution for a density partition in this set. This is either due to an insufficient number of points or hitting the iteration limit.")
              } else {
                current_tpolys <- dplyr::mutate(.data = current_tpolys,
                                                analysis_frame_id = current_density_partition[["analysis_frame_id"]][1],
                                                partition_id = partition,
                                                tpoly_id = paste0("partition_", partition, "_", length(unique(current_density_partitions$partition_id)), "-",
                                                                  tpoly_id),
                                                partition_count = density_partition_count,
                                                tpoly_count = tpoly_count,
                                                seed_id = tpoly_seed_index)
              }
            }
          } else {
            message("Unable to find a solution for a density partition in this set. This is either due to an insufficient number of points or hitting the iteration limit.")
            current_tpolys <- dplyr::select(.data = current_density_partition,
                                            tidyselect::all_of(x = c("analysis_frame_id"))) |>
              dplyr::mutate(.data = _,
                            partition_id = partition,
                            polygon_unique_id = 1,
                            tpoly_id = paste0("partition_", partition, "_", length(unique(current_density_partitions$partition_id)), "-",
                                              "1"),
                            partition_count = density_partition_count,
                            tpoly_count = 1,
                            seed_id = tpoly_seed_index)
          }
          
          current_partition_tpoly_list[[paste0(partition)]] <- current_tpolys
        }
        message(paste0("Completed attempt to find solution ", tpoly_seed_index, " of ", length(tpoly_seeds), " for ",
                       density_partition_count, " density partitions with ",
                       tpoly_count, " tpolys per partition"))
        
        
        if (!all(unique(current_density_partitions$partition_id) %in% names(current_partition_tpoly_list))) {
          message("At least one density partition does not have a valid Thiessen polygon solution due to insufficient points or hitting the iteration limit. Any polygons found for other density partitions for this pass will be dropped and this solution 'skipped'.")
          current_density_tpolys <- NULL
        } else {
          current_density_tpolys <- dplyr::bind_rows(current_partition_tpoly_list)
        }
        
        dp_tpoly_list[[paste0(frame_name, "_frame-",
                              tpoly_seed_index, "_solution-",
                              density_partition_count, "_partitions-",
                              tpoly_count, "_tpolys")]] <- current_density_tpolys
      }
    }
  }
  
  expected_tpoly_list_length <- length(density_partitions_list) * length(dp_tpoly_counts) * length(tpoly_seeds)
  if (expected_tpoly_list_length > length(dp_tpoly_list)) {
    warning(paste0("There were only ", length(dp_tpoly_list), " Thiessen polygon solutions found but you expected ",
                   expected_tpoly_list_length, ". This is almost certainly due to low point counts in one or more density partitions making Thiessen polygons impossible with the required minimum point count or simply so unlikely that the iteration limit was hit."))
  }
  
  ##### Calculating weights -------------------------------------------------
  points_attributed_dp_tpoly_list[[current_frame$uid[1]]] <- lapply(X = names(dp_tpoly_list),
                                                                    dp_tpoly_list = dp_tpoly_list,
                                                                    points_relevant = dplyr::select(headers,
                                                                                                    PrimaryKey),
                                                                    FUN = function(X, dp_tpoly_list, points_relevant){
                                                                      message(paste0("Attributing from tpoly set ", which(names(dp_tpoly_list) == X), " of ", length(dp_tpoly_list)))
                                                                      parameter_vector <- stringr::str_split_1(string = X,
                                                                                                               pattern = "-")
                                                                      names(parameter_vector) <- stringr::str_extract(string = parameter_vector,
                                                                                                                      pattern = "[a-z]+$")
                                                                      parameter_vector <- gsub(x = parameter_vector,
                                                                                               pattern = "_[a-z]+$",
                                                                                               replacement = "")
                                                                      
                                                                      current_tpolys <- dp_tpoly_list[[X]]
                                                                      current_tpolys$area <- as.numeric(sf::st_area(current_tpolys))
                                                                      
                                                                      current_points_attributed <- sf::st_join(x = points_relevant,
                                                                                                               y = dplyr::select(sf::st_transform(current_tpolys,
                                                                                                                                                  crs = sf::st_crs(points_relevant)),
                                                                                                                                 wgtcat_id = tpoly_id),
                                                                                                               left = FALSE)
                                                                      current_weight_summary <- dplyr::summarize(.data = dplyr::group_by(.data = sf::st_drop_geometry(current_points_attributed),
                                                                                                                                         wgtcat_id),
                                                                                                                 n_points = dplyr::n())
                                                                      current_weight_summary <- dplyr::left_join(x = current_weight_summary,
                                                                                                                 y = dplyr::select(sf::st_drop_geometry(current_tpolys),
                                                                                                                                   wgtcat_id = tpoly_id,
                                                                                                                                   area),
                                                                                                                 by = "wgtcat_id")
                                                                      current_weight_summary <- dplyr::mutate(.data = current_weight_summary,
                                                                                                              weight = area / n_points)
                                                                      current_points_attributed <- dplyr::left_join(x = current_points_attributed,
                                                                                                                    y = dplyr::select(current_weight_summary,
                                                                                                                                      wgtcat_id,
                                                                                                                                      weight),
                                                                                                                    by = "wgtcat_id")
                                                                      
                                                                      current_points_attributed[["analysis_frame_id"]] <- current_tpolys[["analysis_frame_id"]][1]
                                                                      current_points_attributed[["weighting"]] <- "Thiessen polygons within density partitions"
                                                                      for (parameter in names(parameter_vector)) {
                                                                        current_points_attributed[[parameter]] <- as.numeric(parameter_vector[parameter])
                                                                      }
                                                                      sf::st_drop_geometry(current_points_attributed)
                                                                    })
}

#### ANALYZING #################################################################
# With the weights in hand, we can run the analyses!
##### Analyze for each set of weights --------------------------------------
# This is for "each set of weights" because it's adapted from code that was
# originally written to compare multiple weighting approaches.
weighting_list <- c(unlist(points_attributed_dp_tpoly_list,
                           recursive = FALSE))

analysis_list <- lapply(X = weighting_list,
                        points = sf::st_drop_geometry(pinyon_cover) |>
                          tidyr::pivot_longer(data = _,
                                              cols = tidyselect::all_of(x = c(continuous_indicators,
                                                                              categorical_indicators)),
                                              names_to = "indicator",
                                              values_to = "value",
                                              values_transform = as.character) |>
                          dplyr::select(.data = _,
                                        PrimaryKey,
                                        indicator,
                                        value) |>
                          dplyr::mutate(.data = _,
                                        value = tidyr::replace_na(data = value,
                                                                  replace = "0")),
                        conf = conf,
                        continuous_indicators = continuous_indicators,
                        categorical_indicators = categorical_indicators,
                        FUN = function(X, points, conf, continuous_indicators, categorical_indicators){
                          message(paste0("Current weights are for the ", stringr::str_replace_all(string = X[["analysis_frame_id"]][1],
                                                                                                  pattern = "_",
                                                                                                  replacement = " "),
                                         " using ", tolower(X[["weighting"]][1]),
                                         " with ", X[["partitions"]][1], " partitions and ",
                                         X[["tpolys"]][1], " tpolys."))
                          points_weighted <- dplyr::inner_join(x = points,
                                                               y = X,
                                                               by = "PrimaryKey",
                                                               relationship = "many-to-one")
                          
                          per_indicator_analyses_list <- lapply(X = split(x = points_weighted,
                                                                          f = points_weighted$indicator),
                                                                conf = conf,
                                                                FUN = function(X, conf){
                                                                  # Dunno why I need this, but something's weird with the join or split
                                                                  X <- dplyr::distinct(X)
                                                                  
                                                                  current_indicator <- X[["indicator"]][1]
                                                                  current_weighting <- X[["weighting"]][1]
                                                                  current_solution_id <- X[["solution"]][1]
                                                                  current_partition_count <- X[["partitions"]][1]
                                                                  current_tpoly_count <- X[["tpolys"]][1]
                                                                  current_frame <- X[["analysis_frame_id"]][1]
                                                                  
                                                                  message(paste0("Working with data for ", current_indicator,
                                                                                 " weighted using ", current_weighting, " with ",
                                                                                 current_partition_count, " partitions and ",
                                                                                 current_tpoly_count, " tpolys. This is solution ",
                                                                                 current_solution_id, "."))
                                                                  
                                                                  current_weights <- dplyr::select(.data = X,
                                                                                                   PrimaryKey,
                                                                                                   weight) |>
                                                                    dplyr::distinct()
                                                                  
                                                                  if (nrow(X) == 0) {
                                                                    return(NULL)
                                                                  } else if (nrow(X) == 1) {
                                                                    analysis <- data.frame(n = nrow(X),
                                                                                           alpha = 1 - conf / 100,
                                                                                           estimate = X$value[1],
                                                                                           sd = 0,
                                                                                           variance = 0,
                                                                                           lower_bound = NA,
                                                                                           upper_bound = NA)
                                                                  } else {
                                                                    # Make sure we calculate using the correct method
                                                                    if (current_indicator %in% continuous_indicators) {
                                                                      analysis <- dplyr::mutate(.data = X,
                                                                                                value = as.numeric(value)) |>
                                                                        analyze_con(data = _,
                                                                                    weights = current_weights,
                                                                                    id_var = "PrimaryKey",
                                                                                    value_var = "value",
                                                                                    wgt_var = "weight",
                                                                                    conf = conf,
                                                                                    verbose = TRUE) |>
                                                                        dplyr::rename(.data = _,
                                                                                      estimate = mean)
                                                                    } else if (current_indicator %in% categorical_indicators) {
                                                                      analysis <- analyze_cat(data = X,
                                                                                              weights = current_weights,
                                                                                              id_var = "PrimaryKey",
                                                                                              cat_var = "value",
                                                                                              wgt_var = "weight",
                                                                                              definitions = c("No pinyon",
                                                                                                              "1-9%",
                                                                                                              "Over 9%"),
                                                                                              conf = conf,
                                                                                              verbose = TRUE) |>
                                                                        dplyr::rename(.data = _,
                                                                                      n = observation_count,
                                                                                      estimate = weighted_observation_proportion,
                                                                                      sd = tidyselect::matches(match = "standard_error"),
                                                                                      cv = tidyselect::matches(match = "coefficient_of_variance"),
                                                                                      lower_bound = tidyselect::matches(match = "lower_bound"),
                                                                                      upper_bound = tidyselect::matches(match = "upper_bound")) |>
                                                                        dplyr::mutate(.data = _,
                                                                                      n = sum(n)) |>
                                                                        dplyr::select(.data = _,
                                                                                      -tidyselect::matches(match = "observation"))
                                                                      
                                                                      # current_indicator <- paste0(current_indicator, "_", analysis$category)
                                                                      
                                                                    } else {
                                                                      stop(paste(current_indicator, "is not currently in vectors categorical_indicators or continuous_indicators and must be assigned to one."))
                                                                    }
                                                                  }
                                                                  
                                                                  analysis$indicator <- current_indicator
                                                                  analysis$weighting <- current_weighting
                                                                  analysis$frame <- current_frame
                                                                  analysis$n_tpolys <- current_tpoly_count
                                                                  analysis$n_partitions <- current_partition_count
                                                                  analysis$solution_id <- current_solution_id
                                                                  analysis$alpha <- (100 - conf) / 100
                                                                  
                                                                  analysis
                                                                })
                          
                          dplyr::bind_rows(per_indicator_analyses_list)
                        })


analysis_names <- sapply(X = analysis_list,
                         FUN = function(X){
                           current_indicator <- X[["indicator"]][1]
                           current_weighting <- X[["weighting"]][1]
                           current_solution_id <- X[["solution_id"]][1]
                           current_partition_count <- X[["n_partitions"]][1]
                           current_tpoly_count <- X[["n_tpolys"]][1]
                           current_alpha <- X[["alpha"]][1]
                           current_frame <- X[["frame"]][1]
                           paste0(current_indicator, "-",
                                  current_frame, "_frame-",
                                  current_solution_id, "_solution-",
                                  current_weighting, "_weighting-",
                                  current_partition_count, "_partitions-",
                                  current_tpoly_count, "_tpolys-",
                                  current_alpha, "_alpha")
                         })
names(analysis_list) <- analysis_names

# Bind everything together!
analysis_results <- dplyr::bind_rows(analysis_list)

##### Combine tpoly results from equivalent solutions ----------------------
# Let's get those tpoly results combined to get a single estimate with a
# confidence interval
tpoly_analyses <- dplyr::filter(.data = analysis_results,
                                # n_tpolys > 0,
                                grepl(weighting,
                                      pattern = "polygons")) |>
  dplyr::mutate(.data = _,
                indicator = dplyr::case_when(!is.na(category) ~ category,
                                             .default = indicator))

tpoly_summaries_list <- lapply(X = split(x = tpoly_analyses,
                                         f = list(tpoly_analyses$indicator,
                                                  tpoly_analyses$frame,
                                                  tpoly_analyses$n_tpolys,
                                                  tpoly_analyses$n_partitions),
                                         drop = TRUE),
                               bootstrap_types = bootstrap_types,
                               bootstrap_count = bootstrap_count,
                               FUN = function(X, bootstrap_types, bootstrap_count){
                                 message(paste0(X$frame) |>
                                           unique())
                                 # message(X$indicator[1])
                                 
                                 # If there's only one solution, we'll skip the
                                 # combination stage.
                                 if (length(unique(X$solution_id)) > 1) {
                                   # If there're variances available, we'll
                                   # go ahead and combine without
                                   # bootstrapping but if we don't have
                                   # variance values (e.g., the estimates
                                   # are categorical) we'll use bootstrapping
                                   if (!any(is.na(X$variance))) {
                                     # nrow(X) is the number of records in the
                                     # data frame which is also the number of
                                     # solutions being combined.
                                     mean_of_means <- mean(X$estimate)
                                     # Each record already has the variance
                                     # calculated so we get to just sum them
                                     # without having to calculate them again.
                                     mean_of_means_variance <- sum(X$variance) / nrow(X)^2
                                     mean_of_means_standard_deviation <- sqrt(mean_of_means_variance)
                                     mean_of_means_standard_error <- mean_of_means_standard_deviation / sqrt(nrow(X))
                                     ci_bounds <- ci_mean(mean = mean_of_means,
                                                          sd = mean_of_means_standard_deviation,
                                                          n = nrow(X),
                                                          alpha = 1 - conf / 100)
                                     
                                     output <- data.frame(indicator = X$indicator[1],
                                                          frame = X$frame[1],
                                                          n_tpolys = X$n_tpolys[1],
                                                          n_partitions = X$n_partitions[1],
                                                          mean_estimate = mean_of_means,
                                                          alpha = 1 - conf / 100,
                                                          n_points = X$n[1],
                                                          n_solutions = nrow(X),
                                                          lower_bound = ci_bounds[["lower_bound"]],
                                                          upper_bound = ci_bounds[["upper_bound"]],
                                                          mean_estimate_variance = mean_of_means_variance,
                                                          mean_estimate_sd = mean_of_means_standard_deviation,
                                                          mean_estimate_se = mean_of_means_standard_error)
                                   } else {
                                     bootstrap_results <- boot::boot(data = X$estimate,
                                                                     # The function special_mean()
                                                                     # is just mean() but with an
                                                                     # added index argument so that
                                                                     # the data are subset appropriately
                                                                     # for each bootstrap replicate.
                                                                     statistic = special_mean,
                                                                     R = bootstrap_count)
                                     
                                     output <- data.frame(indicator = X$indicator[1],
                                                          frame = X$frame[1],
                                                          n_tpolys = X$n_tpolys[1],
                                                          n_partitions = X$n_partitions[1],
                                                          mean_estimate = bootstrap_results$t[1],
                                                          alpha = 1 - conf / 100,
                                                          n_points = X$n[1],
                                                          n_solutions = nrow(X),
                                                          booststrap_replicates = bootstrap_count,
                                                          lower_bound = bootstrap_results$t[1],
                                                          upper_bound = bootstrap_results$t[1],
                                                          ci_bootstrap_type = "none")
                                     if (length(unique(bootstrap_results$t)) != 1) {
                                       bootstrap_cis <- boot::boot.ci(boot.out = bootstrap_results,
                                                                      conf = conf / 100,
                                                                      type = bootstrap_types)
                                       output <- lapply(X = setdiff(x = names(bootstrap_cis),
                                                                    y = c("R", "t0", "call")),
                                                        bootstrap_cis = bootstrap_cis,
                                                        output = output,
                                                        FUN = function(X, bootstrap_cis, output){
                                                          bounds <- bootstrap_cis[[X]][(length(bootstrap_cis[[X]]) - 1):length(bootstrap_cis[[X]])]
                                                          dplyr::mutate(.data = output,
                                                                        lower_bound = bounds[1],
                                                                        upper_bound = bounds[2],
                                                                        ci_bootstrap_type = X)
                                                        }) |>
                                         dplyr::bind_rows()
                                     }
                                   }
                                 } else {
                                   output <- dplyr::select(.data = X,
                                                           tidyselect::all_of(x = c("indicator",
                                                                                    "frame",
                                                                                    "n_tpolys",
                                                                                    "n_partitions",
                                                                                    "mean_estimate" = "estimate",
                                                                                    "n_points" = "n",
                                                                                    "alpha",
                                                                                    "lower_bound",
                                                                                    "upper_bound"))) |>
                                     dplyr::mutate(.data = _,
                                                   n_solutions = 1,
                                                   bootstrap_replicates = 0,
                                                   ci_bootstrap_type = "None")
                                 }
                                 
                                 
                                 # Reorder!
                                 dplyr::select(.data = output,
                                               tidyselect::any_of(x = c("indicator",
                                                                        "frame",
                                                                        "n_tpolys",
                                                                        "n_partitions",
                                                                        "mean_estimate", 
                                                                        "alpha",
                                                                        "n_points",
                                                                        "n_solutions",
                                                                        "booststrap_replicates", 
                                                                        "lower_bound",
                                                                        "upper_bound",
                                                                        "ci_bootstrap_type")))
                               })

tpoly_summaries <- dplyr::bind_rows(tpoly_summaries_list)

##### Add the areal estimates --------------------------------------------------
tpoly_summaries <- dplyr::filter(.data = tpoly_summaries,
                                 !(indicator %in% categorical_indicators)) |>
  dplyr::left_join(x = _,
                   y = sf::st_drop_geometry(x = output_polygons) |>
                     dplyr::select(.data = _,
                                   -polygon_type),
                   relationship = "many-to-one",
                   by = c("frame" = "uid")) |>
  dplyr::mutate(.data = _,
                # Because these were in percentages instead of proportions
                mean_estimate = dplyr::case_when(indicator == "pinyon_cover" ~ mean_estimate / 100,
                                                 .default = mean_estimate),
                acres_estimate = round(mean_estimate * area_acres,
                                       digits = 1),
                acres_lower_bound_acres = round(lower_bound * area_acres,
                                                digits = 1),
                acres_upper_bound_acres = round(upper_bound * area_acres,
                                                digits = 1)) |>
  dplyr::rename(.data = _,
                "proportion_estimate" = "mean_estimate",
                "proportion_lower_bound" = "lower_bound",
                "proportion_upper_bound" = "upper_bound")

#### WRITING ###################################################################
##### Make the outputs ---------------------------------------------------------
###### Widen results -----------------------------------------------------------
# This is super reduced, but it's just to throw on the polygons. The full table
# is still intended to be authoritative.
output_results_wide <- dplyr::select(.data = tpoly_summaries,
                                     tidyselect::all_of(x = c("frame",
                                                              "indicator",
                                                              "proportion_estimate",
                                                              "acres_estimate"))) |>
  tidyr::pivot_longer(data = _,
                      cols = tidyselect::matches(match = "_estimate"),
                      names_to = "type") |>
  dplyr::mutate(.data = _,
                indicator = stringr::str_replace(string = indicator,
                                                 pattern = "%",
                                                 replacement = " percent"),
                indicator = dplyr::case_when(!stringr::str_detect(string = indicator,
                                                                  pattern = "cover") ~ paste0("pinyon_cover_",
                                                                                              indicator),
                                             .default = indicator),
                type = stringr::str_remove(string = type,
                                           pattern = "_estimate"),
                value = round(x = value,
                              digits = 2)) |>
  tidyr::pivot_wider(data = _,
                     names_from = tidyselect::all_of(x = c("indicator",
                                                           "type")),
                     names_glue = "{indicator} ({type})",
                     values_from = value) |>
  dplyr::rename(.data = _,
                "pinyon cover (percent)" = "pinyon_cover (proportion)") |>
  dplyr::mutate(.data = _,
                `pinyon cover (percent)` = `pinyon cover (percent)` * 100)

# output_continuous_results_wide <- dplyr::select(.data = tpoly_summaries_continuous,
#                                                 frame,
#                                                 indicator,
#                                                 mean_estimate) |>
#   tidyr::pivot_wider(data = _,
#                      names_from = indicator,
#                      values_from = mean_estimate,
#                      values_fn = ~ round(x = .x,
#                                          digits = 1)) |>
# dplyr::rename(.data = _,
#               "PJ cover (%)" = "pj_cover")

# output_results_wide <- dplyr::full_join(x = output_categorical_results_wide,
#                                         y = output_continuous_results_wide,
#                                         by = c("frame"))

###### Join the results to the polygon output ----------------------------------
output_polygons <- dplyr::left_join(x = output_polygons,
                                    y = output_results_wide,
                                    by = c("uid" = "frame")) |>
  dplyr::rename_with(.data = _,
                     .fn = ~ stringr::str_replace(string = .x,
                                                  pattern = "-",
                                                  replacement = " to ") |>
                       stringr::str_remove_all(string = _,
                                               pattern = "[\\(\\)]") |>
                       stringr::str_replace_all(string = _,
                                                pattern = " ",
                                                replacement = "_") |>
                       tolower())
# sf::st_geometry(output_polygons) <- "SHAPE"
# sf::st_geometry_type(output_polygons) <- "SHAPE"

###### Join the information to the points --------------------------------------
output_points <- dplyr::left_join(x = headers,
                                  y = pinyon_cover,
                                  by = "PrimaryKey") |>
  sf::st_intersection(x = _,
                      y = polygons_extent) |>
  dplyr::select(.data = _,
                -tidyselect::all_of(x = c("frame_id")))


##### Write! -------------------------------------------------------------------
sf::st_write(obj = output_polygons,
             dsn = file.path(output_path,
                             output_gdb_filename),
             layer = "ecoregion_polygons",
             append = FALSE)

sf::st_write(obj = output_points,
             dsn = file.path(output_path,
                             output_gdb_filename),
             layer = "ecoregion_points",
             append = FALSE)

sf::st_write(obj = analysis_results,
             dsn = file.path(output_path,
                             output_gdb_filename),
             layer = "raw_analysis_results",
             append = FALSE)

sf::st_write(obj = tpoly_summaries,
             dsn = file.path(output_path,
                             output_gdb_filename),
             layer = "analysis_results",
             append = FALSE)

write.csv(x = analysis_results,
          file = paste0(output_path, "/", analysis_name, "_results_", file_date, ".csv"),
          row.names = FALSE)

write.csv(x = tpoly_summaries,
          file = paste0(output_path, "/", analysis_name, "_tpoly_summaries_", file_date, ".csv"),
          row.names = FALSE)
