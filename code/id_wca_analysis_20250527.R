#### OVERVIEW ##################################################################
# This script is designed to make weighted estimates across one or more frames
# (i.e., inference areas) using AIM and LMF data.

# Required packages, which I think is comprehensive:
packages <- c("geosphere",
              "sf",
              # Hopefully not this one, but some functions that I think are
              # deprecated reference it despite it being supplanted by sf.
              "sp",
              "tidyverse",
              "spatstat",
              "boot",
              "BAMMtools",
              "raster",
              "stars",
              "httr")

# The inputs you'll need are:
# A geodatabase with the following feature classes
#  - Polygons describing all frames you want analyses for. These can be
#    multipart polygons, but the attribute *must* have a variable called
#    "frame_id" which contains the unique identifier for each frame. Don't use
#    special characters or spaces in those IDs, e.g., "Deer Valley/Fish Creek"
#    should be changed to something like "deer-valley-fish-creek".
#  - Polygons describing the LMF segments. LMF points are drawn using a two-stage
#    approach, so we need to weight points taking that into account using the
#    segments that represent the first stage of point selection. LMF points are
#    named such that the ID of the segment the point belongs to is everything
#    except the final digit of the name which denotes the point identity within
#    that segment.
#
# AIM and LMF points
#  - Right now, these are written to be pulled from their own geodatabases and
#    assume that the feature classes in those geodatabases are called
#    "BLM_Natl_AIM_TerrADat_Hub" and "BLM_Natl_AIM_LMF_Hub". Feel free to modify
#    that if you need to.

# Other things you'll need:
#  - A copy of functions.R containing the support functions for using Thiessen
#    polygons and density partitions.

# Parameters for the analyses can be set in the SETUP section below.

#### SETUP #####################################################################
##### Parameters ---------------------------------------------------------------
###### Inputs ------------------------------------------------------------------
# The full filepath to the geodatabase containing the frames and LMF segments.
data_path <- "C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/analysis/idaho_wca/data/Owyhee.gdb"

# The name of the feature class containing the frames to use. The feature class
# needs to have a variable called frame_id containing the unique identifiers for
# the frames in the feature class.
frame_feature_class_name <- "huc10s_blm"

# The name of the feature class containing the LMF segments. Each segment must
# have a unique ID in the POLY_ID variable of the attribute table.
lmf_segment_layer_name <- "lmf_segments_dissolve"

# The full filepaths to the geodatabases containing the AIM and LMF data and the
# names of the relevant feature classes in those geodatabases.
# aim_points_path <- "C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/analysis/idaho_wca/data/BLM_Natl_AIM_TerrADat_Public.gdb"
# aim_points_feature_class_name <- "TerrADat"
# lmf_points_path <- "C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/analysis/idaho_wca/data/BLM_Natl_AIM_LMF_Public.gdb"
# lmf_points_feature_class_name <- "LMF"
aim_points_path <- data_path
aim_points_feature_class_name <- "Owyhee_TerrADat"
lmf_points_path <- data_path
lmf_points_feature_class_name <- "Owyhee_LMF"

# The full filepath to the functions.R that contains all the required support
# functions.
functions_filepath <- "C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/code/functions.R"

###### Outputs -----------------------------------------------------------------
# The full path to the folder where outputs will be saved.
output_path <- "C:/Users/Nelson/Documents/Projects/thiessen_polygon_estimates/analysis/idaho_wca/output"

# The base name for the analyses for the frames being analyzed. This will be
# used in the output filenames.
analysis_name <- "id_wca_watersheds"

# The date to append to the output filenames.
file_date <- "20250527"

###### Configuration -----------------------------------------------------------
## General ---
# For most analyses, a confidence of 80% will be sufficient for management
# purposes.
conf <- 80

# The minimum number of points required in a frame in order to move ahead with
# making an estimate for that frame. This probably shouldn't be set below 10
# due to the level of uncertainty with small sample sizes.
min_points <- 10


## Indicators ---
# These are the names of the variables containing the indicators in the AIM and
# LMF point attribute tables. They may not all be in the feature classes in the
# geodatabases, but any that aren't need to be calculated below in the section
# "Calculating bonus indicators".
# It is important that continuous and categorical indicators are assigned
# correctly because the analysis step handles them differently.
continuous_indicators <- c("AH_AnnGrassCover",
                           "AH_PerenGrassCover",
                           "AH_PerenForbGrassCover",
                           "AH_NonNoxAnnForbGrassCover",
                           "AH_NoxAnnGrassCover",
                           "AH_SagebrushCover",
                           "AH_NonNoxTreeCover",
                           "AH_NoxTreeCover",
                           "juniper_cover",
                           "cheatgrass_cover")

categorical_indicators <- c("cheatgrass_present",
                            "cheatgrass_cover_over10",
                            "NoxAnnGrass_cover_over5",
                            "NoxAnnGrass_cover_over10")

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

# How many Thiessen polygons to draw within a frame.
# There will be one estimate per specified Thiessen polygon count for each
# frame, e.g., c(2, 3) would produce estimates with 2 Thiessen polygons and also
# with 3 Thiessen polygons.
tpoly_counts <- c(2, 3)

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
density_partition_counts <- 2:3

# This is the same as tpoly_counts above except that this is the number of
# Thiessen polygons that will be drawn within the density partitions.
# tpoly_counts controls only estimates made without using density partitioning.
dp_tpoly_counts <- c(2, 3)

# This is the sampling distance in meters for calculating density.
# The smaller this number, the finer-grained the density map but that can cause
# serious memory allocation issues when you have larger areas.
# Don't worry about changing this unless you're analyzing very large areas and
# are running into time or memory constraints.
density_sample_spacing <- 100

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

#### INITIALIZATION ############################################################
# Everything beyond here should be automated and does not require user changes.

##### Getting functions set ----------------------------------------------------
# This is necessary because simply calling from the namespace doesn't access the
# support functions in the package (or maybe from support packages it wants to
# attach?)
library(spatstat)

# Bring in the functions from the script specified above.
source(functions_filepath)

##### Reading data, making objects, and storing values -------------------------
# Because things like buffering won't work with degrees as units, we'll use this
# to transform the polygon inputs sometimes.
aea_proj <- "+proj=aea +lat_1=29.5 +lat_2=42.5"

# Just combining the indicators into a single vector for later.
analysis_indicators <- c(continuous_indicators,
                         categorical_indicators)

# Read in the frames from the geodatabase.
frames <- sf::st_read(dsn = data_path,
                      layer = frame_feature_class_name) |>
  sf::st_transform(x = _,
                   crs = aea_proj) |>
  sf::st_make_valid(x = _)

# Make a list of the frames individually so that we can iterate over them below.
frames_list <- lapply(X = setNames(object = frames$frame_id,
                                   nm = frames$frame_id),
                      frames = frames,
                      FUN = function(X, frames){
                        dplyr::filter(.data = frames,
                                      frame_id == X)
                      })

# Read in the AIM points and keep only the indicators we need.
aim_points <- sf::st_read(dsn = aim_points_path,
                          layer = aim_points_feature_class_name) |>
  sf::st_transform(x = _,
                   crs = aea_proj) |>
  dplyr::select(.data = _,
                PrimaryKey,
                tidyselect::any_of(analysis_indicators)) |>
  dplyr::mutate(.data = _,
                source = "AIM") |>
  sf::st_intersection(x = _,
                      y = dplyr::select(.data = frames,
                                        frame_id))

# Read in the LMF points and keep only the indicators we need.
lmf_points <- sf::st_read(dsn = lmf_points_path,
                          layer = lmf_points_feature_class_name) |>
  sf::st_transform(x = _,
                   crs = aea_proj) |>
  dplyr::select(.data = _,
                PrimaryKey,
                tidyselect::any_of(analysis_indicators)) |>
  dplyr::mutate(.data = _,
                source = "LMF") |>
  sf::st_intersection(x = _,
                      y = dplyr::select(.data = frames,
                                        frame_id))

# Combine the AIM and LMF points.
all_points <- dplyr::bind_rows(aim_points,
                               lmf_points)


###### Calculating additional indicators ---------------------------------------
# This particular example also uses some indicators that aren't in TerrADat, so
# this bit calculates them.

# It grabs the LPI data from the LDC for the points in the frames, calculates
# the indicators, and adds them to all_points.

pks <- sf::st_drop_geometry(x = all_points) |>
  dplyr::pull(.data = _,
              var = PrimaryKey) |>
  unique()

lpi_detail <- sf::st_read(dsn = aim_points_path,
                          layer = "tblLPIDetail") |>
  dplyr::filter(.data = _,
                PrimaryKey %in% pks)
lpi_headers <- sf::st_read(dsn = aim_points_path,
                           layer = "tblLPIHeader") |>
  dplyr::filter(.data = _,
                PrimaryKey %in% pks)

aim_lpi_data <- terradactyl::gather_lpi_terradat(tblLPIDetail = lpi_detail,
                                                 tblLPIHeader = lpi_headers)


pintercept <- sf::st_read(dsn = lmf_points_path,
                          layer = "PINTERCEPT") |>
  dplyr::filter(.data = _,
                PrimaryKey %in% pks)

lmf_lpi_data <- terradactyl::gather_lpi_lmf(PINTERCEPT = pintercept)


lpi_data <- dplyr::bind_rows(aim_lpi_data,
                             lmf_lpi_data)

missing_from_raw_data_pks <- setdiff(x = pks,
                                     y = unique(lpi_data$PrimaryKey))


tblnationalplants <- read.csv(file = "C:/Users/Nelson/Documents/Projects/terradactyl/testing/tblNationalPlants.csv")

juniper_codes <- dplyr::filter(.data = tblnationalplants,
                               stringr::str_detect(string = name,
                                                   pattern = "^Juniperus")) |>
  dplyr::pull(.data = _,
              namCode) |>
  unique()

lpi_data <- dplyr::mutate(.data = lpi_data,
                          juniper = dplyr::case_when(code %in% juniper_codes ~ "juniper",
                                                     .default = NA),
                          cheatgrass = dplyr::case_when(code %in% c("BRTE") ~ "cheatgrass",
                                                        .default = NA))

all_points <- dplyr::inner_join(x = all_points,
                                y = terradactyl::pct_cover(lpi_tall = lpi_data,
                                                           tall = FALSE,
                                                           hit = "any",
                                                           by_line = FALSE,
                                                           cheatgrass) |>
                                  dplyr::rename(.data = _,
                                                cheatgrass_cover = tidyselect::matches(match = "cheatgrass")) |>
                                  dplyr::mutate(.data = _,
                                                cheatgrass_present = cheatgrass_cover > 0,
                                                cheatgrass_cover_over10 = cheatgrass_cover >= 10),
                                relationship = "one-to-one",
                                by = "PrimaryKey") |>
  dplyr::left_join(x = _,
                   y = terradactyl::pct_cover(lpi_tall = lpi_data,
                                              tall = FALSE,
                                              hit = "any",
                                              by_line = FALSE,
                                              juniper) |>
                     dplyr::rename(.data = _,
                                   juniper_cover = JUNIPER),
                   relationship = "one-to-one",
                   by = "PrimaryKey") |>
  dplyr::mutate(.data = _,
                NoxAnnGrass_cover_over5 = AH_NoxAnnGrassCover > 5,
                NoxAnnGrass_cover_over10 = AH_NoxAnnGrassCover > 10) |>
  dplyr::select(.data = _,
                PrimaryKey,
                source,
                tidyselect::all_of(x = analysis_indicators))

#### PROCESSING ################################################################
start_time <- Sys.time()

for (current_analysis_frame_id in names(frames_list)) {
  message(toupper(x = paste("Working on frame ID",
                            current_analysis_frame_id,
                            paste0("(", which(names(frames_list) == current_analysis_frame_id), " of ", length(frames_list), ")"))))
  
  ##### Reading and munging ----------------------------------------------------
  points <- sf::st_intersection(x = all_points,
                                y = frames_list[[current_analysis_frame_id]])
  
  if (nrow(points) >= min_points) {
    ##### Reading in polygons --------------------------------------------------
    
    # Read in the LMF segments
    segments <- sf::st_read(dsn = data_path,
                            layer = lmf_segment_layer_name) |>
      dplyr::select(.data = _,
                    segment_id = tidyselect::matches(match = "id$")) |>
      sf::st_transform(x = _,
                       crs = aea_proj) |>
      sf::st_make_valid(x = _) |>
      sf::st_buffer(x = _,
                    dist = 0)
    
    segments$area_m2 <- as.numeric(sf::st_area(x = segments))
    
    
    ##### Attributing points ---------------------------------------------------
    # We need to use the LMF segments to find a "relative weight"
    # For points outside a segment, that relative weight is going to be a 1
    # For points inside a segment, it's going to be the normalized area of the
    # segment (i.e., the largest segment involved will have a value of 1) divided by
    # the number of points evaluated within it.
    points_segments_lut <- sf::st_intersection(x = dplyr::select(.data = points,
                                                                 PrimaryKey,
                                                                 source),
                                               y = dplyr::select(.data = segments,
                                                                 uid_segment = segment_id))
    
    if (nrow(points_segments_lut) > 0) {
      segments <- dplyr::filter(.data = segments,
                                segment_id %in% points_segments_lut$uid_segment)
      
      segments <- dplyr::mutate(.data = segments,
                                # segment_id = dplyr::row_number(),
                                relative_segment_area = as.numeric(area_m2 / max(segments$area_m2)))
      
      points_relative_weights_lut <- apply(X = sf::st_drop_geometry(segments),
                                           MARGIN = 1,
                                           points_segments_lut = points_segments_lut,
                                           FUN = function(X, points_segments_lut){
                                             current_segment_id <- X[["segment_id"]]
                                             
                                             # No clue why, but this is turning into a character value???
                                             current_relative_segment_area <- as.numeric(X[["relative_segment_area"]])
                                             
                                             current_points <- dplyr::filter(.data = points_segments_lut,
                                                                             uid_segment == current_segment_id)
                                             aim_count <- sum(current_points$source == "AIM")
                                             lmf_sampled_count <- sum(current_points$source == "LMF")
                                             # To get the number of evaluated but not sampled points:
                                             # The LMF plot keys end in a digit that represents the intended sampling order within a segment
                                             # 1 and 2 are considered base points and were intended to be sampled
                                             # If a sampled LMF plot's plot key ends in 3, that means that one or both of the base points
                                             # were evaluated and rejected rather than sampled, which brings the evaluated LMF plot count
                                             # to three for the segment.
                                             # This just asks if the third point was used
                                             lmf_suffixes <- dplyr::filter(.data = current_points,
                                                                           source == "LMF")[["PrimaryKey"]] |>
                                               stringr::str_extract(string = _,
                                                                    pattern = "(?<=\\D)[1-3]") |>
                                               unlist() |>
                                               as.numeric()
                                             
                                             lmf_oversample_used <- 3 %in% lmf_suffixes
                                             
                                             # Likewise, if only one LMF plot was sampled in a segment, that means the other two were
                                             # evaluated and rejected rather than sampled, also bringing the total to three.
                                             # So if there was only one sampled or if the oversample was used, there were three evaluated
                                             if (lmf_sampled_count == 1 | lmf_oversample_used) {
                                               lmf_count <- 3
                                             } else {
                                               # This will fire only if there sampled count was 2, but better to be safe here
                                               lmf_count <- lmf_sampled_count
                                             }
                                             
                                             dplyr::mutate(.data = sf::st_drop_geometry(current_points),
                                                           relative_weight = current_relative_segment_area / (aim_count + lmf_count)) |>
                                               dplyr::select(.data = _,
                                                             PrimaryKey,
                                                             uid_segment,
                                                             relative_weight)
                                           }) |>
        dplyr::bind_rows()
      
      points <- dplyr::left_join(x = points,
                                 y = sf::st_drop_geometry(points_relative_weights_lut),
                                 by = "PrimaryKey",
                                 relationship = "one-to-one") |>
        dplyr::mutate(.data = _,
                      relative_weight = tidyr::replace_na(data = relative_weight,
                                                          replace = 1))
    } else {
      message("Insufficient points available for a proper weighted analysis.")
      points <- dplyr::mutate(.data = points,
                              relative_weight = 1)
    }
    
    
    
    ##### Filtering and making polygons ----------------------------------------
    # We want a maximum inference area extent. This can only be used with tpolys
    # and/or density partitions if there are any poststrata that didn't have points.
    analysis_frame <- frames_list[[current_analysis_frame_id]]
    
    
    #### WEIGHTING #################################################################
    ##### Equal ####################################################################
    # No tricks here. We just need this to define equal weights the same way we
    # do the derived weights.
    points_attributed_equal_list <- list()
    points_attributed_equal_list[[1]] <- dplyr::select(.data = sf::st_drop_geometry(points),
                                                       PrimaryKey) |>
      dplyr::mutate(.data = _,
                    analysis_frame_id = current_analysis_frame_id,
                    weightcat_id = "None",
                    weight = 1,
                    weighting = "Equal weighting",
                    solution = 1,
                    partitions = 0,
                    tpolys = 0)
    
    ##### Thiessen polygons ########################################################
    ###### Finding Thiessen polygons ###############################################
    # This is a loop instead of a lapply() so that we get to keep solutions even if
    # one of them fails
    tpoly_list <- list()
    
    # We'll go through this one level of partitioning at a time
    # Get ready for a gnarly nested list situation because it's the most obvious solution
    # to help keep things organized
    for (tpoly_count in tpoly_counts) {
      message(paste0("Beginning to find tpoly solutions with ", tpoly_count, " polygons."))
      current_count_tpoly_list <- list()
      
      # for (frame_name in c("maximum_inference_area", "poststratified_inference_area")) {
      for (tpoly_seed_index in 1:length(tpoly_seeds)) {
        message(paste0("Starting with seed ", tpoly_seed_index, " of ", length(tpoly_seeds)))
        current_tpoly_seed <- tpoly_seeds[tpoly_seed_index]
        
        current_seed_tpoly_list <- list()
        
        if (nrow(points) < min_points_per_tpoly * tpoly_count) {
          message("Inadequate number of points to meet minimum points-per-polygon requirement. Skipping.")
          current_tpolys <- NULL
        } else {
          current_tpolys <- thiessen_polygons_gen_random(frame = analysis_frame,
                                                         n_polygons = tpoly_count,
                                                         points = points,
                                                         points_min = min_points_per_tpoly,
                                                         envelope = sf::st_union(analysis_frame),
                                                         seed_number = current_tpoly_seed,
                                                         seed_increment = 1000,
                                                         iteration_limit = 500,
                                                         use_albers = TRUE,
                                                         verbose = TRUE)
          
          if (is.null(current_tpolys)) {
            message("Unable to find a solution.")
          } else {
            current_tpolys <- dplyr::mutate(.data = current_tpolys,
                                            analysis_frame_id = current_analysis_frame_id,
                                            partition_id = "None",
                                            tpoly_id = paste0(tpoly_id),
                                            partition_count = 0,
                                            tpoly_count = tpoly_count,
                                            seed_id = tpoly_seed_index)
          }
        }
        message(paste0("Completed attempt to find solution ", tpoly_seed_index, " of ", length(tpoly_seeds), " with ",
                       tpoly_count, " tpolys."))
        
        
        
        tpoly_list[[paste0(current_analysis_frame_id, "_frame-",
                           tpoly_seed_index, "_solution-",
                           "0_partitions-",
                           tpoly_count, "_tpolys")]] <- current_tpolys
      }
      # }
    }
    
    expected_tpoly_list_length <- length(tpoly_seeds) * length(tpoly_counts)
    if (expected_tpoly_list_length > length(tpoly_list)) {
      warning(paste0("There were only ", length(tpoly_list), " Thiessen polygon solutions found but you expected ",
                     expected_tpoly_list_length, ". This is almost certainly due to low point counts in one or more density partitions making Thiessen polygons impossible with the required minimum point count or simply so unlikely that the iteration limit was hit."))
    }
    
    ###### Calculating weights -------------------------------------------------
    points_attributed_tpoly_list <- lapply(X = names(tpoly_list),
                                           tpoly_list = tpoly_list,
                                           points_relevant = dplyr::select(points,
                                                                           PrimaryKey),
                                           FUN = function(X, tpoly_list, points_relevant){
                                             message(paste0("Attributing from tpoly set ", which(names(tpoly_list) == X), " of ", length(tpoly_list)))
                                             parameter_vector <- stringr::str_split_1(string = X,
                                                                                      pattern = "-")
                                             names(parameter_vector) <- stringr::str_extract(string = parameter_vector,
                                                                                             pattern = "[a-z]+$")
                                             parameter_vector <- gsub(x = parameter_vector,
                                                                      pattern = "_[a-z]+$",
                                                                      replacement = "")
                                             
                                             current_tpolys <- tpoly_list[[X]]
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
                                             current_points_attributed[["weighting"]] <- "Thiessen polygons"
                                             for (parameter in names(parameter_vector)) {
                                               current_points_attributed[[parameter]] <- as.numeric(parameter_vector[parameter])
                                             }
                                             sf::st_drop_geometry(current_points_attributed)
                                           })
    
    
    ##### Density partitions #######################################################
    ###### Finding density partitions ##############################################
    # We're going to find solutions with varying numbers of tpolys, multiple solutions
    # for each count
    # This is going to be within density partitions, so we'll solve for those first
    
    # This next stretch is basically a vivisected density_polygon_gen_clustered()
    
    density_partitions_list <- lapply(X = density_partition_counts,
                                      current_frame_id = current_analysis_frame_id,
                                      frame = frames_list[[current_analysis_frame_id]],
                                      points = points,
                                      buffer_distance = buffer_distance,
                                      FUN = function(X, current_frame_id, frame, points, buffer_distance){
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
                                        if (length(unique(density_sf$partition_id)) != intended_partition_count) {
                                          NULL
                                        } else {
                                          density_sf
                                        }
                                      })
    
    names(density_partitions_list) <- paste0(current_analysis_frame_id, "_",
                                             density_partition_counts)
    density_partitions_list <- density_partitions_list[!sapply(X = density_partitions_list,
                                                               FUN = is.null)]
    
    ###### Calculating weights -------------------------------------------------
    # Let's get the weights calculated and added for density partitioning only
    points_attributed_density_list <- lapply(X = density_partitions_list,
                                             points_relevant = dplyr::select(points,
                                                                             PrimaryKey),
                                             FUN = function(X, points_relevant){
                                               current_weighting_polygons_list <- lapply(X = unique(X$partition_id),
                                                                                         n_partitions = length(unique(X$partition_id)),
                                                                                         polygons = X,
                                                                                         FUN = function(X, n_partitions, polygons){
                                                                                           current_polygons <- sf::st_as_sf(sf::st_union(x = dplyr::filter(polygons,
                                                                                                                                                           partition_id == X)))
                                                                                           current_polygons$partition_id <- paste0(X, "_", n_partitions)
                                                                                           current_polygons$area <- as.numeric(sf::st_area(current_polygons))
                                                                                           current_polygons
                                                                                         })
                                               
                                               current_weighting_polygons <- dplyr::select(.data = dplyr::bind_rows(current_weighting_polygons_list),
                                                                                           wgtcat_id = partition_id,
                                                                                           area)
                                               
                                               current_points_attributed <- sf::st_join(x = points_relevant,
                                                                                        y = current_weighting_polygons[, "wgtcat_id"],
                                                                                        # Note that there may be lost points
                                                                                        # This is due to points that fall just outside
                                                                                        # the boundaries due to the resolution of the
                                                                                        # density raster used. We could compensate
                                                                                        # by assigning to the nearest polygon,
                                                                                        # but we'll skip that for now
                                                                                        left = FALSE)
                                               current_weight_summary <- dplyr::summarize(.data = dplyr::group_by(.data = sf::st_drop_geometry(current_points_attributed),
                                                                                                                  wgtcat_id),
                                                                                          n_points = dplyr::n())
                                               current_weight_summary <- dplyr::left_join(x = current_weight_summary,
                                                                                          y = sf::st_drop_geometry(current_weighting_polygons),
                                                                                          by = "wgtcat_id")
                                               current_weight_summary <- dplyr::mutate(.data = current_weight_summary,
                                                                                       weight = area / n_points)
                                               current_points_attributed <- dplyr::left_join(x = current_points_attributed,
                                                                                             y = dplyr::select(current_weight_summary,
                                                                                                               wgtcat_id,
                                                                                                               weight),
                                                                                             by = "wgtcat_id")
                                               current_points_attributed[["analysis_frame_id"]] <- X[["analysis_frame_id"]][1]
                                               current_points_attributed[["weighting"]] <- "Density partitions"
                                               current_points_attributed[["solution"]] <- 1
                                               current_points_attributed[["partitions"]] <- nrow(current_weighting_polygons)
                                               current_points_attributed[["tpolys"]] <- 0
                                               sf::st_drop_geometry(current_points_attributed)
                                             })
    
    ##### Within-partition Thiessen polygons -----------------------------------
    ###### Finding Thiessen polygons -------------------------------------------
    # This is a loop instead of a lapply() so that we get to keep solutions even if
    # one of them fails
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
              
              points_available <- sf::st_join(x = points,
                                              y = current_density_partition) |>
                dplyr::filter(.data = _,
                              !is.na(partition_id))
              if (nrow(points_available) < min_points_per_tpoly * tpoly_count) {
                message("Inadequate number of points to meet minimum points-per-polygon requirement. Skipping.")
                current_tpolys <- NULL
              } else {
                current_tpolys <- thiessen_polygons_gen_random(frame = current_density_partition,
                                                               n_polygons = tpoly_count,
                                                               points = points,
                                                               points_min = min_points_per_tpoly,
                                                               envelope = sf::st_union(current_density_partition),
                                                               seed_number = current_tpoly_seed,
                                                               seed_increment = 1000,
                                                               iteration_limit = 500,
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
              current_tpolys <- NULL
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
    
    expected_tpoly_list_length <- length(density_partitions_list) * length(tpoly_counts) * length(tpoly_seeds)
    if (expected_tpoly_list_length > length(dp_tpoly_list)) {
      warning(paste0("There were only ", length(dp_tpoly_list), " Thiessen polygon solutions found but you expected ",
                     expected_tpoly_list_length, ". This is almost certainly due to low point counts in one or more density partitions making Thiessen polygons impossible with the required minimum point count or simply so unlikely that the iteration limit was hit."))
    }
    
    ###### Calculating weights -------------------------------------------------
    points_attributed_dp_tpoly_list <- lapply(X = names(dp_tpoly_list),
                                              dp_tpoly_list = dp_tpoly_list,
                                              points_relevant = dplyr::select(points,
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
    
    
    
    #### ANALYSIS ##############################################################
    ##### Analyze for each set of weights --------------------------------------
    weighting_list <- c(points_attributed_equal_list,
                        points_attributed_tpoly_list,
                        points_attributed_density_list,
                        points_attributed_dp_tpoly_list)
    
    analysis_list <- lapply(X = weighting_list,
                            points = sf::st_drop_geometry(points) |>
                              tidyr::pivot_longer(data = _,
                                                  cols = tidyselect::all_of(analysis_indicators),
                                                  names_to = "indicator",
                                                  values_to = "value") |>
                              dplyr::select(.data = _,
                                            PrimaryKey,
                                            indicator,
                                            value) |>
                              dplyr::mutate(.data = _,
                                            value = tidyr::replace_na(data = value,
                                                                      replace = 0)),
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
                                                                        analysis <- data.frame(n = nrow(points_weighted),
                                                                                               alpha = 1 - conf / 100,
                                                                                               estimate = points_weighted$value[1],
                                                                                               sd = 0,
                                                                                               variance = 0,
                                                                                               lower_bound = NA,
                                                                                               upper_bound = NA)
                                                                      } else {
                                                                        # Make sure we calculate using the correct method
                                                                        if (current_indicator %in% continuous_indicators) {
                                                                          analysis <- analyze_con(data = X,
                                                                                                  weights = current_weights,
                                                                                                  id_var = "PrimaryKey",
                                                                                                  value_var = "value",
                                                                                                  wgt_var = "weight",
                                                                                                  conf = conf,
                                                                                                  verbose = TRUE) |>
                                                                            dplyr::rename(.data = _,
                                                                                          estimate = mean)
                                                                        } else if (current_indicator %in% categorical_indicators) {
                                                                          analysis <- analyze_cat(data = dplyr::mutate(.data = X,
                                                                                                                       value = as.character(as.logical(value))),
                                                                                                  weights = current_weights,
                                                                                                  id_var = "PrimaryKey",
                                                                                                  cat_var = "value",
                                                                                                  wgt_var = "weight",
                                                                                                  definitions = c("TRUE",
                                                                                                                  "FALSE"),
                                                                                                  conf = 80,
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
                                                                            dplyr::filter(.data = _,
                                                                                          category == "TRUE") |>
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
    
    write.csv(x = analysis_results,
              file = paste0(output_path, "/", analysis_name, "_", current_analysis_frame_id, "_results_", file_date, ".csv"),
              row.names = FALSE)
    
    ##### Combine tpoly results from equivalent solutions ----------------------
    # Let's get those tpoly results combined to get a single estimate with a
    # confidence interval
    tpoly_analyses <- dplyr::filter(.data = analysis_results,
                                    # n_tpolys > 0,
                                    grepl(weighting,
                                          pattern = "polygons"))
    
    tpoly_summaries_list <- lapply(X = split(x = tpoly_analyses,
                                             f = list(tpoly_analyses$indicator,
                                                      tpoly_analyses$frame,
                                                      tpoly_analyses$n_tpolys,
                                                      tpoly_analyses$n_partitions),
                                             drop = TRUE),
                                   bootstrap_types = bootstrap_types,
                                   bootstrap_count = bootstrap_count,
                                   FUN = function(X, bootstrap_types, bootstrap_count){
                                     # message(X$indicator[1])
                                     
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
                                     
                                     output
                                   })
    
    tpoly_summaries <- dplyr::bind_rows(tpoly_summaries_list)
    
    write.csv(x = tpoly_summaries,
              file = paste0(output_path, "/", analysis_name, "_", current_analysis_frame_id, "_tpoly_summaries_", file_date, ".csv"),
              row.names = FALSE)
    
    
  } else {
    message("There were no available data for this frame.")
  }
}

completion_time <- Sys.time()

beepr::beep(5)

completion_time - start_time

list.files(path = output_path,
           full.names = TRUE,
           pattern = paste0(analysis_name, ".+tpoly.+", file_date)) |>
  lapply(X = _,
         FUN = function(X){
           tryCatch(read.csv(file = X,
                             stringsAsFactors = FALSE),
                    error = function(e){
                      NULL
                    })
         }) |>
  dplyr::bind_rows() |>
  write.csv(x = _,
            file = file.path(output_path,
                             paste0(analysis_name,
                                    "_tpoly_summaries_",
                                    file_date,
                                    ".csv")))
