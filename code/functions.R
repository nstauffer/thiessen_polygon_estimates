# Required packages (these are called explicitly by namespace and do not need to be attached)
# library(geosphere)
# library(sf)
# library(spsurvey)

#' Create Thiessen/Voronoi polygons from a set of points and bounding polygons
#' @description Generate Thiessen/Voronoi polygons for a set of points and clip the results using a set of polygons
#' @param centroids An sf points object. These points are used as centroids for the Thiessen/Voronoi polygons.
#' @param frame An sf polygon or multipolygon object. This is the clipping boundary which will be applied to the otherise "infinite" Thiessen/Voronoi polygons.
#' @param envelope An sfc polygon object. This will be the outer envelope for the Thiessen polygons before they're clipped to \code{frame}. This will only be applied if it's larger than the default envelope in \code{sf::st_voronoi()}. If \code{NULL} then the default envelope will be used. Defaults to \code{NULL}.
#' @param use_albers Logical. If \code{TRUE} then \code{centroids} and \code{frame} will be reprojected into Albers Equal Area (AEA) and the output will be in AEA. If \code{FALSE} then \code{frame} will be reprojected to match the coordinate reference ssytem (CRS) of \code{centroids} and the output will be in that CRS. CRSs using decimal degrees will throw errors or warnings. Defaults to \code{TRUE}.
#' @return An sf object composed of polygon or multipolygon geometry
thiessen_polygons_gen_fixed <- function(centroids,
                                        frame,
                                        envelope = NULL,
                                        use_albers = TRUE) {
  # Define Alber's Equal Area CRS
  aea_proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  
  # Sanitization
  if (!("sf" %in% class(centroids))) {
    stop("`centroids` must be an sf points object")
  } else if (!all(sf::st_geometry_type(centroids) %in% c("POINT"))) {
    stop("`centroids` must be an sf points object")
  }
  if (!("sf" %in% class(frame))) {
    stop("`frame` must be an sf polygon object")
  } else if (!all(sf::st_geometry_type(frame) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("`frame` must be an sf polygon object")
  }
  
  if (!is.null(envelope)) {
    if (!("sfc" %in% class(envelope))) {
      stop("`envelope` must be an sfc polygon object")
    } else if (!all(sf::st_geometry_type(envelope) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("`envelope` must be an sfc polygon object")
    }
  }
  
  # Remove any Z dimension
  # It screws with the process and is irrelevant
  centroids <- sf::st_zm(centroids,
                         drop = TRUE)
  frame <- sf::st_zm(frame,
                     drop = TRUE)
  
  # Reproject as necessary
  if (use_albers) {
    centroids <- sf::st_transform(x = centroids,
                                  crs = aea_proj)
    frame <- sf::st_transform(x = frame,
                              crs = aea_proj)
    envelope <- sf::st_transform(x = envelope,
                                 crs = aea_proj)
  } else {
    # This just forces the polygons into the same projection as the centroids
    centroids_crs <- sf::st_crs(centroids)
    frame <- sf::st_transform(frame,
                              crs = centroids_crs)
    envelope <- sf::st_transform(envelope,
                                 crs = centroids_crs)
  }
  
  # Draw Thiessen polygons
  # Here's where it gets weird
  # The points need to be a multipoint object, apparently
  points_multipoint <- sf::st_combine(sf::st_geometry(centroids))
  
  # Generate the Thiessen polygons
  # If there's a provided envelope, we attempt to use it
  # This appears to be a list??? It's fine, I promise. We'll convert it in a bit
  if (is.null(envelope)) {
    thiessen_polygons_raw <- sf::st_voronoi(x = points_multipoint)
  } else {
    thiessen_polygons_raw <- sf::st_voronoi(x = points_multipoint,
                                            envelope = envelope)
  }
  
  # This is making sure that we only have polygon features
  # No idea why this is necessary, but without it all kinds of geometry errors pop up at clipping
  thiessen_polygons_raw <- sf::st_collection_extract(thiessen_polygons_raw,
                                                     type = "POLYGON")
  
  # Convert the polygons to an sf object
  # Finally, something comfortingly familiar
  thiessen_polygons <- sf::st_sf(thiessen_polygons_raw)
  
  # Clip to sample frame
  thiessen_polygons_clipped <- sf::st_intersection(x = thiessen_polygons,
                                                   y = frame)
  
  # Add in a unique ID for each polygon
  thiessen_polygons_clipped$polygon_unique_id <- seq_len(nrow(thiessen_polygons_clipped))
  
  # Add in the areas for the polygons
  thiessen_polygons_clipped$area_m2 <- as.vector(sf::st_area(x = thiessen_polygons_clipped))
  
  return(thiessen_polygons_clipped)
}

#' @param frame An sf polygon or multipolygon object. This is the clipping boundary which will be applied to the otherise "infinite" Thiessen/Voronoi polygons.
#' @param points An sf point object. These will be broken into clusters and one Thiessen polygon drawn for the centroid of each cluster.
#' @param n_polygons Numeric value. The number of Thiessen polygons to draw within the frame. This is also the number of clusters to break \code{points} into.
#' @param envelope An sfc polygon object. This will be the outer envelope for the Thiessen polygons before they're clipped to \code{frame}. This will only be applied if it's larger than the default envelope in \code{sf::st_voronoi()}. If \code{NULL} then the default envelope will be used. Defaults to \code{NULL}.
#' @param projection Optional character string or CRS object. The coordinate reference system for the output. May be a PROJ4 string or a CRS object. Defaults to the projection of \code{frame}.
#' @param verbose Logical. If \code{TRUE} then the function will return diagnostic messages as it runs. Defaults to \code(FALSE}.)
thiessen_polygons_gen_clustered <- function(frame,
                                            points,
                                            n_polygons,
                                            envelope = NULL,
                                            projection = NULL,
                                            verbose = FALSE) {
  # Sanitization
  if (!("sf" %in% class(frame))) {
    stop("`frame` must be an sf polygon object")
  } else if (!all(sf::st_geometry_type(frame) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("`frame` must be an sf polygon object")
  }
  if (!("sf" %in% class(points))) {
    stop("`points`` must be an sf points object")
  } else if (!all(sf::st_geometry_type(points) %in% c("POINT"))) {
    stop("`points` must be an sf points object")
  }
  if (n_polygons > nrow(points)) {
    stop("`n_polygons` must be less than the number of observations in `points`")
  }
  if (!is.null(envelope)) {
    if (!("sfc" %in% class(envelope))) {
      stop("`envelope` must be an sfc polygon object")
    } else if (!all(sf::st_geometry_type(envelope) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("`envelope` must be an sfc polygon object")
    }
  }
  
  if (is.null(projection)) {
    projection <- sf::st_crs(frame)
  }
  frame <- sf::st_transform(frame,
                            crs = projection)
  points <- sf::st_transform(points,
                             crs = projection)
  
  # Make an envelope from the frame
  if (is.null(envelope)) {
    envelope <- sf::st_as_sfc(sf::st_bbox(frame,
                                          crs = projection))
  }
  
  
  # Get the point coordinates. We'll need them to calculate distances
  points_coords <- as.data.frame(sf::st_coordinates(points))
  names(points_coords) <- c("x", "y")
  
  # Get a distance matrix
  points_distance_matrix <- geosphere::distm(x = points_coords)
  
  # Do some hierarchical clustering based on the distances
  hierarchical_clusters <- hclust(as.dist(m = points_distance_matrix),
                                  method = "complete")
  
  # Put them into a number of clusters matching the Thiessen polygon count
  cluster_membership <- cutree(tree = hierarchical_clusters,
                               k = n_polygons)
  
  # Write that info into the points object
  points$cluster <- cluster_membership
  points_coords$cluster <- cluster_membership
  
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
  
  tpolys <- thiessen_polygons_gen_fixed(centroids = centroids_sf,
                                        frame = frame,
                                        envelope = envelope)
  
  tpolys
}

#' Create Thiessen/Voronoi polygons from a set of points and bounding polygons
#' @description Generate Thiessen/Voronoi polygons for a set of points and clip the results using a set of polygons
#' @param frame An sf polygon or multipolygon object. This is the clipping boundary which will be applied to the otherise "infinite" Thiessen/Voronoi polygons.
#' @param n_polygons Numeric value. The number of Thiessen polygons to draw within the frame.
#' @param points Optional sf point object. If provided, then the Thiessen polygons will be redrawn with new random seeds until each contains at least \code{points_min} of these points. Defaults to \code{NULL}.
#' @param points_min Optional numeric value. If \code{points} is not \code{NULL} then this is the minimum number of points that each Thiessen polygon will contain. Defaults to \code{2}.
#' @param envelope An sfc polygon object. This will be the outer envelope for the Thiessen polygons before they're clipped to \code{frame}. This will only be applied if it's larger than the default envelope in \code{sf::st_voronoi()}. If \code{NULL} then the default envelope will be used. Defaults to \code{NULL}.
#' @param seed_number Optional numeric value. The seed number to use for generating the polygon centroids. A random seed will be used if this is \code{NULL}. Defaults to \code{NULL}.
#' @param seed_increment Optional numeric value. If attempting to produce polygons with \code{points_min} points from \code{points} in each polygon, this is the step to increment \code{seed_number} by on each attempt. Defaults to \code{100000}.
#' @param use_albers Logical. If \code{TRUE} then \code{centroids} and \code{frame} will be reprojected into Albers Equal Area (AEA) and the output will be in AEA. If \code{FALSE} then everything will be reprojected to match the coordinate reference system (CRS) of \code{frame} and the output will be in that CRS. CRSs using decimal degrees will throw errors or warnings. Defaults to \code{TRUE}.
#' @param verbose Logical. If \code{TRUE} then the function will return diagnostic messages as it runs. Defaults to \code(FALSE}.)
#' @return An sf object composed of polygon or multipolygon geometry
thiessen_polygons_gen_random <- function(frame,
                                         n_polygons,
                                         points = NULL,
                                         points_min = 2,
                                         envelope = NULL,
                                         seed_number = NULL,
                                         seed_increment = 100000,
                                         iteration_limit = 1000,
                                         use_albers = TRUE,
                                         verbose = FALSE) {
  # Define Alber's Equal Area CRS
  projection <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  
  # Sanitization
  if (!("sf" %in% class(frame))) {
    stop("`frame` must be an sf polygon object")
  } else if (!all(sf::st_geometry_type(frame) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("`frame` must be an sf polygon object")
  }
  if (!(class(n_polygons) %in% c("numeric", "integer")) | length(n_polygons) > 1) {
    stop("`n_points` must be a single numeric value")
  }
  if (!is.null(points)) {
    if (!("sf" %in% class(points))) {
      stop("`points` must be an sf points object")
    } else if (!all(sf::st_geometry_type(points) %in% c("POINT"))) {
      stop("`points` must be an sf points object")
    }
    if (nrow(points) < (n_polygons * points_min)) {
      stop(paste0("There are too few points available to put ",
                  points_min,
                  " points in each of ",
                  n_polygons,
                  " Thiessen polygons"))
    }
  }
  if (!is.null(seed_number)) {
    if (!(class(seed_number) %in% c("numeric", "integer")) | length(seed_number) > 1) {
      stop("`seed_number` must be a single numeric value")
    }
  } else {
    seed_number <- sample(x = 1:9999999,
                          size = 1)
  }
  
  # Remove any Z dimension
  # It screws with the process and is irrelevant
  frame <- sf::st_zm(frame,
                     drop = TRUE)
  
  # Dissolve!
  frame <- sf::st_as_sf(sf::st_union(x = frame))
  frame$id <- 1
  
  if (!is.null(points)) {
    points <- sf::st_zm(points,
                        drop = TRUE)
  }
  
  if (verbose) {
    message("Making sure the projections are the same")
  }
  # Reproject as necessary
  if (use_albers) {
    frame <- sf::st_transform(x = frame,
                              crs = projection)
    if (!is.null(points)) {
      points <- sf::st_transform(x = points,
                                 crs = projection)
    }
  } else {
    # This just forces the points into the same projection as the polygons
    projection <- sf::st_crs(frame)
    if (!is.null(points)) {
      points <- sf::st_transform(x = points,
                                 crs = projection)
    }
  }
  
  # # Check to make sure we even have enough points for the request
  # if (!is.null(points)) {
  #   if (verbose) {
  #     message("Checking to see if there are enough points for all polygons to meet the minimum point count.")
  #   }
  #   points <- sf::st_intersection(x = points,
  #                                 y = dplyr::select(.data = frame,
  #                                                   id))
  #   if (nrow(points) < min_points * n_polygons) {
  #     stop(paste0("Insufficient points for the number of polygons requested. There are ",
  #                 nrow(points), " in the frame but ", n_polygons, " with at least ", min_points,
  #                 " points each would require at least ", min_points * n_polygons, " points"))
  #   }
  # }
  
  
  
  # Draw centroids
  if (verbose) {
    message("Drawing centroids")
  }
  centroids <- points_gen(frame = frame[, "id"],
                          sample_type = "simple",
                          n_points = n_polygons,
                          seed_number = seed_number,
                          projection = projection)
  
  # Draw Thiessen polygons
  if (verbose) {
    message("Drawing first set of Thiessen polygons")
  }
  thiessen_polygons <- thiessen_polygons_gen_fixed(centroids = centroids,
                                                   frame = frame,
                                                   envelope = envelope,
                                                   use_albers = FALSE)
  
  # Get the final variables in there
  thiessen_polygons$tpoly_seed <- seed_number
  thiessen_polygons$tpoly_id <- paste0("tpoly_",
                                       thiessen_polygons$tpoly_seed,
                                       "-",
                                       thiessen_polygons$polygon_unique_id)
  
  # if (verbose) {
  #   message("Clipping polygons using frame")
  # }
  # thiessen_polygons_clipped <- sf::st_intersection(x = thiessen_polygons,
  #                                                  y = frame)
  
  # We only care about hitting our minimum number of points per polygon if we have points in the first place
  if (!is.null(points)) {
    ## Check that polygons contain enough points
    if (verbose) {
      message("Attributing points with Thiessen polygon IDs via spatial join")
    }
    points_attributed <- sf::st_join(x = points,
                                     y = thiessen_polygons[, c("tpoly_id")])
    
    tpoly_summary <- data.frame(tpoly_id = names(table(points_attributed$tpoly_id)),
                                n_points = as.vector(table(points_attributed$tpoly_id)),
                                stringsAsFactors = FALSE)
    
    # But what if there were polygons with no points at all?????
    missing_tpoly_ids <- unique(thiessen_polygons$tpoly_id)[!(unique(thiessen_polygons$tpoly_id) %in% tpoly_summary$tpoly_id)]
    
    if (length(missing_tpoly_ids) > 0) {
      missing_tpoly_summary <- data.frame(tpoly_id = missing_tpoly_ids,
                                          n_points = 0,
                                          stringsAsFactors = FALSE)
      
      tpoly_summary <- rbind(tpoly_summary,
                             missing_tpoly_summary)
    }
    
    current_iteration <- 1
    # So, if the polygons didn't have enough points each, increment the seed number and try again
    # over and over until it actually pans out
    while (!all(tpoly_summary[["n_points"]] >= points_min)) {
      current_iteration <- current_iteration + 1
      if (current_iteration > iteration_limit) {
        warning("Iteration limit reached without a solution. Returning NULL.")
        return(NULL)
      }
      seed_number <- seed_number + seed_increment
      if (verbose) {
        message(paste("Not enough points in all thiessen polygons. Drawing new centroids with seed", seed_number))
      }
      
      # Draw centroids
      centroids <- points_gen(frame = frame,
                              sample_type = "simple",
                              n_points = n_polygons,
                              seed_number = seed_number,
                              projection = projection)
      
      if (verbose) {
        message("Generating Thiessen polygons from centroids")
      }
      thiessen_polygons <- thiessen_polygons_gen_fixed(centroids = centroids,
                                                       frame = frame,
                                                       envelope = envelope,
                                                       use_albers = FALSE)
      
      
      # Get the final variables in there
      thiessen_polygons$tpoly_seed <- seed_number
      thiessen_polygons$tpoly_id <- paste0("tpoly_",
                                           thiessen_polygons$tpoly_seed,
                                           "-",
                                           thiessen_polygons$polygon_unique_id)
      
      ## Check that polygons contain enough points
      if (verbose) {
        message("Attributing points with Thiessen polygon IDs via spatial join")
      }
      points_attributed <- sf::st_join(x = points,
                                       y = thiessen_polygons[, c("tpoly_id")])
      
      tpoly_summary <- data.frame(tpoly_id = names(table(points_attributed$tpoly_id)),
                                  n_points = as.vector(table(points_attributed$tpoly_id)),
                                  stringsAsFactors = FALSE)
      # But what if there were polygons with no points at all?????
      missing_tpoly_ids <- unique(thiessen_polygons$tpoly_id)[!(unique(thiessen_polygons$tpoly_id) %in% tpoly_summary$tpoly_id)]
      
      if (length(missing_tpoly_ids) > 0) {
        missing_tpoly_summary <- data.frame(tpoly_id = missing_tpoly_ids,
                                            n_points = 0,
                                            stringsAsFactors = FALSE)
        
        tpoly_summary <- rbind(tpoly_summary,
                               missing_tpoly_summary)
      }
    }
    
    # Add in the point counts while we're here
    output <- merge(thiessen_polygons,
                    tpoly_summary,
                    by = "tpoly_id")
    
    # And why not weights also
    output$weight <- output$area_m2 / output$n_points
    
    return(output)
  } else {
    # Return the polygons
    # return(thiessen_polygons_clipped[, c("tpoly_id", "tpoly_seed")])
    return(thiessen_polygons[, c("tpoly_id", "tpoly_seed")])
  }
}

#' Generate a landscape raster with categorical values
#' @param categories Numeric vector. The possible categories (numeric values) that the raster cells will have.
#' @param ncol Numeric value. The number of columns the raster will have.
#' @param nrow Numeric value. The number of rows the raster will have.
#' @param seed_number Optional numeric value. The seed number to use for generating the raster. A random seed will be used if this is \code{NULL}. Defaults to \code{NULL}.
#' @param projection Optional character string or CRS object. The coordinate reference system for the raster. May be a PROJ4 string or a CRS object. Defaults to Albers Equal Area.
#' @returns Raster object of the specified dimensions with cells randomly assigned the values in \code{categories}
landscape_gen_categorical <- function(categories,
                                      ncol,
                                      nrow,
                                      seed_number = NULL,
                                      projection = NULL) {
  if (class(ncol) != "numeric" | length(ncol) > 1) {
    stop("`ncol` must be a single numeric value")
  }
  if (class(nrow) != "numeric" | length(nrow) > 1) {
    stop("`nrow` must be a single numeric value")
  }
  
  if (!is.null(seed_number)) {
    if (!(class(seed_number) %in% c("numeric", "integer")) | length(seed_number) > 1) {
      stop("`seed_number` must be a single numeric value")
    }
  } else {
    seed_number <- sample(x = 1:9999999,
                          size = 1)
  }
  
  if (is.null(projection)) {
    # Define Alber's Equal Area CRS
    message("No projection provided. Defaulting to Albers Equal Area.")
    projection <- sp::CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  } else if (class(projection) == "character") {
    projection <- sp::CRS(projection)
  } else if (!("CRS" %in% class(projection))) {
    stop("`projection`` must either be a valid PROJ4 string or CRS object")
  }
  
  # Set the seed number
  set.seed(seed_number)
  
  # Generate a vector of values long enough to fill the raster
  raster_values <- sample(x = categories,
                          size = ncol * nrow,
                          replace = TRUE)
  
  # Create a matrix from that vector
  landscape_matrix <- matrix(data = raster_values,
                             nrow = nrow,
                             ncol = ncol)
  
  # Convert the matrix into a raster
  landscape_raster <- raster::raster(x = landscape_matrix,
                                     xmn = 0,
                                     xmx = ncol,
                                     ymn = 0,
                                     ymx = nrow,
                                     crs = projection)
  
  return(landscape_raster)
}

#' Generate a landscape with continuous values
#' @param max Numeric. The maximum value possible for a raster cell.
#' @param min Numeric. The minimum value possible for a raster cell.
#' @param distribution Character string. Determines the type of distribution of values in the raster. Must be either \code{"normal"} or \code{"uniform"}. Defaults to \code{"normal"}.
#' @param ncol Numeric value. The number of columns the raster will have.
#' @param nrow Numeric value. The number of rows the raster will have.
#' @param mean Optional numeric value. The mean for the distribution if \code{distribution} is \code{"normal"}. If not provided for a normal distribution, defaults to the mean of \code{max} and \code{min}.
#' @param sd Optional numeric value. The standard deviation for the distribution if \code{distribution} is \code{"normal"}. If not provided for a normal distribution, defaults to \code{1}.
#' @param seed_number Optional numeric value. The seed number to use for generating the raster. A random seed will be used if this is \code{NULL}. Defaults to \code{NULL}.
#' @param projection Optional character string or CRS object. The coordinate reference system for the raster. May be a PROJ4 string or a CRS object. Defaults to Albers Equal Area.
landscape_gen_continuous <- function(max,
                                     min,
                                     distribution = "normal",
                                     ncol,
                                     nrow,
                                     mean = NULL,
                                     sd = NULL,
                                     seed_number = NULL,
                                     projection = NULL) {
  if (class(max) != "numeric" | length(max) > 1) {
    stop("`max` must be a single numeric value")
  }
  if (class(min) != "numeric" | length(min) > 1) {
    stop("`min` must be a single numeric value")
  }
  if (max < min) {
    stop("`max` must be greater than min")
  }
  
  if (!(distribution %in% c("normal", "uniform"))) {
    stop("`distribution` must be either 'normal' or 'uniform'")
  }
  
  if (distribution == "normal") {
    if (is.null(mean)) {
      message("No mean value provided. Using the mean of max and min")
      mean <- (max + min) / 2
    }
    if (is.null(sd)) {
      message("No standard deviation provided. Using a value of 1")
    }
    sd <- 1
  }
  
  if (class(ncol) != "numeric" | length(ncol) > 1) {
    stop("`ncol` must be a single numeric value")
  }
  if (class(nrow) != "numeric" | length(nrow) > 1) {
    stop("`nrow` must be a single numeric value")
  }
  
  if (!is.null(seed_number)) {
    if (!(class(seed_number) %in% c("numeric", "integer")) | length(seed_number) > 1) {
      stop("`seed_number` must be a single numeric value")
    }
  } else {
    seed_number <- sample(x = 1:9999999,
                          size = 1)
  }
  
  if (is.null(projection)) {
    # Define Alber's Equal Area CRS
    message("No projection provided. Defaulting to Albers Equal Area.")
    projection <- sp::CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  } else if (class(projection) == "character") {
    projection <- sp::CRS(projection)
  } else if (!("CRS" %in% class(projection))) {
    stop("`projection` must either be a valid PROJ4 string or CRS object")
  }
  
  # Set the seed number if specified
  set.seed(seed_number)
  
  # Generate a vector of values long enough to fill the raster
  # This will do it according to the requested distribution method
  raster_values <- switch(distribution,
                          "normal" = {
                            values <- rnorm(n = ncol * nrow,
                                            mean = mean,
                                            sd = sd)
                            values[values > max] <- max
                            values[values < min] <- min
                            return(values)
                          },
                          "uniform" = {
                            runif(n = ncol * nrow,
                                  min = min,
                                  max = max)
                          })
  
  # Convert the vector into a matrix
  landscape_matrix <- matrix(data = raster_values,
                             nrow = nrow,
                             ncol = ncol)
  
  # Convert the matrix into a raster
  landscape_raster <- raster::raster(x = landscape_matrix,
                                     xmn = 0,
                                     xmx = ncol,
                                     ymn = 0,
                                     ymx = nrow,
                                     crs = projection)
  
  return(landscape_raster)
}

#' Generate an area of interest
#' @param xmin Numeric value. Minimum possible x coordinate value.
#' @param xmax Numeric value. Maximum possible x coordinate value.
#' @param ymin Numeric value. Minimum possible y coordinate value.
#' @param ymax Numeric value. Maximum possible y coordinate value.
#' @param n_vertices Numeric value. The number of vertices to initially draw. If \code{convex_hull} is \code{TRUE} then this may not be the final number of vertices. Defaults to \code{3}.
#' @param convex_hull Logical. If \code{TRUE} then the function \code{sf::st_convex_hull()} will be applied to the polygon to ensure that it is convex and does not cross itself. Defaults to \code{TRUE}.
#' @param seed_number Optional numeric value. The seed number to use for generating the raster. A random seed will be used if this is \code{NULL}. Defaults to \code{NULL}.
#' @param projection Optional character string or CRS object. The coordinate reference system for the output. May be a PROJ4 string or a CRS object. Defaults to Albers Equal Area.
#' @returns An sf polygon object with a variable called \code{seed_number} containing the seed number the polygon was generated from.

aoi_gen <- function(xmin,
                    xmax,
                    ymin,
                    ymax,
                    n_vertices = 3,
                    convex_hull = TRUE,
                    seed_number = NULL,
                    projection = NULL) {
  if (class(xmax) != "numeric" | length(xmax) > 1) {
    stop("`xmax` must be a single numeric value")
  }
  if (class(xmin) != "numeric" | length(xmin) > 1) {
    stop("`xmin` must be a single numeric value")
  }
  if (xmax < xmin) {
    stop("`xmax` must be greater than xmin")
  }
  if (class(ymax) != "numeric" | length(ymax) > 1) {
    stop("`ymax` must be a single numeric value")
  }
  if (class(ymin) != "numeric" | length(ymin) > 1) {
    stop("`ymin` must be a single numeric value")
  }
  if (ymax < ymin) {
    stop("`ymax` must be greater than ymin")
  }
  
  if (!is.null(seed_number)) {
    if (!(class(seed_number) %in% c("numeric", "integer")) | length(seed_number) > 1) {
      stop("`seed_number` must be a single numeric value")
    }
  } else {
    seed_number <- sample(x = 1:9999999,
                          size = 1)
  }
  
  if (is.null(projection)) {
    # Define Alber's Equal Area CRS
    message("No projection provided. Defaulting to Albers Equal Area.")
    projection <- sp::CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  } else if (class(projection) == "character") {
    projection <- sp::CRS(projection)
  } else if (!("CRS" %in% class(projection))) {
    stop("`projection` must either be a valid PROJ4 string or CRS object")
  }
  
  # Set the seed for x coordinates first
  set.seed(seed_number)
  
  # Generate the x coordinates
  x_coords <- runif(n = n_vertices,
                    min = xmin,
                    max = xmax)
  # Add the first x coordinate to the end to make this loop around so we can have a closed polygon
  x_coords <- c(x_coords, x_coords[1])
  
  # Now do the same for y
  # We're changing the seed number by 1 so that we get different coordinates than for x
  set.seed(seed_number + 1)
  
  y_coords <- runif(n = n_vertices,
                    min = ymin,
                    max = ymax)
  
  y_coords <- c(y_coords, y_coords[1])
  
  # Combine into a data frame to convert into an sf object
  coords <- data.frame(x_coord = x_coords,
                       y_coord = y_coords)
  
  # Make an sfc object from the coordinates
  aoi_sfc <- sf::st_sfc(sf::st_polygon(x = list(outer = as.matrix(coords))))
  
  aoi_sf <- sf::st_sf(aoi_sfc,
                      crs = projection)
  
  aoi_sf$aoi_id <- paste0("aoi_",
                          seed_number)
  
  aoi_sf$aoi_seed <- seed_number
  
  if (convex_hull) {
    aoi_sf <- sf::st_convex_hull(aoi_sf)
  }
  
  return(aoi_sf[, c("aoi_id", "aoi_seed")])
}

#' Generate stratification polygons
#' @description Generate stratification polygons for a frame. If using \code{type = "partitioned"} and providing \code{landscape_raster} then the strata polygons will represent similar areas of the raster. If \code{type = "random"} then the strata will be randomly drawn Thiessen polygons within the frame.
#' @param frame sf POLYGON or MULTIPOLYGON object. The maximum extent of the strata to return.
#' @param n_strata Numeric. The number of strata to produce.
#' @param type Character string. Must be either \code{"partitioned"} or \code{"random"}. If code\{"partitioned"} then \code{landscape_raster} must be provided. Defaults to \code{"partitioned"}.
#' @param landscape_raster Optional RasterLayer object. Required if \code{type = "partitioned"}. Must have a value range of 0 to 1. Defaults to \code{"NULL"}.
#' @param seed_number Optional numeric value. The seed number to use for generating the polygon stratum centroids if \code{type = "random"}. A random seed will be used if this is \code{NULL}. Defaults to \code{NULL}.
#' @param seed_increment Optional numeric value. If attempting to produce polygons with \code{points_min} points from \code{points} in each polygon, this is the step to increment \code{seed_number} by on each attempt. Defaults to \code{100000}.
strata_gen <- function(frame,
                       n_strata = 2,
                       type = "partitioned",
                       cluster = TRUE,
                       landscape_raster = NULL,
                       seed_number = NULL,
                       seed_increment = 10000) {
  if (!(c("sf") %in% class(frame))) {
    stop("`frame` must be an sf polygon object")
  } else if (!(sf::st_geometry_type(frame) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("`frame` must be an sf polygon object.")
  }
  
  if (n_strata < 1) {
    stop("`strata` must be a number greater than 1.")
  }
  
  if (!(type %in% c("partitioned", "random"))) {
    stop("`type` must be either 'partitioned' or 'random'.")
  }
  
  if (!is.null(landscape_raster)) {
    if (class(landscape_raster) != "RasterLayer") {
      stop("`raster` must be of class RasterLayer.")
    }
  } else if (type == "partitioned") {
    stop("`raster` must be provided when `type` = 'partitioned'")
  }
  
  if (!is.null(seed_number)) {
    if (!(class(seed_number) %in% c("numeric", "integer")) | length(seed_number) > 1) {
      stop("`seed_number` must be a single numeric value")
    }
  } else {
    seed_number <- sample(x = 1:9999999,
                          size = 1)
  }
  
  strata_polygons <- switch(type,
                            "partitioned" = {
                              # Get a raster to mangle
                              landscape_categorical <- landscape_raster
                              
                              if (cluster) {
                                # Figure out where the breaks are for the partitions
                                stratum_breaks <- BAMMtools::getJenksBreaks(var = as.vector(landscape_raster),
                                                                            k = n_strata + 1)
                                
                                for (category in length(stratum_breaks):2) {
                                  # Get the upper and lower cutoff values
                                  current_stratum_max <- stratum_breaks[category]
                                  current_stratum_min <- stratum_breaks[category - 1]
                                  
                                  landscape_categorical[landscape_categorical >= current_stratum_min & landscape_categorical <= current_stratum_max] <- category - 1
                                }
                              } else {
                                # What's the value range for strata?
                                strata_increment <- 1 / n_strata
                                
                                # Convert from continuous values into categories
                                # Note that we go from largest to smallest so that we don't accidentally overwrite
                                # category 1 when we get to the highest number category which has an upper limit of 1
                                for (category in n_strata:1) {
                                  current_stratum_min <- (category - 1) * strata_increment
                                  current_stratum_max <- category * strata_increment
                                  landscape_categorical[landscape_categorical >= current_stratum_min & landscape_categorical <= current_stratum_max] <- category
                                }
                              }
                              
                              landscape_categorical_spdf <- raster::rasterToPolygons(x = landscape_categorical,
                                                                                     dissolve = TRUE)
                              landscape_categorical_sf <- methods::as(landscape_categorical_spdf,
                                                                      "sf")
                              landscape_categorical_sf$stratum_id <- landscape_categorical_sf$layer
                              
                              strata_sf <- sf::st_intersection(x = landscape_categorical_sf,
                                                               y = frame)
                              
                              strata_sf[, c("stratum_id")]
                            },
                            "random" = {
                              
                              boundary_spdf <- methods::as(frame,
                                                           "Spatial")
                              boundary_sfc <- sf::st_as_sfc(boundary_spdf)
                              
                              
                              
                              strata_polygons <- thiessen_polygons_gen_random(frame = frame,
                                                                              n_polygons = n_strata,
                                                                              envelope = boundary_sfc,
                                                                              seed_number = seed_number)
                              
                              # Not necessarily all the strata will be represented in the frame, so we'll
                              # keep drawing until we get a version where that happens
                              not_enough_strata <- nrow(strata_polygons) < n_strata
                              
                              while (not_enough_strata) {
                                seed_number <- seed_number + seed_increment
                                
                                message(paste0("Not enough strata fell inside the frame. Attempting again with seed number ",
                                               seed_number))
                                
                                strata_polygons <- thiessen_polygons_gen_random(frame = frame,
                                                                                n_polygons = n_strata,
                                                                                envelope = raster_boundary_sfc,
                                                                                seed_number = seed_number)
                                
                                not_enough_strata <- nrow(strata_polygons) < n_strata
                              }
                              
                              strata_polygons$stratum_id <- 1:n_strata
                              
                              strata_polygons[, c("stratum_id")]
                            }
  )
  
  strata_polygons
}

#' Make a single-part polygon (more) concave
#' @param polygon An sf polygon object. The must not be multipart or have multiple polygons. It may already be concave.
#' @param seed_number Optional numeric value. The seed number used to generate the coordinates of the new vertex. If \code{NULL} then a random number will be used. Defaults to \code{NULL}.
#' @returns An sf polygon object containing the original polygon, but with a new vertex that adds a concavity to it.
concavify <- function(polygon,
                      seed_number = NULL) {
  if (!("sf" %in% class(polygon))) {
    stop("`polygon` must be an sf object of geometry type 'POLYGON'")
  } else if (sf::st_geometry_type(polygon) != "POLYGON") {
    stop("`polygon` must be an sf object of geometry type 'POLYGON'")
  }
  
  if (!is.null(seed_number)) {
    if (!(class(seed_number) %in% c("numeric", "integer")) | length(seed_number) > 1) {
      stop("`seed_number` must be a single numeric value")
    }
  } else {
    seed_number <- sample(x = 1:9999999,
                          size = 1)
  }
  
  # Get the polygon's projection. We'll use this for the vertices
  projection <- sf::st_crs(polygon)
  
  # Get the bounding box for the polygon
  polygon_bbox <- sf::st_bbox(polygon)
  
  # Generate potential new vertices within the bounding box
  # There's a small chance that none of them will fall in the polygon!
  set.seed(seed_number)
  x_coords <- runif(n = 100,
                    min = polygon_bbox["xmin"],
                    max = polygon_bbox["xmax"])
  set.seed(seed_number + 1)
  y_coords <- runif(n = 100,
                    min = polygon_bbox["ymin"],
                    max = polygon_bbox["ymax"])
  
  # Mash up the x and y coordinates into a single data frame to use
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
  
  # Prepare for a worst-case bad luck scenario
  if (!any(concave_vertex_indices)) {
    stop("There were no valid vertices for that seed number. Try another!")
  }
  
  # Keep only the inside-the-polygon points
  valid_concave_vertices <- potential_vertices[concave_vertex_indices, ]
  
  # Get the polygon coordinates
  polygon_coordinates <- sf::st_coordinates(polygon)
  
  # And get the potential new vertices' coordinates
  valid_concave_vertex_coordinates <- sf::st_coordinates(valid_concave_vertices)
  
  # Get a random insertion point for the new vertex
  # This means that there will never be an insertion made on on side of the polygon because it excludes the final index as a possibility
  # This is just for ease of programming, sorry!
  set.seed(seed_number)
  insertion_index <- sample(x = 1:(nrow(polygon_coordinates) - 1),
                            size = 1)
  
  # Make new vectors for x and y coordinates inserting the first valid concave vertex into the vector after insertion_index
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
  concave_polygon <- sf::st_sf(sf::st_cast(sf::st_combine(new_vertices),
                                           "POLYGON"))
  
  # Write in the original properties
  for (variable in names(polygon)) {
    concave_polygon[[variable]] <- polygon[[variable]]
  }
  
  # Make sure the output is only the original properties
  output <- concave_polygon[, names(polygon)]
  
  return(output)
}

#' Generate weight category polygons
#' @param polygons sf polygon object. Must contain each frame as an observation.
#' @param aoi_index Optional numeric value. The index in \code{polygons} of the polygon acting as the area of interest or boundary for the intersection. The function will only return geometry that overlaps with the polygon at this index, discarding all geometry that consists only of the other polygons. If \code{NULL} then no filtering will happen. Defaults to \code{NULL}.
#' @returns sf polygon object with the variables "wgtcat_id" containing the unique identifier for the weight categories and "area_m2" containing the area in meters squared (assuming that the projection is Albers Equal Area or other with meters as units).
wgtcat_gen <- function(polygons,
                       intersect = TRUE,
                       precision = NULL,
                       aoi_index = NULL) {
  if (!("sf" %in% class(polygons))) {
    stop("`polygons` must be an sf object of geometry type 'POLYGON'")
  } else if (!all(sf::st_geometry_type(polygons) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("`polygons` must be an sf object of geometry type 'POLYGON'")
  }
  if (!is.null(aoi_index)) {
    if (!(aoi_index %in% seq_len(nrow(polygons)))) {
      stop("`aoi_index` must refer to the index of one of the entries in polygons")
    }
  }
  
  if (!is.null(precision)) {
    if (is.numeric(precision)) {
      polygons <- sf::st_set_precision(x = polygons,
                                       precision = precision)
    } else {
      stop("precision must either be NULL or a numeric value.")
    }
  }
  
  
  if (intersect) {
    # Intersect them to create weight category polygons
    polygons <- repair_geometry(polygons = polygons)
    
    wgtcat_inside_intensification <- sf::st_intersection(x = dplyr::filter(polygons,
                                                                           uid != "intensification"),
                                                         y = dplyr::filter(polygons,
                                                                           uid == "intensification"))
    wgtcat_inside_intensification$wgtcat_id <- paste0(wgtcat_inside_intensification$uid,
                                                      "-intensification")
    
    wgtcat_outside_intensification <- sf::st_difference(x = dplyr::filter(polygons,
                                                                          uid != "intensification"),
                                                        y = dplyr::filter(polygons,
                                                                          uid == "intensification"))
    wgtcat_outside_intensification$wgtcat_id <- wgtcat_outside_intensification$uid
    
    
    wgtcat_polygons <- dplyr::bind_rows(wgtcat_outside_intensification,
                                        wgtcat_inside_intensification)
    # Add in a unique ID for them
    # wgtcat_polygons$wgtcat_id <- sapply(X = wgtcat_polygons$origins,
    #                                     FUN = function(X) {
    #                                       paste(X, collapse = "-")
    #                                     })
  } else {
    wgtcat_polygons <- polygons
    wgtcat_polygons$wgtcat_id <- paste0(1:nrow(polygons))
  }
  
  # Add the areas
  wgtcat_polygons$area_m2 <- as.vector(sf::st_area(wgtcat_polygons))
  
  if (any(!(sf::st_geometry_type(wgtcat_polygons) %in% c("MULTIPOLYGON", "POLYGON")))) {
    stop("Hey! Some kind of geometry error(s) occurred in generating the wgtcat polygons that produced non-polygon geometry, which is bad.")
  }
  
  return(wgtcat_polygons)
}

#' Generate points within a given frame by one of three methods
#' @param frame sf polygon object. The sample frame within which points will be drawn. This should probably be a single simple polygon.
#' @param sample_type Character string. The method to draw points by. Valid values are \code{"simple"} (simple random), \code{"balanced"} (spatially-balanced random using GRTS), and \code{"cluster} (two-stage clustered). Defaults to \code{"simple"}.
#' @param n_points Numeric value. The number of points to draw in the sample frame. 
#' @param seed_number Optional numeric value. The seed used to generate random points. A random seed will be used if this is \code{NULL}. Defaults to \code{NULL}.
#' @param projection Optional character string or CRS object. The coordinate reference system for the points. May be a PROJ4 string or a CRS object. Defaults to the projection of \code{frame}.
#' @returns An sf points object with the variables \code{"sample_id"} containing the unique identifiers for the points and \code{sample_seed} containing the random seed number used.
points_gen <- function(frame,
                       sample_type = "simple",
                       n_points,
                       seed_number = NULL,
                       projection = NULL) {
  if (!("sf" %in% class(frame))) {
    stop("`frame` must be an sf object of geometry type 'POLYGON' or 'MULTIPOLYGON'")
  } else if (!any(c("POLYGON", "MULTIPOLYGON") %in% sf::st_geometry_type(frame))) {
    stop("`frame` must be an sf object of geometry type 'POLYGON' or 'MULTIPOLYGON'")
  }
  
  if (!(class(n_points) %in% c("numeric", "integer")) | length(n_points) > 1) {
    stop("`n_points` must be a single numeric value")
  }
  
  if (!(sample_type %in% c("simple", "balanced", "cluster"))) {
    stop("`sample_type` must be 'simple', 'balanced', or 'cluster'")
  }
  
  if (!is.null(seed_number)) {
    if (!(class(seed_number) %in% c("numeric", "integer")) | length(seed_number) > 1) {
      stop("`seed_number` must be a single numeric value")
    }
  } else {
    seed_number <- sample(x = 1:9999999,
                          size = 1)
  }
  
  if (is.null(projection)) {
    projection <- sf::st_crs(frame)
  } else if (class(projection) == "character") {
    projection <- sp::CRS(projection)
  } else if (!("CRS" %in% class(projection))) {
    stop("`projection` must either be a valid PROJ4 string or CRS object")
  }
  # Reproject
  frame <- sf::st_transform(x = frame,
                            crs = projection)
  
  # Do the correct kind of sample draw
  points <- switch(sample_type,
                   "simple" = {
                     set.seed(seed_number)
                     sp::spsample(x = methods::as(frame, "Spatial"),
                                  n = n_points,
                                  type = "random",
                                  iter = 10)
                   },
                   "balanced" = {
                     set.seed(seed_number)
                     spsurvey::grts(sframe = frame,
                                    n_base = n_points,
                                    n_over = 0,
                                    DesignID = "design",
                                    seltype = "equal",
                                    sep = "-")
                   },
                   "cluster" = {
                     set.seed(seed_number)
                     
                   })
  
  # Convert to sf
  if (!("sf" %in% class(points))) {
    points <- methods::as(points, "sf")
  }
  
  # Add our ID in the format *I* want
  points$sample_id <- paste0("sample_",
                             seed_number,
                             "-",
                             seq_len(nrow(points)))
  
  points$sample_seed <- seed_number
  
  output <- points[, c("sample_id", "sample_seed")]
  
  return(output)
}

#' Calculate the upper and lower bounds for a mean given and alpha value
#' @param mean Numeric value. The mean to compute bounds for.
#' @param sd Numeric value. The standard deviation of \code{mean}.
#' @param n Numeric value. The number of observations that were used to calculate \code{mean}.
#' @param alpha Numeric value. The alpha value to use to compute the upper and lower confidence bounds for \code{mean}. Defaults to \code{0.05}.
#' @returns A named list containing the upper and lower bounds for the mean for the given confidence.
ci_mean <- function(mean,
                    sd,
                    n,
                    alpha) {
  if (!is.numeric(mean)) {
    stop("`mean` must be a numeric value")
  }
  if (!is.numeric(sd)) {
    stop("`sd` must be a numeric value")
  }
  if (!is.numeric(n)) {
    stop("`n` must be a numeric value")
  } else if (n <= 1) {
    stop("`n` must be greater than 1")
  }
  if (!is.numeric(alpha)) {
    stop("`alpha` must be a numeric value")
  } else if (alpha <= 0 | alpha >= 1) {
    stop("`alpha` must be a value between 0 and 1")
  }
  
  standard_error <- sd / sqrt(n)
  degrees_freedom <- n - 1
  t_score <- qt(p = alpha / 2,
                df = degrees_freedom)
  margin_error <- abs(standard_error * t_score)
  mean_bound_lower <- mean - margin_error
  mean_bound_upper <- mean + margin_error
  
  list(lower_bound = mean_bound_lower,
       upper_bound = mean_bound_upper)
}

# Literally just from https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
weighted_mean_se <- function(values,
                             weights){
  normalized_weights <- weights / sum(weights)
  variance <- var(x = values)
  sqrt(variance) * sqrt(sum(normalized_weights^2))
}

#' Calculate Goodman's multinomial confidence intervals
#' @description Calculate confidence intervals for multinomial proportions using the method described by Leo Goodman in "On Simultaneous Confidence Intervals for Multinomial Proportions" in Technometrics in 1965. This function can only handle one group of categorical counts at a time, so if you want to calculate confidence intervals for multiple groups, you need to do each separately.
#' @param counts Numeric vector, optionally named. The counts for each of the categories being considered. If there are unequal weights, be sure to adjust these counts by proportional weight with the formula: adjusted count for a category = total observations * sum of weights of observations in the category / sum of all weights. If these values are named, those will be included in the output data frame.
#' @param alpha Numeric value. Must be between 0 and 1. The alpha for the confidence calculation, e.g. for 80 percent confidence, the alpha is 0.2. Defaults to \code{0.2}.
#' @param chisq Character string. This decides which chi squared quantile calculation to use. The accepted values are \code{"A"}, \code{"B"}, or \code{"best"} (use the one which minimizes the confidence intervals). Goodman describes A as his default, calculated as the upper alpha times 100th percentage point of the chi-square distribution with k - 1 degrees of freedom. He also notes the alternative B, calculated as the upper alpha / k times 100th percentage point of the chi-square distribution with one degree of freedom, which will produce tighter intervals when k > 2 and alpha is 0.1, 0.5, or 0.01. Defaults to \code{"best"}
#' @param verbose Logical. If \code{TRUE} then the function will generate additional messages as it executes. Defaults to \code{FALSE}.
#' @export
goodman_cis <- function(counts,
                        alpha = 0.2,
                        chisq = "best",
                        verbose = FALSE) {
  if (!is.numeric(counts) | length(counts) < 2) {
    stop("`counts` must be a numeric vector with at least two values")
  }
  
  if (!(chisq %in% c("A", "B", "best"))) {
    stop("`chisq` must be 'A', 'B', or 'best'.")
  }
  
  # Goodman describes the upper and lower bounds with the equations:
  # Lower estimated pi_i = {A + 2n_i - {A[A + 4n_i(N - n_i) / N]}^0.5} / [2(N + A)]
  # Upper estimated pi_i = {A + 2n_i + {A[A + 4n_i(N - n_i) / N]}^0.5} / [2(N + A)]
  
  # n_i is the "observed cell frequencies in population of size N" (aka count of observations) from a category
  # so that's the incoming argument counts. We'll rename for consistency with the original math (and statistics as a discipline)
  n <- counts
  
  # N is the population those counts make up, or, in lay terms, the total observation count
  N <- sum(counts)
  
  # k is the number of categories the population has been sorted into
  # Useful for degrees of freedom
  k <- length(counts)
  
  # "A is the upper alpha * 100-th percentage point of the chi-square distribution with k - 1 degrees of freedom"
  # and B is an alternative which uses alpha / k and one degree of freedom
  # Goodman states that B should be less than A for situations
  # where k > 2 AND alpha is 0.1, 0.05, or 0.01.
  chisq_quantiles <- c("A" = stats::qchisq(p = 1 - alpha,
                                           df = k - 1),
                       "B" = stats::qchisq(p = 1 - (alpha / k),
                                           df = 1))
  
  
  # According to Goodman, A and B are both valid options for the chi-square quantile
  # So the user can specify which they want or just ask for the one that minimizes the confidence intervals
  chisq_quantile <- switch(chisq,
                           "A" = {
                             chisq_quantiles["A"]
                           },
                           "B" = {
                             chisq_quantiles["B"]
                           },
                           "best" = {
                             pick <- which.min(chisq_quantiles)
                             if (verbose) {
                               switch(names(chisq_quantiles)[pick],
                                      "A" = message("The chi-square quantile calculation that will provide the tighter confidence intervals is A, the upper alpha X 100-th percentage point of the chi-square distribution with k - 1 degrees of freedom"),
                                      "B" = message("The chi-square quantile calculation that will provide the tighter confidence intervals is B, the upper alpha / k X 100-th percentage point of the chi-square distribution with 1 degree of freedom"))
                             }
                             chisq_quantiles[pick]
                           })
  
  # Calculate the bounds!
  # Note that these ARE symmetrical, just not around the proportions.
  # They're symmetrical around A + 2 * n / (2 * (N + A))
  # The variable A has been replaced with chisq_quantile because it may be A or B, depending
  # Since the only multi-value vector involved here is n, these will be vectors of length k,
  # having one value for each of the values in n and in the same order as n
  lower_bounds <- (chisq_quantile + 2 * n - sqrt(chisq_quantile * (chisq_quantile + 4 * n * (N - n) / N))) / (2 * (N + chisq_quantile))
  upper_bounds <- (chisq_quantile + 2 * n + sqrt(chisq_quantile * (chisq_quantile + 4 * n * (N - n) / N))) / (2 * (N + chisq_quantile))
  
  # A proportion can never be greater than 1 or less than 0 (duh)
  # So we'll add bounds any CIs in case that happens
  # That's definitely a thing that can happen if the magnitude of sqrt(A * (A + 4 * n * (N - n) / N))
  # is large enough
  lower_bounds[lower_bounds < 0] <- 0
  upper_bounds[upper_bounds > 1] <- 1
  
  # Build the output
  output <- data.frame(count = n,
                       proportion = n / N,
                       lower_bound = lower_bounds,
                       upper_bound = upper_bounds,
                       stringsAsFactors = FALSE,
                       row.names = NULL)
  
  # What are the categories called? If anything, that is
  k_names <- names(n)
  
  if (!is.null(k_names)) {
    output[["category"]] <- k_names
    output <- output[, c("category", "count", "proportion", "lower_bound", "upper_bound")]
  }
  
  return(output)
}

#' Estimation of weighted proportions of categorical data
#' @description Given categorical data and the weights for the individual observations, calculate estimated proportions by category and Goodman's multinomial confidence intervals.
#' @param data Data frame. Categorical data with the unique identifiers for each observation/row in the variable \code{id_var} and the assigned category for each observation/row in \code{cat_var}. Note that the unique identifiers are the link between \code{data} and \code{weights}
#' @param weights Data frame. This must contain the weighting information using the variables \code{id_var} with a unique identifier for each observation/row and \code{wgt_var} with the relative numeric weight of each observation/row.
#' @param id_var Character string. The name of the variable in \code{data} and \code{weights} that contains the unique identifiers for the observations. All values in \code{data$id_var} must appear in \code{weights$id_var}.
#' @param cat_var Character string. The name of the variable in \code{data} that contains the category values as character strings.
#' @param wgt_var Character string. The name of the variable in \code{weights} that contains the numeric weight values.
#' @param definitions Conditionally optional character vector. The possible categories that the observation could've been classed into. This is NOT optional if there are categories that do not appear in \code{data} because no observations met their criteria because those categories must be included in the calculations. Must contain at least the values in \code{code$cat_var} but should include ALL possible categories.
#' @param conf Numeric. The confidence level in percent. Defaults to \code{80}.
#' @param verbose Logical. If \code{TRUE} then the function will generate additional messages as it executes. Defaults to \code{FALSE}.
#' @return A data frame containing the categories, counts of observations, weighted estimated proportions, and confidence intervals.
#' @export
analyze_cat <- function(data,
                        weights,
                        id_var,
                        cat_var,
                        wgt_var,
                        definitions = NULL,
                        conf = 80,
                        verbose = FALSE){
  # Make sure everything is the right class/length
  if (!("data.frame" %in% class(data))) {
    stop("data must be a data frame")
  }
  if (nrow(data) < 1) {
    stop("There are no values in data")
  }
  if (!("data.frame" %in% class(weights))) {
    stop("weights must be a data frame")
  }
  if (nrow(weights) < 1) {
    stop("There are no values in weights")
  }
  
  if (class(id_var) != "character" | length(id_var) != 1) {
    stop("id_var must be a single character string")
  }
  if (class(cat_var) != "character" | length(cat_var) != 1) {
    stop("cat_var must be a single character string")
  }
  if (class(wgt_var) != "character" | length(wgt_var) != 1) {
    stop("wgt_var must be a single character string")
  }
  if (conf <= 0 | conf >= 100) {
    stop("conf must be a value between 0 and 100")
  }
  
  
  # Make sure all the variables are in place
  required_data_vars <- c(id_var,
                          cat_var)
  missing_data_vars <- required_data_vars[!(required_data_vars %in% names(data))]
  if (length(missing_data_vars) > 0) {
    stop("The following variables are missing from data: ", paste(missing_data_vars, collapse = , ", "))
  }
  data <- data[, required_data_vars]
  category_class <- class(data[[cat_var]])
  # What categories were observed?
  present_categories <- unique(data[[cat_var]])
  
  if (!is.null(definitions)) {
    if (!(category_class %in% class(definitions))) {
      stop("definitions must be the same class as the category values in data")
    }
    if (length(definitions) < 1) {
      stop("There are no values in definitions")
    }
  }
  
  # # Check to make sure the unique identifiers are, in fact, unique
  # non_unique_ids <- any(table(data[[id_var]]) > 1)
  # if (any(non_unique_ids)) {
  #   stop("There are non-unique values in ", id_var, " in data.")
  # }
  
  
  required_weights_vars <- c(id_var, wgt_var)
  missing_weights_vars <- required_weights_vars[!(required_weights_vars %in% names(weights))]
  if (length(missing_weights_vars) > 0) {
    stop("The following variables are missing from weights: ", paste(missing_weights_vars, collapse = , ", "))
  }
  non_unique_ids <- any(table(weights[[id_var]]) > 1)
  if (non_unique_ids) {
    stop("There are non-unique values in ", id_var, " in weights.")
  }
  weights <- weights[, required_weights_vars]
  
  # And what if the user provided definitions?
  # This is important for if there are categories that have no data that qualified!
  if (!is.null(definitions)) {
    missing_categories <- !(present_categories %in% definitions)
    if (any(missing_categories)) {
      stop("The following categories appear in data but not in categories: ",
           paste(present_categories[missing_categories], collapse = ", "))
    }
  }
  
  # Make sure the IDs line up
  data_ids_in_weights_indices <- data[[id_var]] %in% weights[[id_var]]
  if (!all(data_ids_in_weights_indices)) {
    stop("Not all unique IDs in data appear in weights")
  }
  weight_ids_in_weights_indices <- weights[[id_var]] %in% data[[id_var]]
  if (verbose & !all(weight_ids_in_weights_indices)) {
    message("Not all unique IDs in weights appear in data, just so you know.")
  }
  weights <- weights[weight_ids_in_weights_indices, ]
  
  
  # Get each observation with just its category and weight
  weighted_categories <- merge(x = data[, c(id_var, cat_var)],
                               y = weights,
                               by = id_var,
                               all.y = FALSE)
  
  # Calculate the sum of the weights for each of the observed categories
  category_weight_sums <- sapply(X = present_categories,
                                 data = weighted_categories,
                                 cat_var = cat_var,
                                 wgt_var = wgt_var,
                                 USE.NAMES = TRUE,
                                 FUN = function(X,
                                                data,
                                                cat_var,
                                                wgt_var){
                                   relevant_indices <- data[[cat_var]] == X
                                   current_weights <- data[relevant_indices, wgt_var]
                                   weight_sum <- sum(as.numeric(current_weights))
                                   return(weight_sum)
                                 })
  # Calculate the weighted proportions for each category
  category_weighted_proportions <- category_weight_sums / sum(category_weight_sums)
  # Get the pure counts of the categories
  category_counts <- table(weighted_categories[[cat_var]])
  # And the total number of observations. This should be the same as nrow(weighted_categories)
  total_observations <- sum(category_counts)
  # Using the total number of observations and the weighted proportions to calculate "adjusted counts"
  adjusted_counts <- category_weighted_proportions * total_observations
  
  # Here's a tricky bit! Calculating weighted standard error, which should be done
  # for each category as well
  category_weighted_se <- sapply(X = present_categories,
                                 data = weighted_categories,
                                 cat_var = cat_var,
                                 wgt_var = wgt_var,
                                 USE.NAMES = TRUE,
                                 FUN = function(X,
                                                data,
                                                cat_var,
                                                wgt_var){
                                   # For each category, we're going to treat that
                                   # category's records as 1 and the others as 0
                                   weighted_mean_se(values = as.numeric(data[[cat_var]] %in% X),
                                                    weights = data[[wgt_var]])
                                 })
  category_weighted_cv <- sapply(X = present_categories,
                                 data = weighted_categories,
                                 cat_var = cat_var,
                                 wgt_var = wgt_var,
                                 USE.NAMES = TRUE,
                                 FUN = function(X,
                                                data,
                                                cat_var,
                                                wgt_var){
                                   # For each category, we're going to treat that
                                   # category's records as 1 and the others as 0
                                   Weighted.Desc.Stat::w.cv(x = as.numeric(data[[cat_var]] %in% X),
                                                            mu = data[[wgt_var]])
                                 })
  
  # Okay, so if we have definitions to catch categories with zero observations, add those
  # Because it should matter for calculating confidence intervals
  if (!is.null(definitions)) {
    defined_categories <- definitions
    missing_categories <- defined_categories[!(defined_categories %in% present_categories)]
    # Looping because it's easy, not because it's the best solution
    # But we want to populate the 0s for all of these!
    for (category in missing_categories) {
      category_weighted_proportions[[category]] <- 0
      category_weight_sums[[category]] <- 0
      category_counts[[category]] <- 0
      adjusted_counts[[category]] <- 0
    }
  }
  
  # Finally ready to calculate confidence intervals!
  # But first we need the alpha value for our confidence level
  alpha <- 1 - (conf / 100)
  
  confidence_intervals <- goodman_cis(counts = adjusted_counts,
                                      alpha = alpha,
                                      chisq = "best",
                                      verbose = verbose)
  confidence_interval_vars <- c("category", "weighted_observation_count", "weighted_observation_proportion",
                                paste0(c("weighted_observation_proportion_lower_bound", "weighted_observation_proportion_upper_bound"),
                                       "_", conf, "pct"))
  names(confidence_intervals) <- confidence_interval_vars
  
  # And now it's a matter of combining and formatting
  # Yeah, yeah, yeah. It's not """best practice""" to calculate within the data frame construction
  # but I don't care. I'll do math and slicing wherever I want to. Deal with it.
  results <- data.frame(category = names(category_counts),
                        observation_count = as.vector(category_counts[names(category_counts)]),
                        observation_proportion = as.vector(category_counts[names(category_counts)] / total_observations),
                        total_observation_weight = category_weight_sums[names(category_counts)],
                        weighted_observation_proportion = category_weighted_proportions[names(category_counts)],
                        weighted_standard_error = category_weighted_se[names(category_counts)],
                        weighted_coefficient_of_variance = category_weighted_cv[names(category_counts)],
                        row.names = NULL,
                        stringsAsFactors = FALSE)
  
  confidence_interval_keep_vars <- c("category",
                                     paste0(c("weighted_observation_proportion_lower_bound", "weighted_observation_proportion_upper_bound"),
                                            "_", conf, "pct"))
  
  # Combine the results and confidence intervals
  output <- merge(x = results,
                  y = confidence_intervals[, confidence_interval_keep_vars],
                  by = c("category"))
  
  # Get the variables restricted to what we care about and ordered properly
  output_vars <- c("category", "observation_count", "observation_proportion", "total_observation_weight", "weighted_observation_proportion", "weighted_standard_error", "weighted_coefficient_of_variance",
                   paste0(c("weighted_observation_proportion_lower_bound", "weighted_observation_proportion_upper_bound"),
                          "_", conf, "pct"))
  
  output <- output[, output_vars]
  
  return(output)
}

analyze_con <- function(data,
                        weights,
                        id_var,
                        value_var,
                        wgt_var,
                        conf = 80,
                        verbose = FALSE){
  # Make sure everything is the right class/length
  if (!("data.frame" %in% class(data))) {
    stop("data must be a data frame")
  }
  if (nrow(data) < 1) {
    stop("There are no values in data")
  }
  if (!("data.frame" %in% class(weights))) {
    stop("weights must be a data frame")
  }
  if (nrow(weights) < 1) {
    stop("There are no values in weights")
  }
  
  if (class(id_var) != "character" | length(id_var) != 1) {
    stop("id_var must be a single character string")
  }
  if (class(value_var) != "character" | length(value_var) != 1) {
    stop("value_var must be a single character string")
  }
  if (class(wgt_var) != "character" | length(wgt_var) != 1) {
    stop("wgt_var must be a single character string")
  }
  if (conf <= 0 | conf >= 100) {
    stop("conf must be a value between 0 and 100")
  } else {
    alpha <- 1 - conf / 100
  }
  
  # Make sure all the variables are in place
  required_data_vars <- c(id_var,
                          value_var)
  missing_data_vars <- required_data_vars[!(required_data_vars %in% names(data))]
  if (length(missing_data_vars) > 0) {
    stop("The following variables are missing from data: ",
         paste(missing_data_vars,
               collapse = , ", "))
  }
  
  # Just want the bare minimum here.
  data <- dplyr::select(.data = data,
                        tidyselect::all_of(required_data_vars)) |>
    dplyr::rename(.data = _,
                  setNames(object = required_data_vars,
                           nm = c("id", "value")))
  
  # Check to make sure the unique identifiers are, in fact, unique
  non_unique_ids <- any(table(data[[id_var]]) > 1)
  if (non_unique_ids) {
    stop("There are non-unique values in ", id_var, " in data.")
  }
  
  required_weights_vars <- c(id_var,
                             wgt_var)
  missing_weights_vars <- required_weights_vars[!(required_weights_vars %in% names(weights))]
  if (length(missing_weights_vars) > 0) {
    stop("The following variables are missing from weights: ",
         paste(missing_weights_vars,
               collapse = , ", "))
  }
  
  non_unique_ids <- any(table(weights[[id_var]]) > 1)
  if (non_unique_ids) {
    stop("There are non-unique values in ", id_var, " in weights.")
  }
  
  # Paring this down too.
  weights <- dplyr::select(.data = weights,
                           tidyselect::all_of(required_weights_vars)) |>
    dplyr::rename(.data = _,
                  setNames(object = required_weights_vars,
                           nm = c("id", "weight")))
  
  if (!all(data[[id_var]] %in% weights[[id_var]])) {
    warning("Not all data have corresponding weights. They will be dropped from the calculations.")
  }
  if (!all(weights[[id_var]] %in% data[[id_var]])) {
    warning("Not all weights have corresponding data. Depending on your situation, this may be expected or may be indicative of an issue with the unique IDs.")
  }
  
  data <- dplyr::inner_join(x = data,
                            y = weights,
                            by = "id",
                            relationship = "one-to-one") |>
    dplyr::mutate(.data = _,
                  weighted_value = value * weight / sum(weights$weight))
  
  n <- nrow(data)
  # Weighted mean is the sum of the weight-adjusted values divided by the sum of all weights
  mean_weighted <- sum(data$weighted_value)
  # Standard deviation is calculated differently for weighted values than unweighted
  sd_weighted <- sqrt(sum(data$weight * (data$value - mean(data$value))^2) / ((n - 1) / n * sum(data$weight)))
  # So is variance
  variance_weighted <- weighted_variance(values = data$value,
                                         weights = data$weight,
                                         na_remove = FALSE)
  bounds_weighted <- ci_mean(mean = mean_weighted,
                             sd = sd_weighted,
                             n = n,
                             alpha = alpha)
  
  data.frame(n = n,
             alpha = alpha,
             mean = mean_weighted,
             sd = sd_weighted,
             cv = sd_weighted / mean_weighted,
             variance = variance_weighted,
             lower_bound = bounds_weighted$lower_bound,
             upper_bound = bounds_weighted$upper_bound)
}

# This is literally only here for the dang bootstrapping
special_mean <- function(data, indices) {
  mean(data[indices],
       trim = 0)
}

#' Weighted analysis of a continuous variable
#' @param data Any object that can be treated as a data frame. Must contain the variables \code{value} and \code{weight}
#' @param alpha Numeric. The alpha value for calculating confidence intervals.
#' @returns A wide-format data frame with the following values: number of observations, mean, standard deviation, weighted mean, and weighted standard deviation.
continuous_analysis <- function(data,
                                alpha) {
  # Unweighted values to start
  n <- nrow(data)
  mean <- mean(data$value)
  sd <- sd(data$value)
  variance <- var(data$value)
  bounds <- ci_mean(mean = mean,
                    sd = sd,
                    n = n,
                    alpha = alpha)
  
  # Now the weighted stuff
  # Get the weight-adjusted values
  data$weighted_value <- data$value * data$weight / sum(data$weight)
  # Weighted mean is the sum of the weight-adjusted values divided by the sum of all weights
  mean_weighted <- sum(data$weighted_value)
  # Standard deviation is calculated differently for weighted values than unweighted
  sd_weighted <- sqrt(sum(data$weight * (data$value - mean)^2) / ((n - 1) / n * sum(data$weight)))
  # So is variance
  variance_weighted <- weighted_variance(values = data$value,
                                         weights = data$weight,
                                         na_remove = FALSE)
  bounds_weighted <- ci_mean(mean = mean_weighted,
                             sd = sd_weighted,
                             n = n,
                             alpha = alpha)
  
  
  
  output <- data.frame(n = n,
                       alpha = alpha,
                       mean = mean,
                       sd = sd,
                       variance = variance,
                       mean_bound_lower = bounds$lower_bound,
                       mean_bound_upper = bounds$upper_bound,
                       mean_weighted = mean_weighted,
                       sd_weighted = sd_weighted,
                       variance_weighted = variance_weighted,
                       mean_weighted_bound_lower = bounds_weighted$lower_bound,
                       mean_weighted_bound_upper = bounds_weighted$upper_bound)
  
  return(output)
}

#' Calculate a weighted variance
#' @param values Numeric vector. The values to calculate the weighted variance for.
#' @param weights Numeric vector. The weights for the vector \code{values}. They must be in the same order as \code{values}.
#' @param na_remove Logical. If \code{TRUE} then any data with either a value or weight of \code{NA} will be removed before calculating. Defaults to \code{FALSE}.
weighted_variance <- function(values,
                              weights,
                              na_remove = FALSE) {
  # Remove the NAs if asked
  if (na_remove) {
    valid_indices <- !is.na(values) & !is.na(weights)
    values <- values[valid_indices]
    weights <- weights[valid_indices]
  }
  # Get the sum of the weights
  sum_of_weights <- sum(weights)
  # Get the sum of the squares of the weights
  sum_of_weights_squared <- sum(weights^2)
  # Get the weighted mean
  weighted_mean <- sum(values * weights) / sum_of_weights
  # Calculate variance!
  variance <- (sum_of_weights / (sum_of_weights^2 - sum_of_weights_squared)) * sum(weights * (values - weighted_mean)^2,
                                                                                   na.rm = na_remove)
  return(variance)
}

#' Check whether results are within a certain tolerance of the true value
#' @param data Data frame. Must contain the variables \code{variable} and \code{comparison_variable}.
#' @param variable Character string. The name of the variable in \code{data} containing the values to check against the tolerance.
#' @param comparison_variable Character string. The name of the variable in \code{data} containing the values to to use to calculate the tolerance.
#' @param percent_tolerance Numeric. The percent difference from the value in \code{comparison_variable} that its paired value in \code{variable} is allowed to be.
#' @returns Logical vector. The value is \code{TRUE} for any index where the value in \code{variable} was within the permitted tolerance and \code{FALSE} for all other indices.
tolerance_test <- function(data,
                           variable,
                           comparison_variable,
                           percent_tolerance = 5) {
  if (class(data) != "data.frame") {
    stop("`data` must be a data frame")
  }
  if (nrow(data) < 1) {
    stop("`data` must contain at least one row of values")
  }
  if (!(variable %in% names(data))) {
    stop(paste0("The variable ", variable, " is missing from data"))
  }
  if (!(comparison_variable %in% names(data))) {
    stop(paste0("The variable ", comparison_variable, " is missing from data"))
  }
  if (percent_tolerance < 0 | percent_tolerance > 100) {
    stop("'percent_tolerance' must be a value between 0 and 100")
  }
  
  proportion_tolerance <- percent_tolerance / 100
  magnitude_tolerance <- abs(proportion_tolerance * data[[comparison_variable]])
  magnitude_difference <- abs(data[[variable]] - data[[comparison_variable]])
  magnitude_difference < magnitude_tolerance
}


#' Summarize the results of tolerance testing
#' @param data Data frame. Must contain the variables \code{variable} and \code{comparison_variable}.
#' @param variable Character string. The name of the variable in \code{data} containing the values to check against the tolerance.
#' @param comparison_variable Character string. The name of the variable in \code{data} containing the values to to use to calculate the tolerance.
#' @param grouping_variables Vector of character strings. The names of variables to group the summary by. If \code{NULL} then all data will be summarized together. Defaults to \code{NULL}.
#' @param percent_tolerance Numeric. The percent difference from the value in \code{comparison_variable} that its paired value in \code{variable} is allowed to be.
tolerance_summary <- function(data,
                              variable,
                              comparison_variable,
                              grouping_variables = NULL,
                              percent_tolerance = 5) {
  if (class(data) != "data.frame") {
    stop("`data` must be a data frame")
  }
  if (nrow(data) < 1) {
    stop("`data` must contain at least one row of values")
  }
  if (!(variable %in% names(data))) {
    stop(paste0("The variable ", variable, " is missing from data"))
  }
  if (!(comparison_variable %in% names(data))) {
    stop(paste0("The variable ", comparison_variable, " is missing from data"))
  }
  if (percent_tolerance < 0 | percent_tolerance > 100) {
    stop("`percent_tolerance` must be a value between 0 and 100")
  }
  missing_grouping_variables <- grouping_variables[!(grouping_variables %in% names(data))]
  if (length(missing_grouping_variables) > 0) {
    stop(paste0("The following variables are missing from data: ", paste(missing_grouping_variables, collapse = ", ")))
  }
  
  within_tolerance <- tolerance_test(data = data,
                                     variable = variable,
                                     comparison_variable = comparison_variable,
                                     percent_tolerance = percent_tolerance)
  threshold_var_name <- paste0("within_", percent_tolerance, "_percent")
  data[[threshold_var_name]] <- within_tolerance
  
  
  if (is.null(grouping_variables)) {
    # How many sims are we looking at?
    n_observations <- nrow(data)
    # Which variable has the info about whether the tolerance was met or not?
    var_name <- paste0("within_", percent_tolerance, "_percent")
    # How many sims were within the tolerance?
    within_tolerance_count <- sum(data[[var_name]])
    # What's the proportion within the tolerance?
    proportion_within_tolerance <- within_tolerance_count / n_observations
    # Gimme those results!
    output <- data.frame(n_observations = n_observations,
                         percent_tolerance = percent_tolerance,
                         n_within_tolerance = within_tolerance_count,
                         proportion_within_tolerance = proportion_within_tolerance)
  } else {
    split_list <- lapply(X = grouping_variables,
                         data = data,
                         FUN = function(X, data) {
                           data[[X]]
                         })
    data_list <- split(x = data,
                       f = split_list)
    
    tolerance_list <- lapply(X = data_list,
                             percent_tolerance = percent_tolerance,
                             grouping_variables = grouping_variables,
                             FUN = function(X, percent_tolerance, grouping_variables) {
                               # How many sims are we looking at?
                               n_observations <- nrow(X)
                               # Which variable has the info about whether the tolerance was met or not?
                               var_name <- paste0("within_", percent_tolerance, "_percent")
                               # How many sims were within the tolerance?
                               within_tolerance_count <- sum(X[[var_name]])
                               # What's the proportion within the tolerance?
                               proportion_within_tolerance <- within_tolerance_count / n_observations
                               # Gimme those results!
                               output <- data.frame(n_observations = n_observations,
                                                    percent_tolerance = percent_tolerance,
                                                    n_within_tolerance = within_tolerance_count,
                                                    proportion_within_tolerance = proportion_within_tolerance)
                               
                               for (variable in grouping_variables) {
                                 output[[variable]] <- X[[variable]][1]
                               }
                               
                               output[, c(grouping_variables, "n_observations", "percent_tolerance", "n_within_tolerance", "proportion_within_tolerance")]
                             })
    output <- do.call(rbind,
                      tolerance_list)
  }
  return(output)
}

#' Read in results from simulations
#' @description This assumes that the simulations are all kept in the same folder and contain the folder "output/results" with only the tabular results as CSV files. It also assumes that the results files have a unique identifying number at the end of their filenames.
#' @param simulations Vector of character strings. The names of the simulations (that is, the folder names) to read in.
#' @param simulations_path Character string. The filepath to where the folders matching the name(s) in \code{simulations} are stored.
#' @param type Character string. Which set of results to read in, \code{"analysis"} or \code{"points"}. Defaults to \code{"analysis"}.
#' @returns The combined results of all the simulations with additional fields identifying which simulation run and which simulation the results belong to.
read_results <- function(simulations,
                         simulations_path,
                         type = "analysis") {
  if (!is.character(simulations)) {
    stop("`simulations` must be a character vector of at least length one")
  }
  if (!is.character(simulations_path)) {
    stop("`simulations_path` must be a character string")
  }
  if (length(simulations_path) > 1) {
    stop("`simulations_path` must be a character string")
  }
  if (!(type %in% c("analysis", "points"))) {
    stop("'type' must be either 'analysis' or 'points'")
  }
  
  # Read in the results from all the simulations!
  results_list <- lapply(X = simulations,
                         simulations_path = simulations_path,
                         type = type,
                         FUN = function(X, simulations_path, type) {
                           # Just so this is more readable
                           sim_name <- X
                           
                           # Where are the results for this sim living?
                           results_path <- paste0(simulations_path, "/",
                                                  sim_name,
                                                  "/output/results")
                           
                           if (!file.exists(results_path)) {
                             stop("The following filepath does not exist: ", results_path)
                           }
                           
                           # Get a list of all the results files in that location
                           # Find the results files
                           results_files <- switch(type,
                                                   "analysis" = {
                                                     list.files(path = results_path,
                                                                pattern = "results[_]?\\d*\\.(csv|CSV)$",
                                                                full.names = TRUE)
                                                   },
                                                   "points" = {
                                                     list.files(path = results_path,
                                                                pattern = "points[_]?\\d*\\.(csv|CSV)$",
                                                                full.names = TRUE)
                                                   })
                           
                           # Read in the results
                           results_list <- lapply(X = results_files,
                                                  FUN = function(X) {
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
                                                    
                                                    # Just to correct in case the otuput is old-style
                                                    names(current_results)[names(current_results) == "weighted"] <- "weighting"
                                                    
                                                    current_results
                                                  })
                           
                           # Combine all the results!
                           full_results <- do.call(rbind,
                                                   results_list)
                           
                           # Add in the simulation info
                           full_results$sim_name <- sim_name
                           
                           full_results
                         })
  
  results <- do.call(rbind,
                     results_list)
  
  return(results)
}

#' Summarize Wilcoxon results by p value thresholds
#' @param results Data frame. The results of Wilcoxon ranked sign tests, one result per row, with variables for the weighting scheme (\code{"weighting"}) and the p value (\code{"p_value"}).
#' @param thresholds Numeric vector. The thresholds to compare the p values against. All values must be between 0 and 1. Defaults to \code{c(0.2, 0.1, 0.05, 0.01)}.
summarize_wilcoxon <- function(results,
                               thresholds = c(0.2, 0.1, 0.05, 0.01)) {
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


#' Generate polygons partitioning point densities
#' @param points Sf points object.
#' @param frame Sf polygons object.
#' @param quantiles
#' @param projection
#' @param verbose
density_polygon_gen <- function(points,
                                frame,
                                quantiles = seq(from = 0,
                                                to = 1,
                                                by = 1 / 3),
                                projection = NULL,
                                # round_quantile_names = 2,
                                verbose = FALSE){
  if (class(quantiles) != "numeric") {
    stop("quantiles must be numeric: either a single numeric value or a vector of numeric values.")
  }
  
  if (length(quantiles) == 1) {
    if (verbose) {
      message(paste0("Finding probability breaks for ", quantiles, " quantiles."))
    }
    quantiles <- seq(from = 0,
                     to = 1,
                     by = 1 / quantiles)
  }
  
  # Handle reprojecting if necessary
  if (is.null(projection)) {
    projection <- sf::st_crs(points)
  }
  if (!identical(sf::st_crs(points), projection)) {
    points <- sf::st_transform(x = points,
                               crs = projection)
  }
  if (!identical(sf::st_crs(frame), projection)) {
    frame <- sf::st_transform(x = frame,
                              crs = projection)
  }
  
  # Get the point coordinates to feed into spatstat.geom::ppp()
  point_coords <- sf::st_coordinates(points)
  
  # Get the owin object for spatstat.geom::ppp()
  frame_owin <- spatstat.geom::as.owin(frame)
  
  if (verbose) {
    message("Creating point pattern object")
  }
  # Make the ppp object (a point pattern)
  points_ppp <- spatstat.geom::ppp(x = point_coords[, 1],
                                   y = point_coords[, 2],
                                   window = frame_owin)
  
  if (verbose) {
    message("Finding distribution density from point pattern")
  }
  # Get the density info from the point pattern
  points_density <- density(points_ppp)
  
  # Make a data frame of coordinates with the density at each coordinate
  density_df <- expand.grid(y = points_density$yrow,
                            x = points_density$xcol)
  density_df$density <- as.vector(points_density$v)
  
  # Figure out where the quantile breaks are
  quantile_cutoffs <- quantile(x = density_df$density,
                               probs = quantiles,
                               na.rm = TRUE)
  
  if (verbose) {
    message("Classifying area by densities")
  }
  # We're going to make a vector of each quantile "name"
  # e.g., 50%-75%
  # We're also going to assign that to each relevant value
  quantiles <- c()
  for (quantile_id in length(quantile_cutoffs):2) {
    # Get the upper and lower cutoff values
    upper <- quantile_cutoffs[quantile_id]
    lower <- quantile_cutoffs[quantile_id - 1]
    
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
    density_df$quantile_id[applicable_indices] <- quantile_id - 1
    
    # Make a name for the quantile, e.g., "50%-75%"
    quantile <- paste0(names(quantile_cutoffs)[quantile_id - 1],
                       "-",
                       names(quantile_cutoffs)[quantile_id])
    
    # Write that name into the relevant indices
    density_df$quantile[applicable_indices] <- quantile
    
    # Save the name
    quantiles <- c(quantiles,
                   quantile)
  }
  
  # We went from highest value to lowest, but we'd like that reversed
  # So here' we'll reverse that
  quantiles <- quantiles[length(quantiles):1]
  
  # Make a lookup table because our conversion to an sf object will only work
  # with the numeric ID and we'd like the quantile names attached
  quantiles_lut <- data.frame(quantile_id = 2:length(quantile_cutoffs) - 1,
                              quantile = factor(quantiles))
  
  # Get a stars object, which is basically a raster
  density_stars <- stars::st_as_stars(.x = density_df,
                                      # This defaults to 1:2, but that transposes
                                      # the x and y axes, so we'll do 2:1
                                      coords = 2:1)
  
  # Convert the stars object to polygons
  # We have to use quantile_id instead of quantile because this'll only work with
  # numeric values. We can always get quantiles in there later with a join
  density_sf <- sf::st_as_sf(density_stars["quantile_id"],
                             as_points = FALSE,
                             merge = TRUE)
  # Make sure that the CRS is assigned
  sf::st_crs(density_sf) <- sf::st_crs(frame)
  
  # Add in the quantile names
  density_sf <- dplyr::left_join(x = density_sf,
                                 y = quantiles_lut)
  
  density_sf
}

#' Generate polygons partitioning point densities
#' @param points Sf points object.
#' @param frame Sf polygons object.
#' @param n_partitions
#' @param point_check
#' @param projection
#' @param verbose
density_polygon_gen_clustered <- function(points,
                                          frame,
                                          n_partitions = 3,
                                          # If it's important that every density polygon contain
                                          # points, this'll try again with a bonus partition
                                          # then merge the empty polygon into the next partition
                                          point_check = TRUE,
                                          density_sample_spacing = 100,
                                          projection = NULL,
                                          verbose = FALSE){
  if (!(class(n_partitions) %in% c("numeric", "integer"))) {
    stop("n_partitions must be numeric: either a single numeric value or a vector of numeric values.")
  }
  
  # Handle reprojecting if necessary
  if (is.null(projection)) {
    projection <- sf::st_crs(points)
  }
  if (!identical(sf::st_crs(points), projection)) {
    points <- sf::st_transform(x = points,
                               crs = projection)
  }
  if (!identical(sf::st_crs(frame), projection)) {
    frame <- sf::st_transform(x = frame,
                              crs = projection)
  }
  
  # Get the point coordinates to feed into spatstat.geom::ppp()
  point_coords <- sf::st_coordinates(points)
  
  # Get the owin object for spatstat.geom::ppp()
  frame_owin <- spatstat.geom::as.owin(frame)
  
  if (verbose) {
    message("Creating point pattern object")
  }
  # Make the ppp object (a point pattern)
  points_ppp <- spatstat.geom::ppp(x = point_coords[, 1],
                                   y = point_coords[, 2],
                                   window = frame_owin)
  
  if (verbose) {
    message("Finding distribution density from point pattern")
  }
  # Get the density info from the point pattern
  points_density <- density(points_ppp,
                            eps = density_sample_spacing)
  
  # Make a data frame of coordinates with the density at each coordinate
  density_df <- expand.grid(y = points_density$yrow,
                            x = points_density$xcol)
  density_df$density <- as.vector(points_density$v)
  
  # Figure out where the breaks are for the partitions
  partition_breaks <- BAMMtools::getJenksBreaks(var = density_df[["density"]],
                                                # This'll a number of values equal to k
                                                # and the terminal values will be the min and max
                                                # So in order to get partition ranges, we need to add
                                                # 1 to k so we get enough breakpoints
                                                k = n_partitions + 1)
  
  if (verbose) {
    message("Classifying area by densities")
  }
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
  if (point_check) {
    message("Carrying out point check")
    partition_ids <- sf::st_drop_geometry(density_sf)$partition_id
    partitions_with_points <- unique(sf::st_drop_geometry(sf::st_intersection(x = points,
                                                                              y = density_sf[, "partition_id"]))[["partition_id"]])
    
    # And while not all partition IDs are represented, try again
    while (!all(partition_ids %in% partitions_with_points)) {
      warning("Not all density partitions contained points. Attempting to with an additional partition with empty partitions being merged up.")
      n_partitions <- n_partitions + 1
      
      # Figure out where the breaks are for the partitions
      partition_breaks <- BAMMtools::getJenksBreaks(var = density_df[["density"]],
                                                    # This'll a number of values equal to k
                                                    # and the terminal values will be the min and max
                                                    # So in order to get partition ranges, we need to add
                                                    # 2 to k so we get enough breakpoints
                                                    k = n_partitions + 2)
      
      if (verbose) {
        message("Classifying area by densities")
      }
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
      density_partition_list <- lapply(X = 1:length(partition_ids),
                                       current_ids = partition_ids[order(partition_ids)],
                                       density_polygons = density_sf,
                                       FUN = function(X, current_ids, density_polygons){
                                         current_polygons <- density_polygons[density_polygons$partition_id == current_ids[X], ]
                                         
                                         current_polygons <- sf::st_as_sf(x = sf::st_union(x = current_polygons))
                                         
                                         current_polygons[["partition_id"]] <- X
                                         
                                         current_polygons
                                       })
      
      density_sf <- do.call(rbind,
                            density_partition_list)
      
      # And now for the while() loop's benefit
      partition_ids <- sf::st_drop_geometry(density_sf)[["partition_id"]]
      partitions_with_points <- unique(sf::st_drop_geometry(sf::st_intersection(x = points,
                                                                                y = density_sf[, "partition_id"]))[["partition_id"]])
      
    }
  }
  
  density_sf
}

#' Identify and repair invalid geometry in spatial polygons data frames
#' @description Using functions from sf, check the geometry of a set of polygons. If the geometry invalid, it attempts to buffer the polygons with \code{sf::st_buffer(dist = 0)}. If the geometry is corrupt or fine, it does nothing.
#' @param polygons Spatial polygons (either sf or spatial polygons data frame). The polygons to be checked. Note that the function will halt if the geometry is corrupt and not return a value.
#' @param verbose Logical. If \code{TRUE} then the function will produce informative messages as it executes its steps. Useful for debugging. Defaults to \code{FALSE}.
#' @param force Logical. If \code{TRUE} then both valid and invalid polygons will be buffered by 0. This shouldn't be necessary, but is a feature for the paranoid.
#' @param use_spherical_geometry Logical. USE AT YOUR OWN RISK. Controls if sf uses spherical geometry or not. For particularly wonky sf objects, this may be necessary to effect any kind of repair but can have unintended consequences for the resulting geometry. If \code{TRUE} then \code{sf::sf_use_s2()} will be set to \code{TRUE} which is the default for the package. Defaults to \code{TRUE}.
#' @return The spatial polygons data frame \code{polygons1}. This will be unchanged if the geometry was valid or repaired if it was invalid.
#' @export
repair_geometry <- function(polygons,
                            verbose = FALSE,
                            force = FALSE,
                            use_spherical_geometry = TRUE) {
  if("SpatialPolygonsDataFrame" %in% class(polygons)) {
    polygons_sf <- sf::st_as_sf(polygons)
    spdf <- TRUE
  } else if ("sf" %in% class(polygons)) {
    polygons_sf <- polygons
    spdf <- FALSE
  } else {
    stop("polygons must either be a spatial polygons data frame or an sf object")
  }
  
  # Get the current setting for using s2 before setting it to what the user
  # has requested for this
  current_s2_status <- sf::sf_use_s2()
  sf::sf_use_s2(use_spherical_geometry)
  
  validity_check <- sf::st_is_valid(polygons_sf)
  
  if (any(is.na(validity_check))) {
    stop("The geometry of the polygons is corrupt. Unable to repair.")
  }
  
  if (!all(validity_check)) {
    if (verbose) {
      message("Invalid geometry found. Attempting to make the geometry valid.")
    }
    output <- sf::st_make_valid(polygons_sf)
  } else if (force) {
    if (verbose) {
      message("No invalid geometry found. Attempting to make the geometry valid anyway.")
    }
    output <- sf::st_make_valid(polygons_sf)
  } else {
    if (verbose) {
      message("No invalid geometry found.")
    }
    output <- polygons_sf
  }
  
  validity_check <- sf::st_is_valid(output)
  
  if (any(is.na(validity_check))) {
    stop("Attempting to make the polygons valid resulted in corrupted geometry. Unable to repair.")
  }
  
  if (!all(validity_check)) {
    if (verbose) {
      message("Invalid geometry still found. Attempting to repair via buffering.")
    }
    output <- sf::st_buffer(x = output,
                            dist = 0)
    
    postbuffer_validity_check <- sf::st_is_valid(output)
    
    if (!all(postbuffer_validity_check)) {
      warning("Even after buffering, geometry appears to still be invalid.")
    }
    
  } else if (force) {
    if (verbose) {
      message("No invalid geometry found. Attempting to repair via buffering anyway.")
    }
    output <- sf::st_buffer(x = output,
                            dist = 0)
  } else {
    if (verbose) {
      message("No invalid geometry found after attempted repair.")
    }
  }
  
  if (spdf) {
    output <- methods::as(output, "Spatial")
  }
  
  # Make sure we don't leave this in an unexpected state
  sf::sf_use_s2(current_s2_status)
  
  return(output)
}

get_ldc_token <- function(username,
                          password) {
  if (is.character(username)) {
    if (length(username) > 1) {
      stop("Your username must be a single character string.")
    }
  } else {
    stop("Your username must be a single character string.")
  }
  if (is.character(password)) {
    if (length(password) > 1) {
      stop("Your password must be a single character string.")
    }
  } else {
    stop("Your password must be a single character string.")
  }
  
  # Attempt to get an authentication response
  authentication_response <- httr::POST(url = "https://oox5sjuicqhezohcpnbsesp32y0yrcbm.lambda-url.us-east-1.on.aws/",
                                        body = list(username = username, 
                                                    password = password),
                                        encode = "json")
  
  # What if there's an error????
  if (httr::http_error(authentication_response)) {
    stop(paste0("Retrieving authentication token from the API failed with status ",
                authentication_response$status_code))
  }
  
  output_raw_character <- rawToChar(authentication_response[["content"]])
  
  if (grepl(x = output_raw_character, pattern = "^Error")) {
    stop(output_raw_character)
  }
  
  output <- jsonlite::fromJSON(txt = output_raw_character)[["AuthenticationResult"]]
  
  # We'll add an expiration time so we can check the need for a refreshed token
  # without making an API call that gets rejected.
  # This cuts 5 seconds off just as a bit of buffer.
  output[["expiration_time"]] <- Sys.time() + output[["ExpiresIn"]] - 5
  # Turns out this was overengineered, but keeping it for future reference.
  # output[["expiration_time"]] <- lubridate::as_date(lubridate::seconds(Sys.time()) + lubridate::seconds(output[["ExpiresIn"]] - 5))
  
  output
}

#' Fetching data from the Landscape Data Commons via API query
#' @description A function for making API calls to the Landscape Data Commons based on the table, key variable, and key variable values. It will return a table of records of the requested data type from the LDC in which the variable \code{key_type} contains only values found in \code{keys}. See the \href{https://api.landscapedatacommons.org/api-docs}{API documentation} to see which variables (i.e. \code{key_type} values) are valid for each data type.
#' @param keys Optional character vector. A character vector of all the values to search for in \code{key_type}. The returned data will consist only of records where \code{key_type} contained one of the key values, but there may be keys that return no records. If \code{NULL} then the entire table will be returned. Defaults to \code{NULL}.
#' @param key_type Optional character string. The name of the variable in the data to search for the values in \code{keys}. This must be the name of a variable that exists in the requested data type's table, e.g. \code{"PrimaryKey"} exists in all tables, but \code{"EcologicalSiteID"} is found only in some. If the function returns a status code of 500 as an error, this variable may not be found in the requested data type. If \code{NULL} then the entire table will be returned. Defaults to \code{NULL}.
#' @param data_type Character string. The type of data to query. Note that the variable specified as \code{key_type} must appear in the table corresponding to \code{data_type}. Valid values are: \code{'gap'}, \code{'header'}, \code{'height'}, \code{'lpi'}, \code{'soilstability'}, \code{'speciesinventory'}, \code{'indicators'}, \code{'species'}, \code{'dustdeposition'}, \code{'horizontalflux'}, and \code{'schema'}.
#' @param username Optional character string. The username to supply to the Landscape Data Commons API. Some data in the Landscape Data Commons are accessible only to users with appropriate credentials. You do not need to supply credentials, but an API request made without them may return fewer or no data. This argument will be ignored if \code{password} is \code{NULL}. Defaults to \code{NULL}.
#' @param password Optional character string. The password to supply to the Landscape Data Commons API.  Some data in the Landscape Data Commons are accessible only to users with appropriate credentials. You do not need to supply credentials, but an API request made without them may return fewer or no data. This argument will be ignored if \code{username} is \code{NULL}. Defaults to \code{NULL}.
#' @param key_chunk_size Numeric. The number of keys to send in a single query. Very long queries fail, so the keys may be chunked into smaller queries with the results of all the queries being combined into a single output. Defaults to \code{100}.
#' @param timeout Numeric. The number of seconds to wait for a nonresponse from the API before considering the query to have failed. Defaults to \code{300}.
#' @param take Optional numeric. The number of records to retrieve at a time. This is NOT the total number of records that will be retrieved! Queries that retrieve too many records at once can fail, so this allows the process to retrieve them in smaller chunks. The function will keep requesting records in chunks equal to this number until all matching records have been retrieved. If this value is too large (i.e., much greater than about \code{10000}), the server will likely respond with a 500 error. If \code{NULL} then all records will be retrieved in a single pass. Defaults to \code{NULL}.
#' @param delay Optional numeric. The number of milliseconds to wait between API queries. Querying too quickly can crash an API or get you locked out, so adjust this as needed. Defaults to \code{500}.
#' @param exact_match Logical. If \code{TRUE} then only records for which the provided keys are an exact match will be returned. If \code{FALSE} then records containing (but not necessarily matching exactly) the first provided key value will be returned e.g. searching with \code{exact_match = FALSE}, \code{keys = "42"}, and \code{key_type = "EcologicalSiteID"} would return all records in which the ecological site ID contained the string \code{"42"} such as \code{"R042XB012NM"} or \code{"R036XB042NM"}. If \code{FALSE} only the first provided key value will be considered. Using non-exact matching will dramatically increase server response times, so use with caution. Defaults to \code{TRUE}.
#' @param verbose Logical. If \code{TRUE} then the function will report additional diagnostic messages as it executes. Defaults to \code{FALSE}.
#' @returns A data frame of records from the requested \code{data_type} which contain the values from \code{keys} in the variable \code{key_type}.
#' @seealso
#' * To query for data by spatial location, use \code{\link[=fetch_ldc_spatial]{fetch_ldc_spatial()}}.
#' * To retrieve data by ecological site ID from a table that doesn't include ecological site ID use \code{\link[=fetch_ldc_ecosite]{fetch_ldc_ecosite()}}.
#' @examples
#' # To retrieve all sampling location metadata collected in the ecological sites R036XB006NM and R036XB007NM
#' headers <- fetch_ldc(keys = c("R036XB006NM", "R036XB007NM"), key_type = "EcologicalSiteID", data_type = "header")
#' # To retrieve all LPI data collected in ecological sites in the 036X Major Land Resource Area (MLRA)
#' relevant_headers <- fetch_ldc(keys = "036X", key_type = "EcologicalSiteID", data_type = "header", exact_match = FALSE)
#' lpi_data <- fetch_ldc(keys = relevant_headers$PrimaryKey, key_type = "PrimaryKey". data_type = "lpi", take = 10000)
#' @export
fetch_ldc <- function(keys = NULL,
                      key_type = NULL,
                      data_type,
                      username = NULL,
                      password = NULL,
                      token = NULL,
                      key_chunk_size = 100,
                      timeout = 300,
                      take = NULL,
                      delay = 500,
                      exact_match = TRUE,
                      verbose = FALSE) {
  user_agent <- "http://github.com/Landscape-Data-Commons/trex"
  # base_url <- "https://api.landscapedatacommons.org/api/v1/"
  base_url <- "https://devapi.landscapedatacommons.org/api/v1/"
  valid_tables <- data.frame(data_type = c("gap",
                                           "header",
                                           "height",
                                           "lpi",
                                           "soilstability",
                                           "speciesinventory",
                                           "indicators",
                                           "species",
                                           "dustdeposition",
                                           "horizontalflux",
                                           "schema",
                                           "dataGap",
                                           "dataHeader",
                                           "dataHeight",
                                           "dataLPI",
                                           "dataSoilStability",
                                           "dataSpeciesInventory",
                                           "geoIndicators",
                                           "geoSpecies",
                                           "dataDustDeposition",
                                           "dataHorizontalFlux",
                                           "aerosummary"),
                             table_name = c("dataGap",
                                            "dataHeader",
                                            "dataHeight",
                                            "dataLPI",
                                            "dataSoilStability",
                                            "dataSpeciesInventory",
                                            "geoIndicators",
                                            "geoSpecies",
                                            "dataDustDeposition",
                                            "dataHorizontalFlux",
                                            "tbl-schema/latest",
                                            "dataGap",
                                            "dataHeader",
                                            "dataHeight",
                                            "dataLPI",
                                            "dataSoilStability",
                                            "dataSpeciesInventory",
                                            "geoIndicators",
                                            "geoSpecies",
                                            "dataDustDeposition",
                                            "dataHorizontalFlux",
                                            "aerosummary"))
  if (!(data_type %in% valid_tables$data_type)) {
    stop(paste0("data_type must be one of the following character strings (some are aliases of each other): ",
                paste(valid_tables$data_type,
                      collapse = ", "),
                "."))
  }
  
  current_table <- valid_tables[["table_name"]][valid_tables$data_type == data_type]
  
  if (!(class(keys) %in% c("character", "NULL"))) {
    stop("keys must be a character string or vector of character strings or NULL.")
  }
  
  if (!(class(key_type) %in% c("character", "NULL"))) {
    stop("key_type must be a character string or NULL.")
  }
  
  if (!is.null(keys) & is.null(key_type)) {
    stop("Must provide key_type when providing keys.")
  }
  
  if (!is.null(take)) {
    if (!is.numeric(take) | length(take) > 1) {
      stop("take must either be NULL or a single numeric value.")
    }
  }
  
  if (delay < 0) {
    stop("delay must be a positive numeric value.")
  } else {
    # Convert the value from milliseconds to nanoseconds because we'll be using
    # microbenchmark::get_nanotime() which returns the current time in nanoseconds
    delay <- delay * 10^6
  }
  
  # Check the user credentials
  if (!is.null(token)) {
    if (class(token) == "list") {
      if (!("IdToken" %in% names(token))) {
        stop("A valid bearer ID token must be a single character string or a single character string in a list stored at an index named 'IdToken'.")
      } else if (class(token[["IdToken"]]) != "character") {
        stop("A valid bearer ID token must be a single character string or a single character string in a list stored at an index named 'IdToken'.")
      }
      if (!("expiration_time" %in% names(token))) {
        token[["expiration_time"]] <- Sys.time() + 3600
      }
    } else if (class(token) == "character") {
      if (length(token) != 1) {
        stop("A valid bearer ID token must be a single character string or a single character string in a list stored at an index named 'IdToken'.")
      }
      token <- list(IdToken = token,
                    # This is just a rough guess on the validity window of the
                    # token so we don't have to do more complicated handling
                    # later.
                    expiration_time = Sys.time() + 3600)
    }
  } else {
    if (!identical(is.null(username), is.null(password))) {
      if (is.null(username)) {
        warning("No token or username provided. Ignoring provided password and retrieving only data which do not require credentials.")
      }
      if (is.null(password)) {
        warning("No token or password provided. Ignoring provided username and retrieving only data which do not require credentials.")
      }
    } else if (!is.null(username) & !is.null(password)) {
      if (class(username) != "character" | length(username) > 1) {
        stop("Provided username must be a single character string.")
      }
      if (class(password) != "character" | length(password) > 1) {
        stop("Provided username must be a single character string.")
      }
      token <- get_ldc_token(username = username,
                             password = password)
    } else if (verbose) {
      message("No credentials provided. Retrieving only data which do not require credentials.")
    }
  }
  
  
  # If there are no keys, grab the whole table
  if (is.null(keys)) {
    if (verbose) {
      message("No keys provided; retrieving all records.")
    }
    if (!is.null(key_type)) {
      warning("No keys provided. Ignoring key_type and retrieving all records.")
    }
    queries <- paste0(base_url,
                      current_table)
  } else {
    # If there are keys, chunk them then build queries
    # This helps prevent queries so long that they fail
    if (verbose) {
      message("Grouping keys into chunks for queries.")
    }
    # We don't know whether the keys came in as a vector of single keys or if
    # one or more of the character strings contains keys separated by commas
    # so we're going to handle that an get a vector of single-key strings
    keys_vector <- unlist(lapply(X = keys,
                                 FUN = function(X) {
                                   trimws(unlist(stringr::str_split(string = X,
                                                                    pattern = ",")))
                                 }))
    # OKAY! So it turns out that it's not impossible for keys to contain
    # ampersands which will result in malformed API queries, so we'll replace
    # them with the unicode reference %26
    keys_vector_original <- keys_vector
    keys_vector <- gsub(x = keys_vector,
                        pattern = "[&]",
                        replacement = "%26")
    
    if (verbose & !identical(keys_vector_original, keys_vector)) {
      warning("Some keys provided contained illegal characters and have been sanitized. All available data should still be retrieved for all provided keys.")
    }
    
    if (!exact_match) {
      if (verbose) {
        message("Using non-exact matching for the key value.")
      }
      if (length(keys_vector) > 1) {
        warning("There are multiple provided key values. Non-exact matching will only consider the first.")
      }
      keys_vector <- keys_vector[1]
    }
    
    # Figure out how many chunks to break these into based on the max number of
    # keys in a chunk
    key_chunk_count <- ceiling(length(keys_vector) / key_chunk_size)
    
    # Make the key chunks
    # For each chunk, figure out the appropriate indices and paste together the
    # relevant key values into strings that we can use to build per-chunk queries
    keys_chunks <- sapply(X = 1:key_chunk_count,
                          keys_vector = keys_vector,
                          key_chunk_size = key_chunk_size,
                          key_count = length(keys_vector),
                          FUN = function(X, keys_vector, key_chunk_size, key_count) {
                            min_index <- max(c(1, (X - 1) * key_chunk_size + 1))
                            max_index <- min(c(key_count, X * key_chunk_size))
                            indices <- min_index:max_index
                            paste(keys_vector[indices],
                                  collapse = ",")
                          })
    
    if (verbose) {
      if (length(keys_chunks == 1)) {
        message("Building query.")
      } else {
        message("Building queries.")
      }
    }
    
    if (exact_match) {
      queries <- paste0(base_url,
                        current_table,
                        "?",
                        key_type,
                        "=",
                        keys_chunks)
    } else {
      # This adds "Like" to the end of the variable name to do a search for a
      # non-exact match. The object is still called "queries" even though it
      # had better be a single string instead of a vector.
      queries <- paste0(base_url,
                        current_table,
                        "?",
                        key_type,
                        "Like=",
                        keys_chunks)
    }
  }
  
  # Use the queries to snag data
  # This produces a list of results where each index in the list contains the
  # results of one query.
  # It uses a loop instead of a lapply() so that we can check if the token has
  # expired each time we use it.
  data_list <- list()
  
  for (current_query in queries) {
    # The token might expire and need refreshing!
    if (!is.null(token)) {
      if (Sys.time() > token[["expiration_time"]]) {
        if (verbose) {
          message("Current API bearer authorization token has expired. Attempting to request a new one.")
        }
        if (!is.null(username) & !is.null(password)) {
          token <- get_ldc_token(username = username,
                                 password = password)
        } else {
          warning("The API bearer authorization token has expired. Because username and password have not been provided, only data which do not require a token will be retrieved.")
          token <- NULL
        }
      }
    }
    
    # We handle things differently if the data type is header
    # because the header table doesn't have an rid variable
    # and we can't use take or cursor options without that
    
    if (data_type == "header" | is.null(take)) {
      if (verbose) {
        message("Attempting to query LDC with:")
        message(current_query)
      }
      
      # Full query response using the token if we've got one.
      if (is.null(token)) {
        response <- httr::GET(url = current_query,
                              httr::timeout(timeout),
                              httr::user_agent(user_agent))
        
      } else {
        response <- httr::GET(url = current_query,
                              httr::timeout(timeout),
                              httr::user_agent(user_agent),
                              httr::add_headers(Authorization = paste("Bearer",
                                                                      token[["IdToken"]])))
        response <- httr::GET(url = current_query,
                              httr::add_headers(Authorization = paste("Bearer",
                                                                      token[["IdToken"]])))
      }
      # if (!is.null(username) & !is.null(password)) {
      #   response <- httr::GET(current_query,
      #                         config = list(httr::timeout(timeout),
      #                                       httr::user_agent(user_agent),
      #                                       httr::authenticate(user = username,
      #                                                          password = password)))
      # } else {
      #   response <- httr::GET(current_query,
      #                         config = list(httr::timeout(timeout),
      #                                       httr::user_agent(user_agent)))
      # }
      
      
      # What if there's an error????
      if (httr::http_error(response)) {
        if (response$status_code == 500) {
          stop(paste0("Query failed with status ",
                      response$status_code,
                      " which may be due to a very large number of records returned or attempting to query using a variable that doesn't occur in the requested data table. Consider setting the take argument to 10000 or less and consult https://api.landscapedatacommons.org/api-docs to see which variables are in which tables."))
        } else {
          stop(paste0("Query failed with status ",
                      response$status_code))
        }
      }
      
      # Grab only the data portion
      response_content <- response[["content"]]
      # Convert from raw to character
      content_character <- rawToChar(response_content)
      # Convert from character to data frame
      content_df <- jsonlite::fromJSON(content_character)
      
    } else {
      # OKAY! So handling using take and cursor options for
      # anything non-header
      # The first query needs to not specify the cursor position
      # and then after that we'll keep trying with the last
      # rid value + 1 as the cursor until we get an empty
      # response
      if (verbose) {
        message(paste0("Retrieving records in chunks of ", take))
      }
      
      query <- paste0(current_query, "&take=", take)
      
      if (verbose) {
        message("Attempting to query LDC with:")
        message(query)
      }
      
      # Querying with the token if we've got it.
      if (is.null(token)) {
        response <- httr::GET(url = query,
                              httr::timeout(timeout),
                              httr::user_agent(user_agent))
        
      } else {
        response <- httr::GET(url = query,
                              httr::timeout(timeout),
                              httr::user_agent(user_agent),
                              httr::add_headers(Authorization = paste("Bearer",
                                                                      token[["IdToken"]])))
      }
      
      # What if there's an error????
      if (httr::http_error(response)) {
        if (response$status_code == 500) {
          stop(paste0("Query failed with status ",
                      response$status_code,
                      " which may be due to a very large number of records returned or attempting to query using a variable that doesn't occur in the requested data table. Consider setting the take argument to 10000 or less and consult https://api.landscapedatacommons.org/api-docs to see which variables are in which tables."))
        } else {
          stop(paste0("Query failed with status ",
                      response$status_code))
        }
      }
      
      # Grab only the data portion
      response_content <- response[["content"]]
      # Convert from raw to character
      content_character <- rawToChar(response_content)
      # Convert from character to data frame
      current_content_df <- jsonlite::fromJSON(content_character)
      
      content_df_list <- list(current_content_df)
      
      # Here's where we start iterating as long as we're still
      # getting data
      # So while the last returned response wasn't empty,
      # keep requesting the next response where the cursor
      # is set to the rid following the the highest rid in
      # the last chunk
      while (length(content_df_list[[length(content_df_list)]]) > 0) {
        # And to avoid flooding the API server with requests,
        # we'll put in a delay here.
        # This gets the current time then spins its wheels,
        # checking repeatedly to see if enough time has
        # elapsed, at which point it moves on
        start_time <- microbenchmark::get_nanotime()
        repeat {
          current_time <- microbenchmark::get_nanotime()
          elapsed_time <- current_time - start_time
          if (elapsed_time > delay) {
            break
          }
        }
        
        
        last_rid <- max(content_df_list[[length(content_df_list)]][["rid"]])
        
        query <- paste0(current_query, "&take=", take, "&cursor=", last_rid)
        
        if (verbose) {
          message("Attempting to query LDC with:")
          message(query)
        }
        
        # The token might expire and need refreshing!
        # This exists up above too, but we need it here in the while loop
        # because we might be in the while for long enough for a token to expire
        if (Sys.time() > token[["expiration_time"]]) {
          if (verbose) {
            message("Current API bearer authorization token has expired. Attempting to request a new one.")
          }
          if (!is.null(username) & !is.null(password)) {
            token <- get_ldc_token(username = username,
                                   password = password)
          }
        }
        
        # Querying with the token if we've got it.
        if (is.null(token)) {
          response <- httr::GET(url = query,
                                httr::timeout(timeout),
                                httr::user_agent(user_agent))
          
        } else {
          response <- httr::GET(url = query,
                                httr::timeout(timeout),
                                httr::user_agent(user_agent),
                                httr::add_headers(Authorization = paste("Bearer",
                                                                        token[["IdToken"]])))
        }
        
        # What if there's an error????
        if (httr::http_error(response)) {
          stop(paste0("Query failed with status ",
                      response$status_code))
        }
        
        # Grab only the data portion
        response_content <- response[["content"]]
        # Convert from raw to character
        content_character <- rawToChar(response_content)
        # Convert from character to data frame
        current_content_df <- jsonlite::fromJSON(content_character)
        
        # Bind that onto the end of the list
        # The data are wrapped in list() so that it gets added
        # as a data frame instead of as a vector for each variable
        content_df_list <- c(content_df_list, list(current_content_df))
      }
      content_df <- do.call(rbind,
                            content_df_list)
      
      # And another delay for between individual queries
      # that were generated by the key chunking instead of
      # by take
      start_time <- microbenchmark::get_nanotime()
      repeat {
        current_time <- microbenchmark::get_nanotime()
        elapsed_time <- current_time - start_time
        if (elapsed_time > delay) {
          break
        }
      }
    }
    # Append whatever it is that we got back
    data_list <- c(data_list, list(content_df))
  }
  
  # data_list <- lapply(X = queries,
  #                     data_type = data_type,
  #                     timeout = timeout,
  #                     take = take,
  #                     user_agent = user_agent,
  #                     token = token,
  #                     username = username,
  #                     password = password,
  #                     verbose = verbose,
  #                     FUN = function(X, data_type, take, timeout, user_agent, token, username, password, verbose){
  #                       
  #                       
  #                       
  #                       # We handle things differently if the data type is header
  #                       # because the header table doesn't have an rid variable
  #                       # and we can't use take or cursor options without that
  #                       
  #                       if (data_type == "header" | is.null(take)) {
  #                         if (verbose) {
  #                           message("Attempting to query LDC with:")
  #                           message(X)
  #                         }
  #                         
  #                         # Full query response
  #                         if (!is.null(username) & !is.null(password)) {
  #                           response <- httr::GET(X,
  #                                                 config = list(httr::timeout(timeout),
  #                                                               httr::user_agent(user_agent),
  #                                                               httr::authenticate(user = username,
  #                                                                                  password = password)))
  #                         } else {
  #                           response <- httr::GET(X,
  #                                                 config = list(httr::timeout(timeout),
  #                                                               httr::user_agent(user_agent)))
  #                         }
  #                         
  #                         
  #                         # What if there's an error????
  #                         if (httr::http_error(response)) {
  #                           if (response$status_code == 500) {
  #                             stop(paste0("Query failed with status ",
  #                                         response$status_code,
  #                                         " which may be due to a very large number of records returned or attempting to query using a variable that doesn't occur in the requested data table. Consider setting the take argument to 10000 or less and consult https://api.landscapedatacommons.org/api-docs to see which variables are in which tables."))
  #                           } else {
  #                             stop(paste0("Query failed with status ",
  #                                         response$status_code))
  #                           }
  #                         }
  #                         
  #                         # Grab only the data portion
  #                         response_content <- response[["content"]]
  #                         # Convert from raw to character
  #                         content_character <- rawToChar(response_content)
  #                         # Convert from character to data frame
  #                         content_df <- jsonlite::fromJSON(content_character)
  #                         
  #                         content_df
  #                       } else {
  #                         # OKAY! So handling using take and cursor options for
  #                         # anything non-header
  #                         # The first query needs to not specify the cursor position
  #                         # and then after that we'll keep trying with the last
  #                         # rid value + 1 as the cursor until we get an empty
  #                         # response
  #                         if (verbose) {
  #                           message(paste0("Retrieving records in chunks of ", take))
  #                         }
  #                         
  #                         query <- paste0(X, "&take=", take)
  #                         
  #                         if (verbose) {
  #                           message("Attempting to query LDC with:")
  #                           message(query)
  #                         }
  #                         
  #                         # Full query response
  #                         response <- httr::GET(query,
  #                                               config = list(httr::timeout(timeout),
  #                                                             httr::user_agent(user_agent)))
  #                         
  #                         # What if there's an error????
  #                         if (httr::http_error(response)) {
  #                           if (response$status_code == 500) {
  #                             stop(paste0("Query failed with status ",
  #                                         response$status_code,
  #                                         " which may be due to a very large number of records returned or attempting to query using a variable that doesn't occur in the requested data table. Consider setting the take argument to 10000 or less and consult https://api.landscapedatacommons.org/api-docs to see which variables are in which tables."))
  #                           } else {
  #                             stop(paste0("Query failed with status ",
  #                                         response$status_code))
  #                           }
  #                         }
  #                         
  #                         # Grab only the data portion
  #                         response_content <- response[["content"]]
  #                         # Convert from raw to character
  #                         content_character <- rawToChar(response_content)
  #                         # Convert from character to data frame
  #                         current_content_df <- jsonlite::fromJSON(content_character)
  #                         
  #                         content_df_list <- list(current_content_df)
  #                         
  #                         # Here's where we start iterating as long as we're still
  #                         # getting data
  #                         # So while the last returned response wasn't empty,
  #                         # keep requesting the next response where the cursor
  #                         # is set to the rid following the the highest rid in
  #                         # the last chunk
  #                         while (length(content_df_list[[length(content_df_list)]]) > 0) {
  #                           # And to avoid flooding the API server with requests,
  #                           # we'll put in a delay here.
  #                           # This gets the current time then spins its wheels,
  #                           # checking repeatedly to see if enough time has
  #                           # elapsed, at which point it moves on
  #                           start_time <- microbenchmark::get_nanotime()
  #                           repeat {
  #                             current_time <- microbenchmark::get_nanotime()
  #                             elapsed_time <- current_time - start_time
  #                             if (elapsed_time > delay) {
  #                               break
  #                             }
  #                           }
  #                           
  #                           
  #                           last_rid <- max(content_df_list[[length(content_df_list)]][["rid"]])
  #                           
  #                           query <- paste0(X, "&take=", take, "&cursor=", last_rid)
  #                           
  #                           if (verbose) {
  #                             message("Attempting to query LDC with:")
  #                             message(query)
  #                           }
  #                           
  #                           # Full query response
  #                           response <- httr::GET(query,
  #                                                 config = list(httr::timeout(timeout),
  #                                                               httr::user_agent(user_agent)))
  #                           
  #                           # What if there's an error????
  #                           if (httr::http_error(response)) {
  #                             stop(paste0("Query failed with status ",
  #                                         response$status_code))
  #                           }
  #                           
  #                           # Grab only the data portion
  #                           response_content <- response[["content"]]
  #                           # Convert from raw to character
  #                           content_character <- rawToChar(response_content)
  #                           # Convert from character to data frame
  #                           current_content_df <- jsonlite::fromJSON(content_character)
  #                           
  #                           # Bind that onto the end of the list
  #                           # The data are wrapped in list() so that it gets added
  #                           # as a data frame instead of as a vector for each variable
  #                           content_df_list <- c(content_df_list, list(current_content_df))
  #                         }
  #                         content_df <- do.call(rbind,
  #                                               content_df_list)
  #                         
  #                         # And another delay for between individual queries
  #                         # that were generated by the key chunking instead of
  #                         # by take
  #                         start_time <- microbenchmark::get_nanotime()
  #                         repeat {
  #                           current_time <- microbenchmark::get_nanotime()
  #                           elapsed_time <- current_time - start_time
  #                           if (elapsed_time > delay) {
  #                             break
  #                           }
  #                         }
  #                         
  #                         content_df
  #                       }
  #                       
  #                     })
  
  # Combine all the results of the queries
  data <- do.call(rbind,
                  data_list)
  
  # If there aren't data, let the user know
  if (length(data) < 1) {
    warning("No data retrieved. Confirm that your keys and key_type are correct.")
    return(NULL)
  } else {
    # If there are data and the user gave keys, find which if any are missing
    if (!is.null(keys) & exact_match) {
      # Note that we're using keys_vector_original because even if we made
      # alterations to keys_vector, the actual retrieved keys should match the
      # original values despite substituting unicode references for illegal characters
      missing_keys <- keys_vector_original[!(keys_vector_original %in% data[[key_type]])]
      if (length(missing_keys) > 0) {
        warning(paste0("The following keys were not associated with data: ",
                       paste(missing_keys,
                             collapse = ",")))
      }
    }
    return(data)
  }
}

#' Fetching data from the Landscape Data Commons using spatial constraints
#' @description A function for retrieving data from the Landscape Data Commons which fall within a given set of polygons. This is accomplished by retrieving the header information for all points in the LDC, spatializing them, and finding the PrimaryKey values associated with points within the given polygons. Those PrimaryKey values are used to retrieve only the qualifying data from the LDC. Every time this function is called, it retrieves ALL header information via the API, which can be slow. If you plan to do multiple spatial queries back-to-back, it'll be faster to retrieve the headers with \code{\link[=fetch_ldc]{fetch_ldc()}} once, convert them to an sf object with \code{sf::st_as_sf()}, then use \code{sf:st_intersection()} repeatedly on that sf object to find the PrimaryKey values for each set of polygons and query the API using the PrimaryKeys.
#' @param polygons Polygon sf object. The polygon or polygons describing the area to retrieve data from. Only records from sampling locations falling within this area will be returned.
#' @param data_type Character string. The type of data to query. Note that the variable specified as \code{key_type} must appear in the table corresponding to \code{data_type}. Valid values are: \code{'gap'}, \code{'header'}, \code{'height'}, \code{'lpi'}, \code{'soilstability'}, \code{'speciesinventory'}, \code{'indicators'}, \code{'species'}, \code{'dustdeposition'}, \code{'horizontalflux'}, and \code{'schema'}.
#' @param username Optional character string. The username to supply to the Landscape Data Commons API. Some data in the Landscape Data Commons are accessible only to users with appropriate credentials. You do not need to supply credentials, but an API request made without them may return fewer or no data. This argument will be ignored if \code{password} is \code{NULL}. Defaults to \code{NULL}.
#' @param password Optional character string. The password to supply to the Landscape Data Commons API.  Some data in the Landscape Data Commons are accessible only to users with appropriate credentials. You do not need to supply credentials, but an API request made without them may return fewer or no data. This argument will be ignored if \code{username} is \code{NULL}. Defaults to \code{NULL}.
#' @param key_chunk_size Numeric. The number of PrimaryKeys to send in a single query. Very long queries fail, so the keys may be chunked into smaller queries with the results of all the queries being combined into a single output. Defaults to \code{100}.
#' @param timeout Numeric. The number of seconds to wait for a nonresponse from the API before considering the query to have failed. Defaults to \code{300}.
#' @param take Optional numeric. The number of records to retrieve at a time. This is NOT the total number of records that will be retrieved! Queries that retrieve too many records at once can fail, so this allows the process to retrieve them in smaller chunks. The function will keep requesting records in chunks equal to this number until all matching records have been retrieved. If this value is too large (i.e., much greater than about \code{10000}), the server will likely respond with a 500 error. If \code{NULL} then all records will be retrieved in a single pass. Defaults to \code{NULL}.
#' @param delay Optional numeric. The number of milliseconds to wait between API queries. Querying too quickly can crash an API or get you locked out, so adjust this as needed. Defaults to \code{500}.
#' @param return_spatial Logical. If \code{TRUE} then the returned data will be an sf object. Otherwise if this is \code{FALSE} it will be a simple data frame. Defaults to \code{TRUE}.
#' @param verbose Logical. If \code{TRUE} then the function will report additional diagnostic messages as it executes. Defaults to \code{FALSE}.
#' @returns A data frame of records from the requested \code{data_type} which came from locations within \code{polygons}.
#' @seealso
#' * To query for data by key values, use \code{\link[=fetch_ldc]{fetch_ldc()}}.
#' * To retrieve data by ecological site ID from a table that doesn't include ecological site ID use \code{\link[=fetch_ldc_ecosite]{fetch_ldc_ecosite()}}.
#' @examples
#' To retrieve all LPI records for sampling locations found within a given set of polygons provided as an sf object
#' fetch_ldc_spatial(polygons = polygons_sf, data_type = "lpi")
#' @export
fetch_ldc_spatial <- function(polygons,
                              data_type,
                              token = NULL,
                              username = NULL,
                              password = NULL,
                              key_chunk_size = 100,
                              timeout = 300,
                              take = NULL,
                              delay =  500,
                              return_spatial = TRUE,
                              verbose = FALSE) {
  if (!("sf" %in% class(polygons))) {
    stop("polygons must be a polygon sf object")
  }
  
  # Just to get a unique ID in there for sure without having to ask the user
  polygons$unique_id <- 1:nrow(polygons)
  
  if (verbose) {
    message("Fetching the header information from the LDC.")
  }
  headers_df <- fetch_ldc(data_type = "header",
                          username = username,
                          password = password)
  
  # We know that the header info includes coordinates in NAD83, so we can easily
  # convert the data frame into an sf object
  if (verbose) {
    message("Converting header information into an sf point object.")
  }
  headers_sf <- sf::st_as_sf(x = headers_df,
                             coords = c("Longitude_NAD83",
                                        "Latitude_NAD83"),
                             crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs +type=crs")
  
  # We're just after the PrimaryKey values here
  if (verbose) {
    message("Finding points that fall within the polygons.")
  }
  header_polygons_intersection <- sf::st_intersection(x = headers_sf[, "PrimaryKey"],
                                                      y = sf::st_transform(polygons[, "unique_id"],
                                                                           crs = sf::st_crs(headers_sf)))
  
  # What if there're no qualifying data????
  if (nrow(header_polygons_intersection) < 1) {
    warning("No data were located within the given polygons.")
    return(NULL)
  }
  
  # If there were points found, snag the PrimaryKey values
  intersected_primarykeys <- unique(header_polygons_intersection$PrimaryKey)
  
  # Grab only the data associated with the PrimaryKey values we've got
  if (data_type == "header") {
    output <- headers_sf[headers_sf$PrimaryKey %in% intersected_primarykeys, ]
    if (!return_spatial) {
      output <- sf::st_drop_geometry(output)
    }
  } else {
    output <- fetch_ldc(keys = intersected_primarykeys,
                        key_type = "PrimaryKey",
                        data_type = data_type,
                        token = token,
                        username = username,
                        password = password,
                        key_chunk_size = key_chunk_size,
                        timeout = timeout,
                        exact_match = TRUE,
                        delay = delay,
                        verbose = verbose)
    if (return_spatial) {
      output <- dplyr::inner_join(x = dplyr::select(.data = headers_sf,
                                                    PrimaryKey),
                                  y = output)
    }
  }
  
  output
}