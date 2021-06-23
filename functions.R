#' Create Thiessen/Voronoi polygons from a set of points and bounding polygons
#' @description Generate Thiessen/Voronoi polygons for a set of points and clip the results using a set of polygons
#' @param centroids An sf points object. These points are used as centroids for the Thiessen/Voronoi polygons.
#' @param frame An sf polygon or multipolygon object. This is the clipping boundary which will be applied to the otherise "infinite" Thiessen/Voronoi polygons.
#' @param use_albers Logical. If \code{TRUE} then \code{centroids} and \code{frame} will be reprojected into Albers Equal Area (AEA) and the output will be in AEA. If \code{FALSE} then \code{frame} will be reprojected to match the coordinate reference ssytem (CRS) of \code{centroids} and the output will be in that CRS. CRSs using decimal degrees will throw errors or warnings. Defaults to \code{TRUE}.
#' @return An sf object composed of polygon or multipolygon geometry
thiessen_polygons_gen_fixed <- function(centroids,
                                  frame,
                                  seed_number = 420,
                                  use_albers = TRUE){
  # Define Alber's Equal Area CRS
  aea_proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  
  # Sanitization
  if (!("sf" %in% class(centroids))) {
    stop("centroids must be an sf points object")
  } else if (!all(sf::st_geometry_type(centroids) %in% c("POINT"))){
    stop("centroids must be an sf points object")
  }
  if (!("sf" %in% class(frame))) {
    stop("frame must be an sf polygon object")
  } else if (!all(sf::st_geometry_type(frame) %in% c("POLYGON", "MULTIPOLYGON"))){
    stop("frame must be an sf polygon object")
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
  } else {
    # This just forces the polygons into the same projection as the centroids
    centroids_crs <- sf::st_crs(centroids)
    frame <- sf::st_transform(frame,
                              crs = centroids_crs)
  }
  
  # Draw Thiessen polygons
  # Here's where it gets weird
  # The points need to be a multipoint object, apparently
  points_multipoint <- sf::st_combine(sf::st_geometry(centroids))
  
  # Generate the Thiessen polygons, not bothering to try to use the sample frame as an envelope
  # This is because it ignores the envelope argument if the envelope polygons are smaller than the default boundaries
  # (and they probably will be)
  # This appears to be a list??? It's fine, I promise. We'll convert it in a bit
  thiessen_polygons_raw <- sf::st_voronoi(points_multipoint)
  
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
  thiessen_polygons_clipped$polygon_unique_id <- 1:nrow(thiessen_polygons_clipped)
  
  # Add in the areas for the polygons
  thiessen_polygons_clipped$area_m2 <- as.vector(sf::st_area(x = thiessen_polygons_clipped))
  
  return(thiessen_polygons_clipped)
}

#' Create Thiessen/Voronoi polygons from a set of points and bounding polygons
#' @description Generate Thiessen/Voronoi polygons for a set of points and clip the results using a set of polygons
#' @param frame An sf polygon or multipolygon object. This is the clipping boundary which will be applied to the otherise "infinite" Thiessen/Voronoi polygons.
#' @param n_polygons Numeric value. The number of Thiessen polygons to draw within the frame.
#' @param points Optional sf point object. If provided, then the Thiessen polygons will be redrawn with new random seeds until each contains at least \code{points_min} of these points. Defaults to \code{NULL}.
#' @param points_min Optional numeric value. If \code{points} is not \code{NULL} then this is the minimum number of points that each Thiessen polygon will contain. Defaults to \code{2}.
#' @param seed_number Optional numeric value. The seed number to use for generating the polygon centroids. A random seed will be used if this is \code{NULL}. Defaults to \code{NULL}.
#' @param seed_increment Optional numeric value. If attempting to produce polygons with \code{points_min} points from \code{points} in each polygon, this is the step to increment \code{seed_number} by on each attempt. Defaults to \code{100000}.
#' @param use_albers Logical. If \code{TRUE} then \code{centroids} and \code{frame} will be reprojected into Albers Equal Area (AEA) and the output will be in AEA. If \code{FALSE} then everything will be reprojected to match the coordinate reference system (CRS) of \code{frame} and the output will be in that CRS. CRSs using decimal degrees will throw errors or warnings. Defaults to \code{TRUE}.
#' @param verbose Logical. If \code{TRUE} then the function will return diagnostic messages as it runs. Defaults to \code(FALSE}.)
#' @return An sf object composed of polygon or multipolygon geometry
thiessen_polygons_gen_random <- function(frame,
                                         n_polygons,
                                         points = NULL,
                                         points_min = 2,
                                         seed_number = NULL,
                                         seed_increment = 100000,
                                         use_albers = TRUE,
                                         verbose = FALSE){
  # Define Alber's Equal Area CRS
  projection <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  
  # Sanitization
  if (!("sf" %in% class(frame))) {
    stop("frame must be an sf polygon object")
  } else if (!all(sf::st_geometry_type(frame) %in% c("POLYGON", "MULTIPOLYGON"))){
    stop("frame must be an sf polygon object")
  }
  if (!(class(n_polygons) %in% c("numeric", "integer")) | length(n_polygons) > 1) {
    stop("n_points must be a single numeric value")
  }
  if (!is.null(points)) {
    if (!("sf" %in% class(points))) {
      stop("points must be an sf points object")
    } else if (!all(sf::st_geometry_type(points) %in% c("POINT"))){
      stop("points must be an sf points object")
    }
  }
  if (!is.null(seed_number)) {
    if (!(class(seed_number) %in% c("numeric", "integer")) | length(seed_number) > 1) {
      stop("seed_number must be a single numeric value")
    }
  } else {
    seed_number <- sample(x = 1:9999999,
                          size = 1)
  }
  
  # Remove any Z dimension
  # It screws with the process and is irrelevant
  frame <- sf::st_zm(frame,
                     drop = TRUE)
  if (!is.null(points)) {
    points <- sf::st_zm(points,
                        drop = TRUE)
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
    # This just forces the polygons into the same projection as the centroids
    projection <- sf::st_crs(frame)
    if (!is.null(points)) {
      points <- sf::st_transform(x = points,
                                 crs = projection)
    }
  }
  
  # Draw centroids
  centroids <- points_gen(frame = frame,
                          sample_type = "simple",
                          n_points = n_polygons,
                          seed_number = seed_number,
                          projection = projection)
  
  # Draw Thiessen polygons
  thiessen_polygons <- thiessen_polygons_gen_fixed(centroids = centroids,
                                                   frame = frame,
                                                   seed_number = seed_number,
                                                   use_albers = FALSE)
  
  # Get the final variables in there
  thiessen_polygons$tpoly_seed <- seed_number
  thiessen_polygons$tpoly_id <- paste0("tpoly_",
                                       thiessen_polygons$tpoly_seed,
                                       "-",
                                       thiessen_polygons$polygon_unique_id)
  
  if (!is.null(points)) {
    ## Check that polygons contain enough points
    points_attributed <- sf::st_join(x = points,
                                     y = thiessen_polygons[, c("tpoly_id")])
    
    tpoly_summary <- data.frame(tpoly_id = names(table(points_attributed$tpoly_id)),
                                n_points = as.vector(table(points_attributed$tpoly_id)),
                                stringsAsFactors = FALSE)
    
    while (!all(tpoly_summary[["n_points"]] >= points_min)) {
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
      
      thiessen_polygons <- thiessen_polygons_gen_fixed(centroids = centroids,
                                                       frame = frame,
                                                       seed_number = seed_number,
                                                       use_albers = FALSE)
      
      
      # Get the final variables in there
      thiessen_polygons$tpoly_seed <- seed_number
      thiessen_polygons$tpoly_id <- paste0("tpoly_",
                                           thiessen_polygons$tpoly_seed,
                                           "-",
                                           thiessen_polygons$polygon_unique_id)
      
      ## Check that polygons contain enough points
      points_attributed <- sf::st_join(x = points,
                                       y = thiessen_polygons[, c("tpoly_id")])
      
      tpoly_summary <- data.frame(tpoly_id = names(table(points_attributed$tpoly_id)),
                                  n_points = as.vector(table(points_attributed$tpoly_id)),
                                  stringsAsFactors = FALSE)
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
    return(thiessen_polygons_clipped[, c("tpoly_id", "tpoly_seed")])
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
                                      projection = NULL){
  if (class(ncol) != "numeric" | length(ncol) > 1) {
    stop("ncol must be a single numeric value")
  }
  if (class(nrow) != "numeric" | length(nrow) > 1) {
    stop("nrow must be a single numeric value")
  }
  
  if (!is.null(seed_number)) {
    if (!(class(seed_number) %in% c("numeric", "integer")) | length(seed_number) > 1) {
      stop("seed_number must be a single numeric value")
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
    stop("projection must either be a valid PROJ4 string or CRS object")
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
                                     projection = NULL){
  if (class(max) != "numeric" | length(max) > 1) {
    stop("max must be a single numeric value")
  }
  if (class(min) != "numeric" | length(min) > 1) {
    stop("min must be a single numeric value")
  }
  if (max < min) {
    stop("max must be greater than min")
  }
  
  if (!(distribution %in% c("normal", "uniform"))) {
    stop("distribution must be either 'normal' or 'uniform'")
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
    stop("ncol must be a single numeric value")
  }
  if (class(nrow) != "numeric" | length(nrow) > 1) {
    stop("nrow must be a single numeric value")
  }
  
  if (!is.null(seed_number)) {
    if (!(class(seed_number) %in% c("numeric", "integer")) | length(seed_number) > 1) {
      stop("seed_number must be a single numeric value")
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
    stop("projection must either be a valid PROJ4 string or CRS object")
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
#' @param projection Optional character string or CRS object. The coordinate reference system for the raster. May be a PROJ4 string or a CRS object. Defaults to Albers Equal Area.
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
    stop("xmax must be a single numeric value")
  }
  if (class(xmin) != "numeric" | length(xmin) > 1) {
    stop("xmin must be a single numeric value")
  }
  if (xmax < xmin) {
    stop("xmax must be greater than xmin")
  }
  if (class(ymax) != "numeric" | length(ymax) > 1) {
    stop("ymax must be a single numeric value")
  }
  if (class(ymin) != "numeric" | length(ymin) > 1) {
    stop("ymin must be a single numeric value")
  }
  if (ymax < ymin) {
    stop("ymax must be greater than ymin")
  }
  
  if (!is.null(seed_number)) {
    if (!(class(seed_number) %in% c("numeric", "integer")) | length(seed_number) > 1) {
      stop("seed_number must be a single numeric value")
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
    stop("projection must either be a valid PROJ4 string or CRS object")
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

#' Make a single-part polygon (more) concave
#' @param polygon An sf polygon object. The must not be multipart or have multiple polygons. It may already be concave.
#' @param seed_number Optional numeric value. The seed number used to generate the coordinates of the new vertex. If \code{NULL} then a random number will be used. Defaults to \code{NULL}.
#' @returns An sf polygon object containing the original polygon, but with a new vertex that adds a concavity to it.
concavify <- function(polygon,
                      seed_number = NULL) {
  if (!("sf" %in% class(polygon))) {
    stop("polygon must be an sf object of geometry type 'POLYGON'")
  } else if (sf::st_geometry_type(polygon) != "POLYGON") {
    stop("polygon must be an sf object of geometry type 'POLYGON'")
  }
  
  if (!is.null(seed_number)) {
    if (!(class(seed_number) %in% c("numeric", "integer")) | length(seed_number) > 1) {
      stop("seed_number must be a single numeric value")
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
                       projection = NULL){
  if (!("sf" %in% class(frame))) {
    stop("frame must be an sf object of geometry type 'POLYGON'")
  } else if (sf::st_geometry_type(frame) != "POLYGON") {
    stop("frame must be an sf object of geometry type 'POLYGON'")
  }
  
  if (!(class(n_points) %in% c("numeric", "integer")) | length(n_points) > 1) {
    stop("n_points must be a single numeric value")
  }
  
  if (!(sample_type %in% c("simple", "balanced", "cluster"))) {
    stop("sample_type must be 'simple', 'balanced', or 'cluster'")
  }
  
  if (!is.null(seed_number)) {
    if (!(class(seed_number) %in% c("numeric", "integer")) | length(seed_number) > 1) {
      stop("seed_number must be a single numeric value")
    }
  } else {
    seed_number <- sample(x = 1:9999999,
                          size = 1)
  }
  
  if (is.null(projection)) {
    projection <- sf::st_crs(polygon)
  } else if (class(projection) == "character") {
    projection <- sp::CRS(projection)
  } else if (!("CRS" %in% class(projection))) {
    stop("projection must either be a valid PROJ4 string or CRS object")
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
                                  type = "random")
                   },
                   "balanced" = {
                     # grts() requires a design object with this format
                     # The "None" is because there are no strata to name
                     design <- list(None = list(panel = c("1" = n_points),
                                                seltype = "Equal",
                                                over = 0))
                     set.seed(seed_number)
                     spsurvey::grts(design = design,
                                    DesignID = "",
                                    type.frame = "area",
                                    src.frame = "sf.object",
                                    sf.object = frame,
                                    shapefile = FALSE)
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
                             1:nrow(points))
  
  points$sample_seed <- seed_number
  
  output <- points[, c("sample_id", "sample_seed")]
  
  return(output)
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
                        verbose = FALSE){
  if (!is.numeric(counts) | length(counts) < 2) {
    stop("counts must be a numeric vector with at least two values")
  }
  
  if (!(chisq %in% c("A", "B", "best"))) {
    stop("The only valid values for chisq are 'A', 'B', and 'best'.")
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
                           "A" = {chisq_quantiles["A"]},
                           "B" = {chisq_quantiles["B"]},
                           "best" = {
                             pick <- which.min(chisq_quantiles)
                             if (verbose){
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