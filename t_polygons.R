# Consider shifting to spatstat::dirichlet() for tessellation?
# https://mgimond.github.io/Spatial/interpolation-in-r.html

t_polygons <- function(coords,
                       extent,
                       projection) {
  if (class(coords) != "data.frame" | !all(c("X", "Y") %in% names(coords))) {
    stop("coords must be a data frame with the variables 'X' and 'Y'")
  }
  
  if (class(projection) == "character") {
    projection <- sp::CRS(projection)
  } else if (class(projection) != "CRS") {
    stop("projection must either be a character string to pass to sp::CRS() or a CRS object")
  }
  
  if (class(extent) == "Extent") {
    extent_vector <- sapply(X = 1:4,
                            extent = extent,
                            FUN = function(X, extent){
                              extent[X]
                            })
    names(extent_vector) <- c("xmin", "xmax", "ymin", "ymax")
  } else if (class(extent) == "numeric") {
    if (!all(c("xmin", "xmax", "ymin", "ymax") %in% names(extent))) {
      stop("extent as a numeric vector must have four values named xmin, xmax, ymin, and ymax")
    } else {
      extent_vector <- extent[c("xmin", "xmax", "ymin", "ymax")]
    }
  } else {
    stop("extent must be either a named numeric vector or an extent class object")
  }
  
  # What if the points don't fall within the extent?
  if (min(coords[["X"]]) < extent_vector["xmin"]) {
    stop ("The minimum X coordinate value is outside the extent")
  }
  if (max(coords[["X"]]) > extent_vector["xmax"]) {
    stop ("The maximum X coordinate value is outside the extent")
  }
  if (min(coords[["Y"]]) < extent_vector["ymin"]) {
    stop ("The minimum Y coordinate value is outside the extent")
  }
  if (max(coords[["Y"]]) > extent_vector["ymax"]) {
    stop ("The maximum Y coordinate value is outside the extent")
  }
  
  # Triangulate
  tessellation <- deldir::deldir(x = coords[["X"]],
                                 y = coords[["Y"]],
                                 rw = extent_vector)
  
  # Convert to sets of coordinates where bisection lines intersect
  tessellation_tiles <- deldir::tile.list(tessellation)
  
  # For each set of vertex coordinates, make the polygon
  polygons_list <- lapply(X = 1:length(tessellation_tiles),
                          tiles = tessellation_tiles,
                          FUN = function(X, tiles){
                            current_tile <- tiles[[X]]
                            # Get the coordinates
                            vertex_coordinates <- cbind(current_tile[["x"]],
                                                        current_tile[["y"]])
                            # Create a closed polygon by adding the first coordinates to the end of the data frame
                            vertex_coordinates <- rbind(vertex_coordinates,
                                                        vertex_coordinates[1, ])
                            # Make polygons
                            sp::Polygons(list(sp::Polygon(vertex_coordinates)),
                                         ID = as.character(X))
                          })
  
  # Make polygons!
  polygons <- sp::SpatialPolygons(polygons_list,
                                  proj4string = projection)
  # Get the IDs so that we can associate a data frame with the polygons
  row_names <- sapply(slot(polygons, 'polygons'),
                      function(x){
                        slot(x, 'ID')
                      })
  # Convert to an SPDF
  thiessen_polygons <- sp::SpatialPolygonsDataFrame(polygons,
                                                    data = data.frame(id = row_names,
                                                                      x = coords[["X"]],
                                                                      y = coords[["Y"]],
                                                                      row.names = row_names
                                                    )
  )
}