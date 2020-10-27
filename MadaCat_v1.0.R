#################################################################################
### MadaCAT: Interactive Preliminary IUCN Red List Assessments for Madagascar ###
#################################################################################

# author: simon.crameri@env.ethz.ch, July 2020

## Install (if needed) and load libraries
load.pkgs <- c("shiny",
               "shinythemes",
               "shinyWidgets",
               "sf",
               "stars",
               "smoothr",
               "rgeos",
               "raster",
               "data.table",
               "parallel",
               "snow",
               "KernSmooth",
               "MASS",
               "ConR",
               "rnaturalearthdata",
               "alphahull",
               "rasterVis",
               "latticeExtra",
               "RColorBrewer",
               "ggplot2",
               "forecast",
               "plotly",
               "leaflet",
               "leaflet.extras",
               "htmlwidgets"
)
new.pkgs <- load.pkgs[!load.pkgs %in% installed.packages()]
if (length(new.pkgs) >= 1) {
  for (j in new.pkgs) {
    cat("installing", j, "\n") ; install.packages(j)
  }
}
for (j in load.pkgs) suppressPackageStartupMessages(library(j, character.only = TRUE))


## Define basic app parameters
# user interface
now <- as.numeric(strsplit(as.character(Sys.time()), split = "-")[[1]][1]) # current year
plot.height = 1230 # 1100
plot.width = 10
sidebar.width = 2
input.height = 65
altrangeslider = FALSE
aooslider = FALSE
periodcheckbox = TRUE

# map appearance
col.historic = "#FF0000"
col.extant = "#0000FF"
contour.margin = 50 # needed to close contour lines of sampling density
contour.sqrt = TRUE # if TRUE, contour lines will be based on sqrt-transformed kernel densities
providers <- c("Esri.WorldImagery","OpenStreetMap.Mapnik","Esri.WorldTopoMap","Esri.WorldGrayCanvas") #≤#≤# leaflet map background tiles
palette.altitude <- c("deepblue", rev(brewer.pal(n = 11, name = "RdYlGn")), "grey37","grey50","grey62","grey75","white")
palette.breaks <- c(seq(-100 ,1000, by = 100), seq(1200, 2000, by = 200), 3000)
scale.position <- "topright" # "topleft" # "bottomleft"

# behaviour
verbose = TRUE
print.date = FALSE
nb.cells = FALSE
l.maxpixels = Inf # for leaflet
p.maxpixels = 5E05 # for lattice
buffer.method = "xy" # "xy" allows separate definition of x and y margins, or "x" (same margin for x and y, quadratic extent)
round.percent <- 2 #≤#≤# display this number of digits for percent habitat change
add.layers = TRUE # if TRUE, will add ecoregions and protected areas onto the leaflet map
add.sampling.effort = TRUE # if TRUE, will add black circle markers for Malagasy collections (227,884) inside the chosen extent e()
get.subpop <- TRUE # if TRUE, will add subpopulation radii
get.aoo.grid <- TRUE # if TRUE, will compute and display AOO grids
get.loc.grid <- TRUE # if TRUE, will compute and display location grids

## Define ConR parameters (IUCN AOO and EOO calcluation)
var.extant <- "occurrenceRemarks" # extant / historic information is looked for in this variable (can be NULL)
str.extant <- c("Extant: 1", "1", "extant", "Extant") # these strings are considered to denote extant collections
var.date <- "eventDate" # column with collection year information (can be NULL)

# EOO
country_map <- "border"

# AOO
nbe.rep.rast.AOO <- 200 # the larger, the easier to find the minimum number of occupied grid cells (AOO), the longer the calculations

# locations
method_locations <- "fixed_grid" # "sliding scale" (grid size = Rel_cell_size * max interoccurrence distance) OR "fixed_grid"

# map
DrawMap = FALSE
add.legend = TRUE
file_name = "IUCN_"
map_pdf = FALSE #TRUE

# general
showWarnings = TRUE
export_shp = TRUE
write_shp = FALSE
write_results = FALSE
write_file_option = "excel"


## Define geographical parameters
raster.folder <- "original_1360943464_utm" #≤#≤# folder with .tif raster files
for.pattern <- "^for[0-9]+.tif" #≤#≤# forest cover (forXXXX) file pattern
dif.pattern <- "^for[0-9]+_vs_for[0-9]+.tif" #≤#≤# forest cover change (forXXXX_vs_forXXXX) file pattern
raster.altitude <- "MDG_alt_1360943464_utm.tif" #≤#≤# path to altitude raster
project.raster <- "+proj=utm +zone=38 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" #≤#≤# projection of forest cover (should be utm)
project.coords <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #≤#≤# projection of input coordinates (should be longlat)

unit <- "kha"  #≤#≤# km^2 or Mha
fac <- 30*30/1E04/1E03  #≤#≤# factor to convert raster pixels to <unit>
projectForLeaflet <- TRUE #≤#≤# if TRUE, will project the raster(s) for leaflet maps, needed if they are not Web Mercator Projection (epsg3857)
jitter.amount <- 0.0003


## Read forest cover, forest cover difference, and altitude rasters
raster.for <- list.files(path = raster.folder, pattern = for.pattern, full.names = TRUE)
if (length(raster.for) == 0) stop("No raster with pattern ", for.pattern, " found in folder ", raster.folder)
raster.dif <- list.files(path = raster.folder, pattern = dif.pattern, full.names = TRUE)
if (length(raster.dif) == 0) stop("No raster with pattern ", dif.pattern, " found in folder ", raster.folder)
forvars <- tools::file_path_sans_ext(basename(raster.for))
difvars <- tools::file_path_sans_ext(basename(raster.dif))
for (i in raster.for) {assign(forvars[which(raster.for == i)], raster(i))}
for (i in raster.dif) {assign(difvars[which(raster.dif == i)], raster(i))}
altitude <- raster(file.path(raster.folder, raster.altitude))
test <- try(compareRaster(get(difvars[1]), altitude, extent = T, rowcol = T, crs = T, res = T, orig = T), silent = TRUE)
if (inherits(test, "try-error")) {
  f <- ceiling(floorsqrt(ncell(get(difvars[1]))) / sqrt(ncell(altitude)))
  if (!all.equal(res(altitude), res(get(difvars[1]))) & f > 1) {
    if (verbose) cat(paste0(Sys.time(), " disaggregating altitude raster...\n"))
    altitude <- disaggregate(altitude, fact = f)
  }
  if (! all.equal(attributes(altitude)$crs, attributes(get(difvars[1]))$crs)) {
    if (verbose) cat(paste0(Sys.time(), " projecting altitude raster...\n"))
    raster.altitude2 <- paste0(file.path(raster.folder, tools::file_path_sans_ext(basename(raster.altitude))), "_proj.tif")
    projectRaster(from = altitude, to = get(difvars[1]), 
                  filename = raster.altitude2)
    altitude <- raster(raster.altitude2)
  }
  test2 <- try(compareRaster(get(difvars[1]), altitude, extent = T, rowcol = T, crs = T, res = T, orig = T), silent = TRUE)
  if (inherits(test2, "try-error")) stop("check conformity of forest and altitude rasters")
}
dvars <- paste0("d.", forvars)
choices <- gsub("^for", "", forvars)

## Read Vieilledent 2018 data (tables with forest cover loss and deforestation rates)
load("d.vieilledent.rda")

## Read protected areas and ecoregions
load("MDG_eco_PA.rda")

## Load sampling effort data
if (add.sampling.effort) load("20200806_MadagascarAngios.rda")

## Define helperfunctions
habitat.stats <- function(x, S, alt, altmin, altmax, e, o = -1, fcc = TRUE, large = 2E+09,
                          n = 0, l = -2, c = 1, g = 2) {
  
  ## Parameters
  # x      rasterLayer of forest cover difference (fcc = TRUE) or forest cover (fcc = FALSE)
  # S      SpatialPolygons defining the extent of interest
  # alt    rasterLayer of altitude (needs to have the same dim/res/ext/crs as x)
  # altmin minimum altitude
  # altmax maximum altitude
  # e      Extent object defining the square to crop (S must lie within e, e must be inside x)
  # o      value assigned to region outside extent of interest (must be a unique value)
  # fcc    if TRUE, will also count number of cells that are l or g
  # n      value for non-habitat
  # l      value for habitat loss
  # c      value for constant habitat
  # g      value for habitat gain
  
  ## crop to extent of interest
  r <- crop(x, e)
  
  
  ## mask forest using custom function (and layers)
  if (ncell(r) >= large) mode <- "clusterR" else mode <- "simple"
  switch(mode,
         simple = {
           ## arithmetic solution (for small rasters)  
           r[r != n & (alt > altmax | alt < altmin)] <- o
         },
         clusterR = {
           # clusterR solution (seems to work)
           mask.habitat <- function(x) {
             if (x[1] != n & !is.na(x[2]) & (x[2] > altmax | x[2] < altmin)) {
               o
             } else {
               x[1]
             }
           }
           n <<- n
           o <<- o
           altmin <<- altmin
           altmax <<- altmax
           beginCluster()
           r <- clusterR(stack(r, alt), calc, args = list(fun = mask.habitat), export = c("n","o","altmin","altmax"))
           endCluster()
         },
         overlay = {
           ## overlay solution (does not work)
           mask.habitat.overlay <- function(x, y, n, o, altmin, altmax, ...) {
             if (x != n & !is.na(y) & (y > altmax | y < altmin)) {
               o
             } else {
               x
             }
           }
           r <- overlay(x = r, y = alt, fun = function(x, y, ...) {mask.habitat.overlay(x, y, n = n, o = o, altmin = altmin, altmax = altmax)})
         },
         calc = {
           ## calc solution (not stable)
           mask.habitat.calc <- function(x, n, o, altmin, altmax, ...) {
             if (x[1] != n & !is.na(x[2]) & (x[2] > altmax | x[2] < altmin)) {
               o
             } else {
               x[1]
             }
           }
           r <- calc(stack(r, alt), fun = function(x, ...) {mask.habitat.calc(x, n = n, o = o, altmin = altmin, altmax = altmax)})
         }
  )
  
  
  ## mask forest outside extent
  m <- mask(r, S, updatevalue = NA)
  
  ## calculate habitat change
  if (fcc) {
    
    loss  <- extract(m, extent(m), fun = function(x,...) sum(x == l, na.rm = TRUE))
    non1  <- extract(m, extent(m), fun = function(x,...) sum(x == n, na.rm = TRUE))
    non2  <- extract(m, extent(m), fun = function(x,...) sum(x == o, na.rm = TRUE))
    const <- extract(m, extent(m), fun = function(x,...) sum(x == c, na.rm = TRUE))
    gain  <- extract(m, extent(m), fun = function(x,...) sum(x == g, na.rm = TRUE))
    
    res <- list(raster = r, HabitatLoss = loss, HabitatConstant = const, HabitatGain = gain, outForest = non1, outAltitude = non2)
  } else {
    
    const <- extract(m, extent(m), fun = function(x,...) sum(x == c, na.rm = TRUE))
    non1  <- extract(m, extent(m), fun = function(x,...) sum(x == n, na.rm = TRUE))
    non2  <- extract(m, extent(m), fun = function(x,...) sum(x == o, na.rm = TRUE))
    
    res <- list(HabitatLoss = NA, HabitatConstant = const, HabitatGain = NA, outForest = non1, outAltitude = non2)
  }
  
  return(res)
}
habitat.wrapper <- function(year, r, S, alt, altmin, altmax, e, period1, period2, verbose = TRUE) {
  if (year %in% c(period1, period2)) {
    if (verbose) cat(paste0(Sys.time(), " using habitat change instead of ", year, " raster...\n"))
    return(r[-1])
  } else {
    if (verbose) cat(paste0(Sys.time(), " processing ", year, " raster...\n"))
    habitat.stats(x = get(paste0("for", year)), S = S, alt = alt, 
                  altmin = altmin, altmax = altmax, e = e, fcc = FALSE)
  }
}
get.stats <- function(d.res, factor = 1, unit = "Pixels", period1, period2, refyear = 1, extent = "extent", use.pairwise = TRUE) {
  d.res <- data.frame(t(matrix(do.call(cbind, d.res), nrow = 5, dimnames = attr(d.res, "dimnames"))))
  d.res$Period <- gsub("^d.", "", rownames(d.res))
  d.res$Year <- as.numeric(gsub("^for", "", d.res$Period))
  habitat <- nonhabitat <- numeric(length = nrow(d.res))
  names(habitat) <- names(nonhabitat) <- rownames(d.res)
  for (i in seq(nrow(d.res))) {
    if (use.pairwise & d.res[i,"Year"] == period1) {
      ls1 <- d.res[i,"HabitatConstant"] + d.res[i,"HabitatLoss"]
      ls2 <- d.res[i,"outForest"] + d.res[i,"outAltitude"] + d.res[i,"HabitatGain"]
      d.res[i,"HabitatGain"] <- NA
    } else if (use.pairwise & d.res[i,"Year"] == period2) {
      ls1 <- d.res[i,"HabitatConstant"] + d.res[i,"HabitatGain"]
      ls2 <- d.res[i,"outForest"] + d.res[i,"outAltitude"] + d.res[i,"HabitatLoss"]
      d.res[i,"HabitatLoss"] <- NA
    } else {
      ls1 <- d.res[i,"HabitatConstant"]
      ls2 <-d.res[i,"outForest"] + d.res[i,"outAltitude"]
    }
    habitat[rownames(d.res)[i]] <- ls1
    nonhabitat[rownames(d.res)[i]] <- ls2
  }
  d.hab <- data.frame(Year = d.res$Year, Extent = extent,
                      NonHabitat_pixels = nonhabitat, value1 = nonhabitat*factor,
                      Habitat_pixels = habitat, value2 = habitat*factor,
                      HabitatLoss_pixels = d.res$HabitatLoss, value3 = d.res$HabitatLoss*factor,
                      HabitatGain_pixels = d.res$HabitatGain, value4 = d.res$HabitatGain*factor,
                      HabitatFraction_perc = 100*habitat/(habitat+nonhabitat),
                      Remaining_perc = 100*habitat/habitat[refyear],
                      Reduction_perc = 100 - (100*habitat/habitat[refyear]))
  names(d.hab)[4] <- paste0("NonHabitat_", unit)
  names(d.hab)[6] <- paste0("Habitat_", unit)
  names(d.hab)[8] <- paste0("HabitatLoss_", unit)
  names(d.hab)[10] <- paste0("HabitatGain_", unit)
  return(d.hab)
}
get.iucn <- function(DATA, country_map, protec.areas, ID_shape_PA,
                     method.range, exclude.area, alpha, buff.alpha, title = "",
                     Cell_size_AOO, nbe.rep.rast.AOO, get.aoo.grid = FALSE,
                     Resol_sub_pop, get.loc.grid = FALSE,
                     method_protected_area, method_locations, Cell_size_locations, Rel_cell_size,
                     DrawMap, add.legend, file_name, map_pdf, showWarnings, export_shp, write_shp, 
                     write_results, write_file_option) {
  
  # check
  if (exclude.area & !is.null(country_map)) {
    check <- rowSums(sapply(1:length(country_map@polygons[[1]]@Polygons), FUN = function(x) {
      point.in.polygon(point.x = DATA[,2], point.y = DATA[,1], 
                       pol.x = country_map@polygons[[1]]@Polygons[[x]]@coords[,1], 
                       pol.y = country_map@polygons[[1]]@Polygons[[x]]@coords[,2])
    }))
    
    dout <- sum(check == 0)
    din <- sum(check == 1)
    dna <- sum(is.na(check))
    msg <- paste0(dout, " / ", length(check), " (", round(100*dout/length(check),2), "%) ",
                  "of occurrences lie outside of the <country_map> polygon")
    if (any(check == 0)) {
      warning(msg)
    } 
    if (all(check == 0)) {
      warning("reverting to <country_map> = land")
      country_map <- land
    }
  }
  
  # get IUCN.eval
  sink("/dev/null")
  u <- try(IUCN.eval(DATA,
                 # input shape files
                 country_map = country_map, protec.areas = protec.areas, ID_shape_PA = ID_shape_PA,
                 # EOO
                 method.range = method.range, exclude.area = exclude.area,
                 alpha = alpha, buff.alpha = buff.alpha,
                 # AOO
                 Cell_size_AOO = Cell_size_AOO, nbe.rep.rast.AOO = nbe.rep.rast.AOO,
                 # subpopulations
                 Resol_sub_pop = Resol_sub_pop,
                 # locations
                 method_locations = method_locations,
                 method_protected_area = method_protected_area,
                 Cell_size_locations = Cell_size_locations,
                 Rel_cell_size = Rel_cell_size,
                 # map
                 DrawMap = DrawMap, add.legend = add.legend, file_name = file_name, map_pdf = map_pdf,
                 # general
                 showWarnings = showWarnings,
                 export_shp = export_shp, write_shp = write_shp,
                 write_results = write_results, write_file_option = write_file_option)[[1]], silent = TRUE)
  if (inherits(u, "try-error") & !is.null(protec.areas) & buff.alpha >= 0.01) {
    warning(gsub("  ", " ", gsub("\n", " ", u[1])), "\nreverting to <protec.areas = NULL>")
    protec.areas = NULL

    u <- try(IUCN.eval(DATA,
                       # input shape files
                       country_map = country_map, protec.areas = protec.areas, ID_shape_PA = ID_shape_PA,
                       # EOO
                       method.range = method.range, exclude.area = exclude.area,
                       alpha = alpha, buff.alpha = buff.alpha,
                       # AOO
                       Cell_size_AOO = Cell_size_AOO, nbe.rep.rast.AOO = nbe.rep.rast.AOO,
                       # subpopulations
                       Resol_sub_pop = Resol_sub_pop,
                       # locations
                       method_locations = method_locations,
                       method_protected_area = method_protected_area,
                       Cell_size_locations = Cell_size_locations,
                       Rel_cell_size = Rel_cell_size,
                       # map
                       DrawMap = DrawMap, add.legend = add.legend, file_name = file_name, map_pdf = map_pdf,
                       # general
                       showWarnings = showWarnings,
                       export_shp = export_shp, write_shp = write_shp,
                       write_results = write_results, write_file_option = write_file_option)[[1]], silent = TRUE)
  }
  if (inherits(u, "try-error") & buff.alpha < 0.01) {
    warning(gsub("  ", " ", gsub("\n", " ", u[1])), "\nreverting to <buff.alpha = 0.01>")
    buff.alpha = 0.01
    
    u <- try(IUCN.eval(DATA,
                       # input shape files
                       country_map = country_map, protec.areas = protec.areas, ID_shape_PA = ID_shape_PA,
                       # EOO
                       method.range = method.range, exclude.area = exclude.area,
                       alpha = alpha, buff.alpha = buff.alpha,
                       # AOO
                       Cell_size_AOO = Cell_size_AOO, nbe.rep.rast.AOO = nbe.rep.rast.AOO,
                       # subpopulations
                       Resol_sub_pop = Resol_sub_pop,
                       # locations
                       method_locations = method_locations,
                       method_protected_area = method_protected_area,
                       Cell_size_locations = Cell_size_locations,
                       Rel_cell_size = Rel_cell_size,
                       # map
                       DrawMap = DrawMap, add.legend = add.legend, file_name = file_name, map_pdf = map_pdf,
                       # general
                       showWarnings = showWarnings,
                       export_shp = export_shp, write_shp = write_shp,
                       write_results = write_results, write_file_option = write_file_option)[[1]], silent = TRUE)

  }
  if (inherits(u, "try-error")) {
    warning(gsub("  ", " ", gsub("\n", " ", u[1])), "\nreverting to <protec.areas> = NULL & <buff.alpha = 0.01>")
    protec.areas = NULL
    buff.alpha = 0.01
    
    u <- try(IUCN.eval(DATA,
                       # input shape files
                       country_map = country_map, protec.areas = protec.areas, ID_shape_PA = ID_shape_PA,
                       # EOO
                       method.range = method.range, exclude.area = exclude.area,
                       alpha = alpha, buff.alpha = buff.alpha,
                       # AOO
                       Cell_size_AOO = Cell_size_AOO, nbe.rep.rast.AOO = nbe.rep.rast.AOO,
                       # subpopulations
                       Resol_sub_pop = Resol_sub_pop,
                       # locations
                       method_locations = method_locations,
                       method_protected_area = method_protected_area,
                       Cell_size_locations = Cell_size_locations,
                       Rel_cell_size = Rel_cell_size,
                       # map
                       DrawMap = DrawMap, add.legend = add.legend, file_name = file_name, map_pdf = map_pdf,
                       # general
                       showWarnings = showWarnings,
                       export_shp = export_shp, write_shp = write_shp,
                       write_results = write_results, write_file_option = write_file_option)[[1]], silent = TRUE)
    
  }
  if (inherits(u, "try-error")) {
    sink()
    stop(gsub("  ", " ", gsub("\n", " ", u[1])), "\nplease use a different parameter combination")
  }
  sink()
  
  # get AOO polygons
  if (get.aoo.grid) {
    u$spatialPoly_AOO <- AOO.computing(XY = DATA, Cell_size_AOO = Cell_size_AOO,
                                       nbe.rep.rast.AOO = 0,
                                       export_shp = TRUE, show_progress = FALSE)[[2]][[1]]
  }
  if (get.loc.grid) {
    u$spatialPoly_Loc <- locations.comp(XY = DATA, method = method_locations, nbe_rep = 0, 
                                        method_protected_area = method_protected_area, 
                                        protec.areas = protec.areas, ID_shape_PA = ID_shape_PA,
                                        Cell_size_locations = Cell_size_locations, Rel_cell_size = Rel_cell_size, 
                                        show_progress = FALSE)[[1]][[1]]
  }
  
  # add number of subpopulations
  u$spatialPoly_subpop@data$Subpopulation <- 1:u$Results[[1]][[4]]
  
  # compile results
  iucn.res <- data.frame(t(u$Results), stringsAsFactors = FALSE)
  if (is.null(protec.areas)) {
    iucn.res <- data.frame(iucn.res[,c(1:5)], lPA = NA, iucn.res[,c(6:7)], pPA = NA, iucn.res[,c(8:ncol(iucn.res))])
  }
  names(iucn.res) <- c("EOO","AOO","# unique occurrences", "# subpopulations", "# locations", "# locations in PA", 
                       "IUCN Cat B", "IUCN Code", "Occurrences in PA", "IUCN Cat AOO", "IUCN Cat EOO")
  iucn.res[,c(1,2)] <- paste0(iucn.res[,c(1,2)], " km<sup>2</sup>")
  iucn.res[,9] <- paste0(iucn.res[,9], "%")
  
  # create popups
  popup.eoo <- paste0("<b>", title, "</b>",
                      "<br/>", get.popup(iucn.res, names(iucn.res)[c(1:2)]),
                      "Min. # occupied AOO cells: ", as.numeric(u$Results[[1]][2])/Cell_size_AOO^2, "<br/>",
                      "<br/>", get.popup(iucn.res, names(iucn.res)[c(3:6,9)]),
                      "<br/>", get.popup(iucn.res, names(iucn.res)[c(7:8,10:11)]),
                      "<br/>Parameters:",
                      "<br/>Country map: ", ifelse(is.null(country_map), NA, paste(round(as.vector(country_map@bbox), 2), collapse = ",")),
                      "<br/>PA map: ", ifelse(is.null(protec.areas), NA, length(protec.areas)), " ", ID_shape_PA,
                      "<br/>EOO exclude sea: ", exclude.area,
                      "<br/>EOO method: ", method.range,
                      "<br/>EOO alpha (radius): ", alpha,
                      "<br/>EOO buffer: ", buff.alpha, " degrees",
                      "<br/>AOO cell width: ", Cell_size_AOO, " km",
                      "<br/>AOO iterations: ", nbe.rep.rast.AOO,
                      "<br/>Subpopulation radius: ", Resol_sub_pop, " km",
                      "<br/>PA method: ", method_protected_area,
                      "<br/>Location method: ", method_locations,
                      "<br/>Location cell size: ", Cell_size_locations, " km",
                      "<br/>Location multiplyer: ", Rel_cell_size, " * max. distance km<br/>")
  popup.sp <- paste0("<b>", title, "</b>",
                     "<br/>", gsub("([A-Za-z]+): ([0-9]+)", paste("\\1: \\2 /", length(u$spatialPoly_subpop)), 
                                   get.popup(u$spatialPoly_subpop@data, names(u$spatialPoly_subpop))))
  popup.aoo <- paste0("<b>", title, "</b>",
                      "<br/>Min. AOO: ", u$Results[[1]][2], " km<sup>2</sup>", 
                      "<br/>Min. # occupied AOO cells: ",  as.numeric(u$Results[[1]][2])/Cell_size_AOO^2,
                      "<br/><br/>Max. AOO: ", length(u$spatialPoly_AOO)*Cell_size_AOO^2, " km<sup>2</sup>",
                      "<br/>Max. # occupied AOO cells: ", length(u$spatialPoly_AOO),
                      "<br/><br/>Cell width: ", Cell_size_AOO, " km<br/># iterations: ", nbe.rep.rast.AOO)
  
  # return results
  res <- list(DATA = DATA,
              Results = data.frame(t(u$Results), stringsAsFactors = FALSE),
              iucn.res = iucn.res,
              spatialPoly_EOO = u$spatialPoly_EOO,
              spatialPoly_AOO = u$spatialPoly_AOO,
              spatialPoly_Loc = u$spatialPoly_Loc,
              spatialPoly_subpop = u$spatialPoly_subpop,
              popup.eoo = popup.eoo,
              popup.aoo = popup.aoo,
              popup.sp = popup.sp)
  return(res)
}
get.DATA <- function(x, var.date = "eventDate", var.extant = "occurrenceRemarks", str.extant = c("Extant: 1", "1", "extant", "Extant")) {
  if (is.null(var.date)) {
    x[,"eventDate"] <- NA
  } else {
    if (var.date %in% names(x)) {
      x[,"eventDate"] <- x[, var.date]
    } else {
      x[,"eventDate"] <- NA
    }
  }
  if (is.null(var.extant)) {
    x[,"occurrenceRemarks"] <- "1"
  } else {
    if (var.extant %in% names(x)) {
      x[,"occurrenceRemarks"] <- x[, var.extant]
    } else {
      x[,"occurrenceRemarks"] <- "1"
    }
  }
  x[,"extant"] <- factor(ifelse(x[,"occurrenceRemarks"] %in% str.extant, "1", "0"), levels = c("0","1"))
  XY <- suppressWarnings(
    data.frame(ddlat = as.numeric(as.character(x[,1])), ddlon = as.numeric(as.character(x[,2])),
               tax = "Taxon", family = "Family",
               coly = as.numeric(sapply(strsplit(as.character(x[,"eventDate"]), split = "/"), "[", 3)),
               occurrenceRemarks = x[,"occurrenceRemarks"],
               extant = x[,"extant"],
               x[,-c(1:2)][,!names(x[-c(1:2)]) %in% c("occurrenceRemarks","extant")])
  )
  
  # remove observations without coordinates
  nocoords <- which(is.na(XY[,1]) | is.na(XY[,2]))
  if (length(nocoords) > 0) {
    msg <- paste0("removed ", length(nocoords), " / ", nrow(XY), " (", round(100*length(nocoords)/nrow(XY), 2), "%) ",
                  "observation(s) with missing coordinates")
    warning(msg)
    XY <- XY[-nocoords,]
  }
  
  # warn if the coordinates are not between -180 <=  lon <= 180 and -90 <= lat <= 90
  wronglat <- which(abs(XY[,2]) > 90)
  wronglon <- which(abs(XY[,1]) > 180)
  wrong <- unique(c(wronglon, wronglat))
  if (length(wronglat) > 0) {
    msg <- paste0("removed ", length(wronglat), " / ", nrow(XY), " (", round(100*length(wronglat)/nrow(XY), 2), "%) ",
                  " observation(s) with abs(latitude)  > 90")
    warning(msg)
  }
  if (length(wronglon) > 0) {
    msg <- paste0("removed ", length(wronglon), " / ", nrow(XY), " (", round(100*length(wronglon)/nrow(XY), 2), "%) ",
                  " observation(s) with abs(longitude)  > 180")
    warning(msg)
  }
  if (length(wrong) > 0) XY <- XY[-wrong,]
  if (nrow(XY) == 0) {
    stop("Zero observations remain")
  }
  return(XY)
}
get.popup <- function(x, names) {
  popups <- character()
  for (i in seq(nrow(x))) {
    popup <- character()
    for (name in names) {
      popup <- paste0(popup, name, ": ", x[i,name], sep = "<br/>")
    }
    popups[i] <- popup
  }
  popups
}
get.popups <- function(x, method = "darwin", verbose = FALSE) {
  popups <- character()
  link <- unique(gsub("[0-9]+$", "", names(unlist(sapply(names(x), function(y) grep("http://", x[,y]))))))[1]
  switch(method,
         darwin = {
           for (i in seq(nrow(x))) {
             if (verbose) cat(i, "/", nrow(x), "\n")
             popups[i] <- paste(
               # collector
               paste0("Collector: ", x[i, "collector"]),
               paste0("Code: ", x[i, "institutionCode"]),
               "",
               
               # coordinates
               paste0("Latitude: ", x[i,1]),
               paste0("Longitude: ", x[i,2]),
               paste0("Coordinate Method: ", x[i,"coordinateuncertaintyinmeters"]),
               "",
               paste0("Elevation: ", x[i,"verbatimElevation"]),
               paste0("Date: ", x[i,"eventDate"]),
               "",
               
               # locality
               paste0("Locality: ", x[i,"locality"]),
               "",
               
               # url
               paste0("Remarks: ", x[i,"occurrenceRemarks"]),
               ifelse(!is.null(x[i,link]), paste0("<a href='", paste0(x[i,link]), "', target=\"_blank\">", "visit URL", "</a>"), "no ULR"),
               sep = "<br/>")
           }
         },
         unknown = {
           for (i in seq(nrow(x))) {
             if (verbose) cat(i, "/", nrow(x), "\n")
             b <- character()
             for (j in colnames(x)[!colnames(x) == link]) {
               b <- paste(b, paste0(j, ": ", x[i,j]), sep = "<br/>")
             }
             popups[i] <- paste(b,"",
                                ifelse(!is.null(x[i,link]), paste0("<a href='", paste0(x[i,link]), "', target=\"_blank\">", "visit URL", "</a>"), "no ULR"),
                                sep = "<br/>")
           }
         })
  
  return(popups)
}
contour.popups <- function(x, mat, b, thr, n) {
  get.freqs <- function(y, x, mat, b, l) {
    dy <- as.data.frame(table(cut(x[mat[,y]], breaks = b, labels = l)))
    return(as.character(apply(dy[dy$Freq != 0,], 1, paste, collapse = ": ")))
  }
  l <- character()
  for (i in seq(length(b))[-length(b)]) {
    if (b[i] < thr) {
      l <- c(l, paste(b[i]+1, b[i+1], sep = "-"))
    } else {
      l <- c(l, b[i]+1)
    }
  }
  ages <- sapply(seq(length(n)), get.freqs, x = x, mat = mat, b = b, l = l)
  popups <- character()
  for (i in seq(length(n))) {
    popups[i] <- paste(paste0("N = ", n[i]), "",
                       paste(ages[[i]], collapse = "<br/>"),
                       sep = "<br/>")
  }
  return(popups)
}
get.contour <- function(sp, margin, nlevels = 10, sqrt = FALSE, nb = 5000, resolution = 1000) {
  
  # define bandwidth for kernel density estimate (MASS package: Sheather & Jones method)
  h <- c(width.SJ(sp@coords[,1], method = "dpi", nb = nb),
         width.SJ(sp@coords[,2], method = "dpi", nb = nb))
  
  # estimate 2D kernel density (KernSmooth package)
  kde <- bkde2D(sp@coords[,1:2], bandwidth = h,
                gridsize = c(round(resolution*diff(range(sp@coords[,1]))),
                             round(resolution*diff(range(sp@coords[,2])))))
  
  # define a margin around the kernel density extent (prevents edge effects)
  dx1 <- diff(head(kde$x1,2))
  dx2 <- diff(head(kde$x2,2))
  l1 <- seq((kde$x1[1]-margin*dx1), (kde$x1[1]-dx1), by = dx1)
  t1 <- seq((kde$x1[length(kde$x1)]+dx1), (kde$x1[length(kde$x1)]+margin*dx1), by = dx1)
  l2 <- seq((kde$x2[1]-margin*dx2), (kde$x2[1]-dx2), by = dx2)
  t2 <- seq((kde$x2[length(kde$x2)]+dx2), (kde$x2[length(kde$x2)]+margin*dx2), by = dx2)
  hmat <- matrix(data = 0, nrow = margin, ncol = ncol(kde$fhat))
  vmat <- matrix(data = 0, nrow = nrow(kde$fhat)+2*margin, ncol = margin)
  if (sqrt) {
    nmat <- cbind(vmat, rbind(hmat, sqrt(kde$fhat), hmat), vmat)
  } else {
    nmat <- cbind(vmat, rbind(hmat, kde$fhat, hmat), vmat)
  }
  
  # get contour lines
  CL <- contourLines(c(l1,kde$x1,t1), c(l2,kde$x2,t2) , nmat, nlevels = nlevels)
  return(list(h = h, kde = kde, CL = CL, margin = margin, nb = nb, resolution = resolution))
}
contour.polygons <- function(CL, sp, crs) {
  
  # extract contour line levels
  # https://gis.stackexchange.com/questions/168886/r-how-to-build-heatmap-with-the-leaflet-package
  LEVS <- as.factor(sapply(CL, "[[", "level")) # density level for each contour polygon
  NLEV <- length(levels(LEVS)) # number of contour lines
  
  # convert contour lines to polygons
  pgons <- lapply(1:length(CL), function(i) {Polygons(list(Polygon(cbind(CL[[i]]$x, CL[[i]]$y))), ID=i)})
  spgons <- SpatialPolygons(pgons)
  attributes(spgons)$proj4string <- CRS(crs)
  
  # count samples per contour polygon
  pinpgon <- gContains(spgons, sp, byid = TRUE)
  nperpgon <- colSums(pinpgon)
  return(list(CL = CL, LEVS = LEVS, NLEV = NLEV, spgons = spgons, pinpgon = pinpgon, nperpgon = nperpgon, crs = crs))
}
get.raster.contour <- function(r, level, method = "continuous", preprocessed = FALSE, large = FALSE, nchunks = 6,
                               drop.crumbs = TRUE, thr = 100, project = FALSE, proj4string = r@crs@projargs,
                               verbose = FALSE) {
  
  library(stars)
  stopifnot(method %in% c("continuous","discrete"))
  if (verbose) cat(paste0("Polygonizing level ", ifelse(method == "continuous", ">= ", ""), level, "...\n"))
  
  # create raster of level
  if (preprocessed) {
    l <- r
  } else {
    if (large) {
      replace.val <- function(x, method = "discrete", value) { 
        switch(method,
               continuous = {
                 if(x >= value) {
                   return(1)
                 } else {
                   return(NA)
                 }
               },
               discrete = {
                 if(x == value) {
                   return(1)
                 } else {
                   return(NA)
                 }
               })
      }
      
      l <- calc(r, fun = function(x){replace.val(x, method = method, value = level)}, 
                filename = paste0("tmp_level_", level))
    } else {
      l <- raster(extent(r), nrows = nrow(r), ncols = ncol(r), crs = crs(r))
      values(l) <- NA
      
      switch(method,
             continuous = {
               l[r[] >= level] <- 1
             },
             discrete = {
               l[r[] == level] <- 1
             })
    }
  }
   
  
  # create collection of POLYGONs (st_as_sf on stars object is much faster than rasterToPolygons)
  if (large) {
    
    # create chunks
    xmin <- l@extent@xmin
    xmax <- l@extent@xmax
    ymin <- l@extent@ymin
    ymax <- l@extent@ymax
    delta <- floor((ymax-ymin)/(sqrt(nchunks)))
    
    es <- list()
    for (n in 1:nchunks) {
      for (m in 1:nchunks) {
        e <- c(min(xmin+delta*(n-1),xmax), min(xmin+delta*n,xmax), min(ymin+delta*(m-1),ymax), min(ymin+delta*m, ymax))
        if (!(e[1] == e[2] | e[3] == e[4])) {
          es <- c(es, list(e))
        }
      }
    }
    
    # polygonize using st_as_sf(merge = TRUE)
    if (exists("p")) rm(p)
    for (i in seq(length(es))) {
      if (verbose) cat("cropping chunk", i, "/", length(es), "...")
      c <- crop(r, es[[i]])
      if (any(!is.na(values(c)))) {
        if (verbose) cat("polygonizing chunk", i, "\n")
        x <- st_as_stars(c) %>%
          st_as_sf(merge = TRUE)
        if (!exists("p")) {
          p <- x
        } else {
          p <- rbind(p, x)
        }
      } else {
        if (verbose) cat("only NA in chunk", i, "\n")
      }
      invisible(file.remove(c@file@name)) ; rm(c) # prevent accumulation of heavy tmp files
    }
    x <- p ; rm(p)
  } else {
    x <- st_as_stars(l) %>%
      st_as_sf(merge = TRUE)
  }
  
  # drop crumbs
  if (nrow(x) > 0) {
    if (drop.crumbs) {
      library(smoothr)
      if (verbose) cat("dropping crumbs in", nrow(x), "polygons...\n")
      x <- drop_crumbs(x, thr) # in coordinate units
    }
    
    # smooth lines
    if (verbose) cat("smoothing", nrow(x), "polygons...\n")
    s <- smoothr::smooth(x, method = "ksmooth") # ksmooth
    
    # drop NA values
    if (verbose) cat("removing NA coordinates...\n")
    for (i in seq(length(s$geometry))) {
      for (j in seq(length(s$geometry[[i]]))) {
        s$geometry[[i]][[j]] <- na.omit(s$geometry[[i]][[j]])
      }
    }
    
    sp <- try(as_Spatial(s), silent = TRUE)
    if (inherits(sp, "try-error")) {
      if (as.character(sp) == "Error in s$geometry[[i]] : subscript out of bounds\n") {
        warning("level ", level, " not found")
      } else {
        warning(as.character(sp))
      }
      sp <- NA
    } else {
      if (project) {
        sp <- suppressWarnings(spTransform(sp, crs(proj4string)))
      }
      sp@data$level <- level
    }
  } else {
    sp <- NA
  }
  return(sp)
}
raster2SPDF <- function(r, breaks, method = "continuous", preprocessed = FALSE, large = FALSE, nchunks = 6,
                        drop.crumbs = TRUE, thr = 100, 
                        project = FALSE, proj4string = r@crs@projargs,
                        verbose = FALSE) {
  
  # apply get.raster.contour to all levels (breaks[from:to])
  l <- sapply(breaks, FUN = get.raster.contour, r = r, method = method, 
              preprocessed = preprocessed, large = large, nchunks = nchunks,
              drop.crumbs = drop.crumbs, thr = thr, 
              project = project, proj4string = proj4string, verbose = verbose)
  
  # create list of objects "pg" of class "Polygons"
  pgons <- list() 
  for (i in seq(length(l))) {
    if (class(l[[i]]) != "logical") {
      pg <- list()
      for (j in seq(length(l[[i]]@polygons))) {
        pg <- c(pg, Polygon(l[[i]]@polygons[[j]]@Polygons[[1]]@coords))
      }
      pgons <- c(pgons, Polygons(pg, ID = as.character(unique(l[[i]]$level))))
    }
  }
  
  # create  SpatialPolygonsDataFrame from list of Polygons
  cl <- seq(length(l))[unlist(lapply(l, function(x) {ifelse(class(x) != "logical", TRUE, FALSE)}))]
  levs <- sapply(cl, function(i) as.character(unique(l[[i]]$level)))
  spgons <- SpatialPolygonsDataFrame(SpatialPolygons(pgons),
                                     data = data.frame(alt = as.numeric(levs), row.names = levs))
  attributes(spgons)$proj4string <- l[[1]]@proj4string
  
  # return result
  return(spgons)
}
write.raster.level <- function(r, method = "discrete", level, filename = "stdout", format = "GTiff",
                               maxmemory = 25e+09, chunksize = 1e+08, memfrac = 0.9) {
  library(raster)
  library(rgdal)
  rasterOptions(maxmemory = maxmemory, chunksize = chunksize+08, memfrac = memfrac)
  replace.val <- function(x, method = "discrete", level) { 
    switch(method,
           continuous = {
             if(x >= level) {
               return(1)
             } else {
               return(NA)
             }
           },
           discrete = {
             if(x == level) {
               return(1)
             } else {
               return(NA)
             }
           })
  }
  calc(r, fun = function(x){replace.val(x, method = method, level = level)},
       filename = filename, format = format)
}

####################################################################################################

######################
### User Interface ###
######################
ui <- fluidPage(theme = shinythemes::shinytheme("flatly"),
  
  # App title
  titlePanel("MadaCAT: Interactive Preliminary IUCN Red List Assessments for Madagascar"),
  
  fluidPage(
  
    # Sidebar layout with input definitions
    sidebarLayout(
      
      # Sidebar panel for inputs
      sidebarPanel(width = sidebar.width,
        
        ## ICUN Cat A+B options (habitat change)
        h4("Upload observations:", noWS = "before"),
        
        # Input: read .csv
        div(style = paste0("height: ", input.height*0.85, "px;"),
        fileInput(inputId = "filedata", 
                  label = NULL,
                  placeholder = "Select .csv", width = "255px", 
                  accept = c(".csv"))),
        
        # Action Button: calculations are initiated only upon clicking
        div(style = paste0("height: ", input.height*0.75, "px;"),
        actionButton(inputId = "submit",
                     label = "Submit",
                     width = "100%")),
        
        h4("\nIUCN Category A:"),
        
        ## ICUN Cat A Options
        # Input: Period 1 to show
        div(style = paste0("height: ", input.height, "px;"),
        selectInput(inputId = "period1",
                    label = "Define year 1:",
                    choices = gsub("^for", "", forvars),
                    selected = gsub("^for", "", forvars)[1])),
        
        # Input: Period 2 to show
        div(style = paste0("height: ", input.height, "px;"),
        selectInput(inputId = "period2",
                    label = "Define year 2:",
                    choices = gsub("^for", "", forvars),
                    selected = rev(gsub("^for", "", forvars))[1])),
        
        if (periodcheckbox) {
          # Input: Checkbox for selected periods
          checkboxGroupInput(inputId = "periods", 
                             label = "Years for time series:", 
                             choices = choices, inline = TRUE,
                             selected = choices[c(1,3,5,9)])
        },
        
        # Input: Slider for number of forecasting years
        div(style = paste0("height: ", input.height, "px;"),
            sliderInput(inputId = "cast",
                        label = "Forecast (years):",
                        min = 0,
                        max = 100,
                        value = 36,
                        step = 1,
                        ticks = FALSE)),
        # Input: Slider for the altitude range
        if (altrangeslider) {
          div(style = paste0("height: ", input.height*1.2, "px;"),
          sliderInput(inputId = "altrange",
                      label = "Altitude range (m):",
                      min = 0,
                      max = 2800,
                      value = c(0, 2800),
                      step = 10))
        },
        
        if (!altrangeslider) {
          
          # Input: minimum altitude
          div(style = paste0("height: ", input.height*1.2, "px;"),
          textInput(inputId = "altmin",
                    label = "Min. altitude (m):",
                    value = 0,
                    placeholder = "any zero-positive number"))
        },
        
        if (!altrangeslider) {
          
          # Input: maximum altitude
          div(style = paste0("height: ", input.height*1.2, "px;"),
          textInput(inputId = "altmax",
                    label = "Max. altitude (m):",
                    value = 2800,
                    placeholder = "any zero-positive number"))
        },
        
        ## IUCN Cat A + B options (EOO)
        br(),
        h4("IUCN Category A + B:"),
        
        ## EOO
        div(style = paste0("height: ", input.height, "px;"),
            selectInput(inputId = "method.range",
                        label = "EOO method:",
                        choices = c("alpha.hull", "convex.hull"),
                        selected = "convex.hull")),
        
        div(style = paste0("height: ", input.height, "px;"),
            selectInput(inputId = "exclude.area",
                        label = "EOO exclude sea:",
                        choices = c(TRUE, FALSE),
                        selected = FALSE)),
        
        div(style = paste0("height: ", input.height, "px;"),
            sliderTextInput(inputId = "alpha",
                            label = "EOO alpha (radius):",
                            choices = c(seq(0, 1, by = 0.1), seq(2, 100, by = 1)),
                            selected = 2,
                            grid = FALSE)),
        
        div(style = paste0("height: ", input.height, "px;"),
            sliderTextInput(inputId = "buff.alpha",
                            label = "EOO buffer (degrees):",
                            choices = c(seq(0.001, 0.009, by = 0.001), seq(0.01, 0.09, by = 0.01), seq(0.1, 0.95, by = 0.05), seq(1, 10, by = 0.5)),
                            selected = 0.01,
                            grid = FALSE)),
        
        ## IUCN Cat B options (ConR)
        br(),
        h4("IUCN Category B:"),
        
        # AOO
        if (aooslider) {
          div(style = paste0("height: ", input.height, "px;"),
          sliderTextInput(inputId = "Cell_size_AOO",
                          label = "AOO cell size (km):",
                          choices = c(seq(0.1, 0.9, by = 0.1), seq(1, 50, by = 1)),
                          selected = 2,
                          grid = FALSE))
        },
        
        # subpopulations
        div(style = paste0("height: ", input.height, "px;"),
        sliderTextInput(inputId = "Resol_sub_pop",
                        label = "Subpopulation radius (km):",
                        choices = c(seq(0.1, 0.9, by = 0.1), seq(1, 50, by = 1)),
                        selected = 5,
                        grid = FALSE)),
        
        # locations
        div(style = paste0("height: ", input.height, "px;"),
            selectInput(inputId = "protec.areas",
                        label = "Protected areas:",
                        choices = c("NULL", "Madagascar.protec", "PA2", "PA123"),
                        selected = "PA2")),
        
        div(style = paste0("height: ", input.height, "px;"),
            selectInput(inputId = "method_protected_area",
                        label = "Location method:",
                        choices = c("no_more_than_one","other"),
                        selected = "no_more_than_one")),
        
        switch(method_locations,
               "sliding scale" = {
                 div(style = paste0("height: ", input.height, "px;"),
                 sliderTextInput(inputId = "Rel_cell_size",
                                 label = "Location multiplyer:",
                                 choices = c(seq(0.01, 1, by = 0.01)),
                                 selected = 0.05,
                                 grid = FALSE))
               },
               "fixed_grid" = {
                 div(style = paste0("height: ", input.height, "px;"),
                 sliderTextInput(inputId = "Cell_size_locations",
                                 label = "Location threat scale (km):",
                                 choices = c(seq(0.1, 0.9, by = 0.1), seq(1, 19, by = 1), seq(20, 95, by = 5), seq(100, 1000, by = 50)),
                                 selected = 10,
                                 grid = FALSE))
               }),

        # Input: Slider to highlight genus
        div(style = paste0("height: ", input.height, "px;"),
            selectInput(inputId = "genus",
                        label = "Sampling effort genus:",
                        choices = sort(names(sort(table(dpu$Genus), decreasing = TRUE))[1:1000]),
                        selected = "Dalbergia")),
        
        # bandwidth for sampling effort kernel density estimation
        div(style = paste0("height: ", input.height, "px;"),
        sliderTextInput(inputId = "sdl",
                        label = "Sampling density levels:",
                        choices = seq(1,20, by = 1),
                        selected = 10,
                        grid = FALSE)),
      
        ## Display options:
        h4("\nDisplay options:"),
        
        # Input: Slider for the extent buffer (area outside Polygon to display)
        switch(strsplit(project.raster, " ")[[1]][1],
               "+proj=utm" = {
                 div(style = paste0("height: ", input.height, "px;"),
                 sliderInput(inputId = "bufferx",
                             label = ifelse(buffer.method == "xy", "Raster x margin (m):", "Raster margin (m):"),
                             min = 0,
                             max = 1E06,
                             value = 20000,
                             step = 10000,
                             ticks = FALSE))
               },
               "+proj=longlat" = {
                 div(style = paste0("height: ", input.height, "px;"),
                 sliderInput(inputId = "bufferx",
                             label = ifelse(buffer.method == "xy", "Raster x margin (degrees):", "Raster margin (degrees):"),
                             min = 0,
                             max = 10,
                             value = 0.5,
                             step = 0.5,
                             ticks = FALSE))
               }),
        if (buffer.method == "xy") {
          switch(strsplit(project.raster, " ")[[1]][1],
                 "+proj=utm" = {
                   div(style = paste0("height: ", input.height, "px;"),
                   sliderInput(inputId = "buffery",
                               label = "Raster y margin (m):",
                               min = 0,
                               max = 1E06,
                               value = 20000,
                               step = 10000,
                               ticks = FALSE))
                 },
                 "+proj=longlat" = {
                   div(style = paste0("height: ", input.height, "px;"),
                   sliderInput(inputId = "buffery",
                               label = "Raster y margin (degrees):",
                               min = 0,
                               max = 10,
                               value = 0.5,
                               step = 0.5,
                               ticks = FALSE))
                 })
        },

        if (nb.cells) {
          # Input: Slider for the maximum number of pixels (maps)
          div(style = paste0("height: ", input.height, "px;"),
          sliderInput(inputId = "maxpixels",
                      label = "Number of Cells:",
                      min = 1e05,
                      max = 5e06,
                      value = 4.2e06, # ~default in addRasterImage
                      step = 1e04,
                      ticks = FALSE))
        },
                
        # Input: Slider for the color opacity (interactive map)
        div(style = paste0("height: ", input.height, "px;"),
        sliderInput(inputId = "opacity",
                    label = "Raster opacity:",
                    min = 0,
                    max = 1,
                    value = 0.5,
                    step = 0.1,
                    ticks = FALSE)),
        
        # Input: Altitude raster selection
        div(style = paste0("height: ", input.height, "px;"),
        selectInput(inputId = "altitude",
                    label = "Show altitude raster:",
                    choices = c("raster", "polygons", "neither (fast)"),
                    selected = "neither (fast)")),
        
        # Input: Slider for plot sizes
        div(style = paste0("height: ", input.height, "px;"),
        sliderInput(inputId = "psize",
                    label = "Plot size (pixels):",
                    min = 500,
                    max = 5000,
                    value = 800,
                    step = 100,
                    ticks = FALSE)),
        
        
        ## Download options
        h4("\nDownload options:"),
        
        # Input: Slider for the size of the downloaded plot
        div(style = paste0("height: ", input.height, "px;"),
        sliderInput(inputId = "size",
                    label = "Download size (in):",
                    min = 5,
                    max = 50,
                    value = 15,
                    step = 1,
                    ticks = FALSE)),
        
        # Input: Prefix of the downloaded files
        div(style = paste0("height: ", input.height*1.2, "px;"),
        textInput(inputId = "name",
                  label = "Download prefix:",
                  value = "MadaCAT")),
        
        # Export: export key variables from R environment to disk (.rda)
        div(style = paste0("height: ", input.height*0.75, "px;"),
            downloadButton("export", "Export", style = "width:100%;"))

      ),
      
      mainPanel(
        tabsetPanel(
          ## Own data
          tabPanel("Cat A + B (map)",
                   leafletOutput("leafletAB", height = plot.height),
                   downloadButton("down1", "Download")
          ),
          tabPanel("Cat A (static)",
                   downloadButton("down2", "Download"),
                   plotOutput("staticplot")
          ),
          tabPanel("Cat A (forecast)",
                   downloadButton("down3", "Download"),
                   plotlyOutput("forecast")
          ),
          tabPanel("Cat A (table)",
                   downloadButton("down4", "Download"),
                   DT::dataTableOutput("table1")
          ),
          tabPanel("Cat B (map)",
                   leafletOutput("leafletB", height = plot.height),
                   downloadButton("down5", "Download")
                   # https://stackoverflow.com/questions/38713809/how-to-dynamically-change-the-size-of-leaflet-map-in-shiny-r
          ),
          tabPanel("Cat B (table)",
                   downloadButton("down6", "Download"),
                   DT::dataTableOutput("table2")
          ),
          
          ## Vieilledent 2018 data
          tabPanel("Forest Cover (MAD)",
                   plotlyOutput("ForestCover")
          ),
          tabPanel("Deforestation (MAD)",
                   plotlyOutput("Deforestation")
          ),
          tabPanel("Deforestation Rate (MAD)",
                   plotlyOutput("DeforestationRate")
          )
        ),
      width = plot.width)
    )
  )
)


##############
### SERVER ###
##############
server <- function(input, output) {
  
  ## Get input parameters
  if (altrangeslider) {
    altmin <- reactive({input$altrange[1]})
    altmax <-  reactive({input$altrange[2]})
  } else {
    altmin <- reactive({as.numeric(input$altmin)})
    altmax <- reactive({as.numeric(input$altmax)})
  }
  fn <- renderText({input$name})
  psize <- reactive({input$psize})
  
  # ConR
  method.range <- reactive({input$method.range})
  exclude.area <- reactive({as.logical(input$exclude.area)})
  alpha <- reactive({as.numeric(input$alpha)})
  buff.alpha <- reactive({as.numeric(input$buff.alpha)})
  Cell_size_AOO <- reactive({ifelse(aooslider, as.numeric(input$Cell_size_AOO), 2)})
  Resol_sub_pop <- reactive({as.numeric(input$Resol_sub_pop)})
  protec.areas <- reactive({if (input$protec.areas == "NULL") return(NULL) else return(get(input$protec.areas))})
  ID_shape_PA <- reactive({names(protec.areas())[names(protec.areas()) %in% c("NOM_AP_2","WDPA_PID")][1]})
  method_protected_area <- reactive({input$method_protected_area})
  Cell_size_locations <- reactive({as.numeric(input$Cell_size_locations)})
  Rel_cell_size <- reactive({as.numeric(input$Rel_cell_size)})
  
  
  ## Get input coordinates (CSV)
  dd <- reactive({
    req(input$filedata)
    
    dd <- read.csv(input$filedata$datapath) # 1 = Latitude; 2 = Longitude, header expected. Optional: Darwin format columns with additional info
    return(get.DATA(x = dd, var.date = var.date, var.extant = var.extant, str.extant = str.extant))
  })
  
  
  ## Get historic EOO / AOO
  R1 <- reactive({
    req(dd())
    
    if (verbose) cat(paste0(Sys.time(), " computing historic EOO / AOO / Subpopulations / Locations...\n"))
    
    DATA <- dd()
    
    # make assessment
    R1 <- get.iucn(DATA = DATA[,c(1:4)], title = paste0("HISTORIC (N = ", nrow(DATA), ")"), get.aoo.grid = get.aoo.grid, get.loc.grid = get.loc.grid,
                   country_map = get(country_map), protec.areas = protec.areas(), ID_shape_PA = ID_shape_PA(),
                   method.range = method.range(), exclude.area = exclude.area(), alpha = alpha(), buff.alpha = buff.alpha(),
                   Cell_size_AOO = Cell_size_AOO(), nbe.rep.rast.AOO = nbe.rep.rast.AOO,
                   Resol_sub_pop = Resol_sub_pop(),
                   method_protected_area = method_protected_area(), method_locations = method_locations,
                   Cell_size_locations = ifelse(method_locations == "fixed_grid", Cell_size_locations(), 10), 
                   Rel_cell_size = ifelse(method_locations == "sliding scale", Rel_cell_size(), 0.05),
                   DrawMap = DrawMap, add.legend = add.legend, file_name = file_name, map_pdf = map_pdf,
                   showWarnings = showWarnings, export_shp = export_shp, write_shp = write_shp,
                   write_results = write_results, write_file_option = write_file_option)
    
    if (verbose) {
      # print(data.frame(t(R1$Results)))
      cat(paste0(Sys.time(), " done.\n"))
    }
    return(R1)
  })
  
  
  ## Get extant EOO / AOO
  R2 <- reactive({
    req(dd())
    
    if (verbose) cat(paste0(Sys.time(), " computing extant EOO / AOO / Subpopulations / Locations...\n"))
    
    DATA <- subset(dd(), occurrenceRemarks %in% str.extant)
    
    # make assessment
    R2 <- get.iucn(DATA = DATA[,c(1:4)], title = paste0("EXTANT (N = ", nrow(DATA), ")"), get.aoo.grid = get.aoo.grid, get.loc.grid = get.loc.grid,
                   country_map = get(country_map), protec.areas = protec.areas(), ID_shape_PA = ID_shape_PA(),
                   method.range = method.range(), exclude.area = exclude.area(), alpha = alpha(), buff.alpha = buff.alpha(),
                   Cell_size_AOO = Cell_size_AOO(), nbe.rep.rast.AOO = nbe.rep.rast.AOO,
                   Resol_sub_pop = Resol_sub_pop(),
                   method_protected_area = method_protected_area(), method_locations = method_locations,
                   Cell_size_locations = ifelse(method_locations == "fixed_grid", Cell_size_locations(), 10), 
                   Rel_cell_size = ifelse(method_locations == "sliding scale", Rel_cell_size(), 0.05),
                   DrawMap = DrawMap, add.legend = add.legend, file_name = file_name, map_pdf = map_pdf,
                   showWarnings = showWarnings, export_shp = export_shp, write_shp = write_shp,
                   write_results = write_results, write_file_option = write_file_option)
    
    if (verbose) {
      # print(data.frame(t(R2$Results)))
      cat(paste0(Sys.time(), " done.\n"))
    }
    return(R2)
  })

  
  ## Get projected collection polygon
  S <- reactive({
    req(dd()) ; req(R1())
    
    #   S <- SpatialPolygons(list(Polygons(list(Polygon(dd()[chull(dd()),c(2,1)])),1)))
    #   attributes(S)$proj4string <- crs(project.coords)
    #   
    #   # project to raster
    #   if (project.coords != project.raster) {
    #     S <- spTransform(S, CRS(project.raster))
    #   }
    S <- spTransform(R1()$spatialPoly_EOO, CRS(project.raster))
    return(S)
  })
  
  
  ## Get projected collection points
  s <- reactive({
    req(dd())
    
    s <- list(longlat = SpatialPoints(dd()[c(2,1)]))
    attributes(s[["longlat"]])$proj4string <-  crs(project.coords)
    
    # project to raster
    if (project.coords != project.raster) {
      s[["utm"]] <- spTransform(s[["longlat"]], CRS(project.raster))
    }
    return(s)
  })

  
  ## Crop (forest-projected) altitude by polygon
  alt <- reactive({
    req(input$filedata) ; req(dd()) ; req(e()) ; req(R1()) ; req(R2())
    
    if (verbose) cat(paste0(Sys.time(), " cropping altitude raster...\n"))
    alt <- crop(altitude, e()$utm)
    attributes(alt)$proj4string <- attributes(get(difvars[1]))$crs
    return(alt)
  })

  
  ## Define extent to display habitat
  e <- reactive({
    req(S())
    
    dX <- extent(S())@xmax - extent(S())@xmin
    dY <- extent(S())@ymax - extent(S())@ymin
    
    e <- list()
    
    switch(buffer.method,
           x = {
             bufferx <- input$bufferx
             takeX <- ifelse(dX > dY, TRUE, FALSE)
             e[["utm"]] <- extent(c(ifelse(takeX, extent(S())@xmin - bufferx, extent(S())@xmax - dX/2 - dY/2 - bufferx),
                                    ifelse(takeX, extent(S())@xmax + bufferx, extent(S())@xmax - dX/2 + dY/2 + bufferx),
                                    ifelse(takeX, extent(S())@ymax - dY/2 - dX/2 - bufferx, extent(S())@ymin - bufferx),
                                    ifelse(takeX, extent(S())@ymax - dY/2 + dX/2 + bufferx, extent(S())@ymax + bufferx)))
             x <- SpatialPoints(data.frame(lon = c(e[["utm"]]@xmin, e[["utm"]]@xmax, e[["utm"]]@xmin, e[["utm"]]@xmax), 
                                           lat = c(e[["utm"]]@ymin, e[["utm"]]@ymin, e[["utm"]]@ymax, e[["utm"]]@ymax)))
             attributes(x)$proj4string <-  crs(project.raster)
             x <- spTransform(x, CRS(project.coords))
             e[["longlat"]] <- extent(x)
           },
           xy = {
             bufferx <- input$bufferx
             buffery <- input$buffery
             e[["utm"]] <- extent(c(extent(S())@xmin - bufferx,
                                    extent(S())@xmax + bufferx,
                                    extent(S())@ymin - buffery),
                                    extent(S())@ymax + buffery)
             
             x <- SpatialPoints(data.frame(lon = c(e[["utm"]]@xmin, e[["utm"]]@xmax, e[["utm"]]@xmin, e[["utm"]]@xmax), 
                                           lat = c(e[["utm"]]@ymin, e[["utm"]]@ymin, e[["utm"]]@ymax, e[["utm"]]@ymax)))
             attributes(x)$proj4string <-  crs(project.raster)
             x <- spTransform(x, CRS(project.coords))
             e[["longlat"]] <- extent(x)
           })

    return(e)
  })
  
  
  ## Sampling effort data
  ds <- reactive({
    req(e())
    
    setkey(dpu, LatitudeDecimal, LongitudeDecimal)
    sel <- which(dpu$LatitudeDecimal %between% c(e()[["longlat"]]@ymin, e()[["longlat"]]@ymax) & 
                   dpu$LongitudeDecimal %between% c(e()[["longlat"]]@xmin, e()[["longlat"]]@xmax))

    if (length(sel) > 0) {
      
      # SpatialPoints of jittered sampling effort coordinates
      dsub <- dpu[sel, nomatch=0L] 
      dsub[, jitx := jitter(LongitudeDecimal, factor = 0, amount = jitter.amount)]
      dsub[, jity := jitter(LatitudeDecimal, factor = 0, amount = jitter.amount)]
      x <- SpatialPoints(data.frame(dsub[,c("jitx","jity")]))
      attributes(x)$proj4string <-  crs(project.coords)
      
      # year distribution
      dyr <- sapply(strsplit(sapply(strsplit(pu[sel], split = "<br/>Minimum.Date:"), "[", 2), split = "<br/>Habitat:"), "[", 1)
      dyr <- as.numeric(sapply(strsplit(dyr, split = "[/]"), "[", 3))
      dyr[dyr < 1700] <- NA # discard likely errors
      
      # subset genus
      dsub$highlight <- 0
      dsub$highlight[which(dsub$Genus == input$genus)] <- 1
      x <- SpatialPointsDataFrame(x, data = data.frame(highlight = dsub$highlight))
    } else {
      x <- NA
      dyr <- NA
    }
    
    return(list(sp = x, popups = pu[sel], years = dyr))
  })
  
  
  ## Get habitat difference raster
  fcc <- reactive({
    req(dd())
    
    obj <- paste0("for", input$period1, "_vs_for", input$period2)
    get(obj)
  })
  
  
  ## Compute habitat change between periods
  r <- reactive({
    req(alt()) ; req(fcc()) ; req(S())
    
    if (verbose) cat(paste0(Sys.time(), " computing habitat change between ", input$period1, " and ", input$period2, "...\n"))
    habitat.stats(x = fcc(), S = S(), alt = alt(), 
                  altmin = altmin(), altmax = altmax(), e = e()$utm, fcc = TRUE)
  })
  
  d.for1953 <- reactive({req(r()) ; habitat.wrapper(year = "1953", r = r(), S = S(), alt = alt(), altmin = altmin(), altmax = altmax(), e = e()$utm, period1 = input$period1, period2 = input$period2, verbose = verbose)})
  d.for1973 <- reactive({req(r()) ; habitat.wrapper(year = "1973", r = r(), S = S(), alt = alt(), altmin = altmin(), altmax = altmax(), e = e()$utm, period1 = input$period1, period2 = input$period2, verbose = verbose)})
  d.for1990 <- reactive({req(r()) ; habitat.wrapper(year = "1990", r = r(), S = S(), alt = alt(), altmin = altmin(), altmax = altmax(), e = e()$utm, period1 = input$period1, period2 = input$period2, verbose = verbose)})
  d.for2000 <- reactive({req(r()) ; habitat.wrapper(year = "2000", r = r(), S = S(), alt = alt(), altmin = altmin(), altmax = altmax(), e = e()$utm, period1 = input$period1, period2 = input$period2, verbose = verbose)})
  d.for2005 <- reactive({req(r()) ; habitat.wrapper(year = "2005", r = r(), S = S(), alt = alt(), altmin = altmin(), altmax = altmax(), e = e()$utm, period1 = input$period1, period2 = input$period2, verbose = verbose)})
  d.for2010 <- reactive({req(r()) ; habitat.wrapper(year = "2010", r = r(), S = S(), alt = alt(), altmin = altmin(), altmax = altmax(), e = e()$utm, period1 = input$period1, period2 = input$period2, verbose = verbose)})
  d.for2014 <- reactive({req(r()) ; habitat.wrapper(year = "2014", r = r(), S = S(), alt = alt(), altmin = altmin(), altmax = altmax(), e = e()$utm, period1 = input$period1, period2 = input$period2, verbose = verbose)})
  d.for2015 <- reactive({req(r()) ; habitat.wrapper(year = "2015", r = r(), S = S(), alt = alt(), altmin = altmin(), altmax = altmax(), e = e()$utm, period1 = input$period1, period2 = input$period2, verbose = verbose)})
  d.for2017 <- reactive({req(r()) ; habitat.wrapper(year = "2017", r = r(), S = S(), alt = alt(), altmin = altmin(), altmax = altmax(), e = e()$utm, period1 = input$period1, period2 = input$period2, verbose = verbose)})

  ## Project habitat change raster for leaflet
  rp <- reactive({
    req(r())
    
    if (projectForLeaflet) {
      if (verbose) cat(paste0(Sys.time(), " projecting habitat change raster for leaflet...\n"))
      
      r <- projectRasterForLeaflet(r()$raster, method = "ngb")
    } else {
      r <- r()
    }
    return(r)
  })
  
  ## Project habitat rasters for leaflet
  hp <- reactive({
    req(alt())
    
    if (projectForLeaflet & input$altitude == "raster") {
      if (verbose) cat(paste0(Sys.time(), " projecting altitude raster for leaflet...\n"))
      alt <- projectRasterForLeaflet(alt(), method = "ngb")
    } else {
      alt <- alt()
    }
    return(alt)
  })

  
  ### PLOTS 
  ## Interactive leaflet map (IUCN A and B)
  plotInput1 <- eventReactive(input$submit, {
    req(rp()) ; req(hp()) ; req(dd()) ; req(R1()) ; req(R2())
    
    # get extant
    dd <- dd()
    pal.ext <- colorFactor(palette = c(col.historic, col.extant), domain = dd[,"extant"])
    
    # jitter coordinates
    dd.jitter <- dd
    set.seed(4065)
    dd.jitter[,2] <- jitter(dd[,2], factor = 0, amount = jitter.amount)
    set.seed(4165)
    dd.jitter[,1] <- jitter(dd[,1], factor = 0, amount = jitter.amount)
    dd.jitter <- SpatialPoints(dd.jitter[,c(2,1)])
    
    # define raster category colors and labels
    r <- rp()
    alt <- hp()
    
    cpal <- c("red", "lightgray", "white", "lightgreen", "darkgreen")
    at <- c(-2,-1,0,1,2)
    lab <- c("forest loss (in habitat)",     # -2
             "forest (gain/loss/constant, outside habitat)",     # -1
             "no forest",                    #  0
             "forest constant (in habitat)", #  1
             "forest gain (in habitat)"      #  2
    )
    
    # drop unused color levels
    choose <- sapply(at, function(x) {any(x %in% values(r))})
    cpal <- cpal[choose]
    lab <- lab[choose]
    
    if (verbose) cat(paste0(Sys.time(), " creating leaflet map...\n"))
    
    # decrease raster resolution for display
    # nc <- ncell(r)
    # if (nc > l.maxpixels) {
    #   if (verbose) cat(paste0(Sys.time(), " aggregating raster...\n"))
    #   # aggregate raster number of cells exceeds threshold (TO DO)
    #   r.maxpixels <- 5E06
    #   ratio <- ncell(r)/r.maxpixels
    #   r <- raster::aggregate(x = r, fact = f, fun = modal, ties = "highest")
    # }
    
    ## create leaflet map
    tit <- paste0("Change in habitat [", input$period1, " - ", input$period2, 
                  ", at ", altmin(), " - ", altmax(), "m a.s.l.]")
    
    # basic leaflet tiles
    M <- leaflet()
    for (i in 1:length(providers)) {
      M <- M %>% 
        addProviderTiles(providers[i], group = providers[i])
    }
    
    # bounds
    coords <- na.omit(dd[,c(1,2)])
    buffer <- 1
    M <- M  %>%
      fitBounds(min(coords[,2])*buffer, min(coords[,1])*buffer, max(coords[,2])*buffer, max(coords[,1])*buffer)
    
    # add altitude raster
    overlayGroups <- c("Habitat")
    hideGroups <- character()
    
    if (input$altitude == "raster") {
   
      from <- min(values(alt), na.rm = TRUE)
      to <-  max(values(alt), na.rm = TRUE)
      palette.range <- palette.breaks[(max(which(palette.breaks <= from))):(max(which(palette.breaks <= to)))]
      pal.alt <- palette.altitude[palette.breaks %in% palette.range]
      pal.alt <- colorNumeric(palette = pal.alt, domain = values(alt))
      
      M <- M %>%
        addRasterImage(alt, 
                       project = FALSE,
                       colors = pal.alt, 
                       opacity = input$opacity,
                       maxBytes = ifelse(nb.cells, input$maxpixels, l.maxpixels),
                       group = "Altitude") %>%
        addLegend(pal = pal.alt,
                  opacity = input$opacity,
                  values = values(alt),
                  title = "Altitude (m)",
                  group = "Altitude")
      
      overlayGroups <- c(overlayGroups, "Altitude")
      hideGroups <- c(hideGroups, "Altitude")
    }
    
    if (input$altitude == "polygons") {
      
      if (verbose) cat(paste0(Sys.time(), " computing smoothed altitude polygons...\n"))
      from <- min(values(alt), na.rm = TRUE)
      to <-  max(values(alt), na.rm = TRUE)
      palette.range <- palette.breaks[(max(which(palette.breaks <= from))):(max(which(palette.breaks <= to)))]
      pal.alt <- palette.altitude[palette.breaks %in% palette.range]
      pal.alt <- colorFactor(palette = pal.alt, domain = palette.range)
      
      spgons <- raster2SPDF(r = alt, breaks = palette.range, method = "continuous",
                            drop.crumbs = TRUE, thr = units::set_units(0.5, "km^2"), 
                            project = TRUE, proj4string = project.coords)

      popup.alt <- paste(paste(c(from, palette.range[-1]), c(palette.range[-1], to), sep = " - "), "m")
      
      M <- M %>%
        addPolygons(data = spgons, 
                    color = pal.alt(spgons$alt),
                    opacity = input$opacity,
                    weight = 1,
                    fillColor = pal.alt(spgons$alt),
                    fillOpacity = input$opacity,
                    popup = popup.alt,
                    group = "Altitude")
      
      overlayGroups <- c(overlayGroups, "Altitude")
      hideGroups <- c(hideGroups, "Altitude")
    }
    
    # add ecological layers
    if (add.layers) {
      pal.eco <- colorFactor(palette = c("yellow","purple","lightgreen","pink","red","darkgreen","orange"), domain = eco2017$ECO_NAME)
      pal.vege <- colorFactor(palette = c("#2D82AF","#2D82AF","#F16667","orange","#6F9E4C","#F16667","darkgreen","#825D99","pink","#B89B74","#F0EB99"), domain = vege$VEGETATION)
      
      M <- M %>%
        addPolygons(data = eco2017,
                    color = ~pal.eco(eco2017$ECO_NAME),
                    weight = 2,
                    opacity = input$opacity,
                    popup = get.popup(eco2017@data, c("ECO_NAME", "BIOME_NAME", "REALM", "NNH_NAME")),
                    group = "Ecoregions") %>%
        addPolygons(data = vege,
                    color = ~pal.vege(vege$VEGETATION),
                    weight = 2,
                    opacity = input$opacity,
                    popup = get.popup(vege@data, c("VEGETATION", "GEOLOGY_NA", "ACRES", "HA")),
                    group = "Primary Vegetation (1996)") %>%
        addPolygons(data = PA1,
                    color = "orange",
                    opacity = 1,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    group = "Protected areas") %>% # popups are similar to PA2
        addPolygons(data = PA2,
                    color = "orange",
                    opacity = 1,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = get.popup(PA2@data, names(PA2)),
                    group = "Protected areas") %>%
        addPolygons(data = PA3,
                    color = "orange",
                    opacity = 1,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = get.popup(PA3@data, c("DATAADMIN","NOM","HECTARES","DATE_CREAT","STATUT_A_1","SHORT_NAME","PROVINCE","REGION","DISTRICT")),
                    group = "Protected areas")
      
      overlayGroups <- c(overlayGroups, "Ecoregions", "Primary Vegetation (1996)", "Protected areas")
      hideGroups <- c(hideGroups, "Ecoregions", "Primary Vegetation (1996)", "Protected areas")
    }
    
    if (input$protec.areas != "NULL") {
      M <- M %>%
        addPolygons(data = PA1,
                    color = "orange",
                    opacity = 1,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = get.popup(protec.areas()@data, names(protec.areas())),
                    group = input$protec.areas)
      
      overlayGroups <- c(overlayGroups, input$protec.areas)
      # hideGroups <- c(hideGroups, input$protec.areas)
        
    }
    
    # add habitat, EOO, AOO, Subpopulations and (jittered) occurrences
    M <- M %>%
      addRasterImage(r,
                     project = FALSE,
                     colors = cpal, 
                     opacity = input$opacity,
                     maxBytes = ifelse(nb.cells, input$maxpixels, l.maxpixels),
                     group = "Habitat") %>%
      addLegend(colors = cpal,
                opacity = input$opacity,
                labels = lab,
                values = values(r),
                title = tit,
                group = "Habitat") %>%
      addPolygons(data = R1()$spatialPoly_EOO,
                  color = col.historic,
                  opacity = input$opacity,
                  fillColor = NA,
                  fillOpacity = 0,
                  weight = 5,
                  popup = R1()$popup.eoo,
                  group = "EOO (historic)") %>%
      addPolygons(data = R2()$spatialPoly_EOO,
                  color = col.extant,
                  opacity = input$opacity,
                  fillColor = NA,
                  fillOpacity = 0,
                  weight = 5,
                  popup = R2()$popup.eoo,
                  group = "EOO (extant)")
    
    overlayGroups <- c(overlayGroups, "EOO (historic)", "EOO (extant)")
    hideGroups <- c(hideGroups)

    # add AOO grid, subpopulation polygon and location grid
    if (get.aoo.grid) {
      M <- M %>%
        addPolygons(data = R1()$spatialPoly_AOO,
                    color = col.historic,
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = NULL,
                    group = "AOO (historic)") %>%
        addPolygons(data = R2()$spatialPoly_AOO,
                    color = col.extant,
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = NULL,
                    group = "AOO (extant)")
      
      overlayGroups <- c(overlayGroups, "AOO (historic)", "AOO (extant)")
      hideGroups <- c(hideGroups, "AOO (historic)", "AOO (extant)")
    }
    
    if (get.subpop) {
      M <- M  %>%
        addPolygons(data = R1()$spatialPoly_subpop,
                    color = col.historic,
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = R1()$popup.sp,
                    group = "Subpopulations (historic)") %>%
        addPolygons(data = R2()$spatialPoly_subpop,
                    color = col.extant,
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = R2()$popup.sp,
                    group = "Subpopulations (extant)")
      
      overlayGroups <- c(overlayGroups, "Subpopulations (historic)", "Subpopulations (extant)")
      hideGroups <- c(hideGroups, "Subpopulations (historic)")
    }
   
    if (get.loc.grid) {
      M <- M %>%
        addPolygons(data = R1()$spatialPoly_Loc,
                    color = col.historic,
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = NULL,
                    group = "Locations (historic)") %>%
        addPolygons(data = R2()$spatialPoly_Loc,
                    color = col.extant,
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = NULL,
                    group = "Locations (extant)")
      
      overlayGroups <- c(overlayGroups, "Locations (historic)", "Locations (extant)")
      hideGroups <- c(hideGroups, "Locations (historic)", "Locations (extant)")
    }
    
    # add sampling effort
    if (add.sampling.effort) {
      
      if (class(ds()$sp) != "logical") {
        # get contour lines, contour polygons and count table popups
        dcontour <- get.contour(sp = ds()$sp, margin = contour.margin, nlevels = input$sdl, sqrt = contour.sqrt)
        cpgons <- contour.polygons(CL = dcontour$CL, sp = ds()$sp, crs = project.coords)
        popups <- contour.popups(x = ds()$years, mat = cpgons$pinpgon, 
                                 b = c(1699, 1799, 1899, seq(1909, 1999, by = 10), seq(2000, now, by = 1)), 
                                 thr = 1999, n = cpgons$nperpgon)
        pal.effort <- colorFactor(palette = c("gray75","orange"), domain = ds()$sp$highlight)
        
        M <- M %>%
          addPolygons(data = cpgons$spgons,
                      color = heat.colors(cpgons$NLEV, NULL)[cpgons$LEVS],
                      fillOpacity = 0.1,
                      opacity = input$opacity,
                      weight = 1,
                      group = "Sampling Density",
                      popup = popups) %>%
        addCircleMarkers(data = ds()$sp,
                         radius = ~ 1.5*(ds()$sp$highlight+1)^2,
                         color = "black",
                         fillColor = pal.effort(ds()$sp$highlight),
                         fillOpacity = 0.8,
                         weight = 1,
                         popup = ds()$popups,
                         group = "Sampling Effort")
        
        overlayGroups <- c(overlayGroups, "Sampling Density", "Sampling Effort")
        hideGroups <- c(hideGroups, "Sampling Density", "Sampling Effort")
      }
    }
    
    # add layers control and observations
    overlayGroups <- c(overlayGroups, "Observations")
    hideGroups <- c(hideGroups)
    
    M <- M %>%
      addCircleMarkers(data = dd.jitter, 
                       radius = 6,
                       color = "black",
                       fillColor = pal.ext(dd[,"extant"]),
                       fillOpacity = 0.8, # input$opacity,
                       weight = 1,
                       popup = get.popups(dd),
                       group = "Observations") %>%
      addScaleBar(position = scale.position, options = scaleBarOptions(maxWidth = 200, imperial = FALSE)) %>%
      addLayersControl(baseGroups = providers,
                       overlayGroups = overlayGroups,
                       position = "topleft") %>%
      hideGroup(hideGroups) %>%
      # addDrawToolbar(editOptions = editToolbarOptions(selectedPathOptions = selectedPathOptions())) %>% 
      addMeasure(position = "topright", 
                 primaryLengthUnit = "kilometers", secondaryLengthUnit = NULL,
                 primaryAreaUnit = "hectares", secondaryAreaUnit = NULL)
    
    if (verbose) cat(paste0(Sys.time(), " done.\n"))
    return(M)
  })
  output$leafletAB <- renderLeaflet({plotInput1()})
  
  
  ## Interactive leaflet map (IUCN B only)
  plotInput5 <- eventReactive(input$submit, {
    req(dd()) ; req(R1()) ; req(R2())
    
    # get extant
    dd <- dd()
    pal.ext <- colorFactor(palette = c(col.historic, col.extant), domain = dd[,"extant"])
    
    # jitter coordinates
    dd.jitter <- dd
    set.seed(4065)
    dd.jitter[,2] <- jitter(dd[,2], factor = 0, amount = jitter.amount)
    set.seed(4165)
    dd.jitter[,1] <- jitter(dd[,1], factor = 0, amount = jitter.amount)
    dd.jitter <- SpatialPoints(dd.jitter[,c(2,1)])
    
    if (verbose) cat(paste0(Sys.time(), " creating leaflet map...\n"))
    
    ## create leaflet map
    # basic leaflet tiles
    M <- leaflet()
    for (i in 1:length(providers)) {
      M <- M %>% 
        addProviderTiles(providers[i], group = providers[i])
    }
    
    # bounds
    coords <- na.omit(dd[,c(1,2)])
    buffer <- 1
    M <- M  %>%
      fitBounds(min(coords[,2])*buffer, min(coords[,1])*buffer, max(coords[,2])*buffer, max(coords[,1])*buffer)
    
    # add altitude raster
    overlayGroups <- character()
    hideGroups <- character()
    
    # add ecological layers
    if (add.layers) {
      pal.eco <- colorFactor(palette = c("yellow","purple","lightgreen","pink","red","darkgreen","orange"), domain = eco2017$ECO_NAME)
      pal.vege <- colorFactor(palette = c("#2D82AF","#2D82AF","#F16667","orange","#6F9E4C","#F16667","darkgreen","#825D99","pink","#B89B74","#F0EB99"), domain = vege$VEGETATION)
      
      M <- M %>%
        addPolygons(data = eco2017,
                    color = ~pal.eco(eco2017$ECO_NAME),
                    weight = 2,
                    opacity = input$opacity,
                    popup = get.popup(eco2017@data, c("ECO_NAME", "BIOME_NAME", "REALM", "NNH_NAME")),
                    group = "Ecoregions") %>%
        addPolygons(data = vege,
                    color = ~pal.vege(vege$VEGETATION),
                    weight = 2,
                    opacity = input$opacity,
                    popup = get.popup(vege@data, c("VEGETATION", "GEOLOGY_NA", "ACRES", "HA")),
                    group = "Primary Vegetation (1996)") %>%
        addPolygons(data = PA1,
                    color = "orange",
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    group = "Protected areas") %>% # popups are similar to PA2
        addPolygons(data = PA2,
                    color = "orange",
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = get.popup(PA2@data, names(PA2)),
                    group = "Protected areas") %>%
        addPolygons(data = PA3,
                    color = "orange",
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = get.popup(PA3@data, c("DATAADMIN","NOM","HECTARES","DATE_CREAT","STATUT_A_1","SHORT_NAME","PROVINCE","REGION","DISTRICT")),
                    group = "Protected areas")
      
      overlayGroups <- c(overlayGroups, "Ecoregions", "Primary Vegetation (1996)", "Protected areas")
      hideGroups <- c(hideGroups, "Ecoregions", "Primary Vegetation (1996)", "Protected areas")
    }
    
    if (input$protec.areas != "NULL") {
      M <- M %>%
        addPolygons(data = PA1,
                    color = "orange",
                    opacity = 1,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = get.popup(protec.areas()@data, names(protec.areas())),
                    group = input$protec.areas)
      
      overlayGroups <- c(overlayGroups, input$protec.areas)
      # hideGroups <- c(hideGroups, input$protec.areas)
      
    }
    
    # add EOO
    M <- M %>%
      addPolygons(data = R1()$spatialPoly_EOO,
                  color = col.historic,
                  opacity = input$opacity,
                  fillColor = NA,
                  fillOpacity = 0,
                  weight = 5,
                  popup = R1()$popup.eoo,
                  group = "EOO (historic)") %>%
      addPolygons(data = R2()$spatialPoly_EOO,
                  color = col.extant,
                  opacity = input$opacity,
                  fillColor = NA,
                  fillOpacity = 0,
                  weight = 5,
                  popup = R2()$popup.eoo,
                  group = "EOO (extant)")
    
    overlayGroups <- c(overlayGroups, "EOO (historic)", "EOO (extant)")
    hideGroups <- c(hideGroups)
    
    # add AOO grid, subpopulation polygon and location grid
    if (get.aoo.grid) {
      M <- M %>%
        addPolygons(data = R1()$spatialPoly_AOO,
                    color = col.historic,
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = NULL,
                    group = "AOO (historic)") %>%
        addPolygons(data = R2()$spatialPoly_AOO,
                    color = col.extant,
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = NULL,
                    group = "AOO (extant)")
      
      overlayGroups <- c(overlayGroups, "AOO (historic)", "AOO (extant)")
      hideGroups <- c(hideGroups, "AOO (historic)", "AOO (extant)")
    }
    
    if (get.subpop) {
      M <- M  %>%
        addPolygons(data = R1()$spatialPoly_subpop,
                    color = col.historic,
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = R1()$popup.sp,
                    group = "Subpopulations (historic)") %>%
        addPolygons(data = R2()$spatialPoly_subpop,
                    color = col.extant,
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = R2()$popup.sp,
                    group = "Subpopulations (extant)")
      
      overlayGroups <- c(overlayGroups, "Subpopulations (historic)", "Subpopulations (extant)")
      hideGroups <- c(hideGroups, "Subpopulations (historic)")
    }
    
    if (get.loc.grid) {
      M <- M %>%
        addPolygons(data = R1()$spatialPoly_Loc,
                    color = col.historic,
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = NULL,
                    group = "Locations (historic)") %>%
        addPolygons(data = R2()$spatialPoly_Loc,
                    color = col.extant,
                    opacity = input$opacity,
                    fillColor = NA,
                    fillOpacity = 0,
                    weight = 2,
                    popup = NULL,
                    group = "Locations (extant)")
      
      overlayGroups <- c(overlayGroups, "Locations (historic)", "Locations (extant)")
      hideGroups <- c(hideGroups, "Locations (historic)", "Locations (extant)")
    }
    
    # add sampling effort
    if (add.sampling.effort) {
      
      if (class(ds()$sp) != "logical") {
        # get contour lines, contour polygons and count table popups
        dcontour <- get.contour(sp = ds()$sp, margin = contour.margin, nlevels = input$sdl, sqrt = contour.sqrt)
        cpgons <- contour.polygons(CL = dcontour$CL, sp = ds()$sp, crs = project.coords)
        popups <- contour.popups(x = ds()$years, mat = cpgons$pinpgon, 
                                 b = c(1699, 1799, 1899, seq(1909, 1999, by = 10), seq(2000, now, by = 1)), 
                                 thr = 1999, n = cpgons$nperpgon)
        pal.effort <- colorFactor(palette = c("gray75","orange"), domain = ds()$sp$highlight)
        
        M <- M %>%
          addPolygons(data = cpgons$spgons,
                      color = heat.colors(cpgons$NLEV, NULL)[cpgons$LEVS],
                      fillOpacity = 0.1,
                      opacity = input$opacity,
                      weight = 1,
                      group = "Sampling Density",
                      popup = popups) %>%
          addCircleMarkers(data = ds()$sp,
                           radius = ~ 1.5*(ds()$sp$highlight+1)^2,
                           color = "black",
                           fillColor = pal.effort(ds()$sp$highlight),
                           fillOpacity = 0.8,
                           weight = 1,
                           popup = ds()$popups,
                           group = "Sampling Effort")
        
        overlayGroups <- c(overlayGroups, "Sampling Density", "Sampling Effort")
        hideGroups <- c(hideGroups, "Sampling Density", "Sampling Effort")
      }
    }
    
    
    # add layers control and observations
    overlayGroups <- c(overlayGroups, "Observations")
    hideGroups <- c(hideGroups)
    
    M <- M %>%
      addCircleMarkers(data = dd.jitter, 
                       radius = 6,
                       color = "black",
                       fillColor = pal.ext(dd[,"extant"]),
                       fillOpacity = 0.8, # input$opacity,
                       weight = 1,
                       popup = get.popups(dd),
                       group = "Observations") %>%
      addScaleBar(position = scale.position, options = scaleBarOptions(maxWidth = 200, imperial = FALSE)) %>%
      addLayersControl(baseGroups = providers,
                       overlayGroups = overlayGroups,
                       position = "topleft") %>%
      hideGroup(hideGroups) %>%
      # addDrawToolbar(editOptions = editToolbarOptions(selectedPathOptions = selectedPathOptions())) %>% 
      addMeasure(position = "topright", 
                 primaryLengthUnit = "kilometers", secondaryLengthUnit = NULL,
                 primaryAreaUnit = "hectares", secondaryAreaUnit = NULL)
    
    if (verbose) cat(paste0(Sys.time(), " done.\n"))
    return(M)
  })
  output$leafletB <- renderLeaflet({plotInput5()})
  
  
  ## Static lattice map
  plotInput2 <- eventReactive(input$submit, {
    req(r()) ; req(alt()) ; req(dd()) ; req(R1()) ; req(R2())
    
    if (verbose) cat(paste0(Sys.time(), " creating static (lattice) map...\n"))
    
    # get extant
    dd <- dd()
    pal.ext <- colorFactor(palette = c(col.historic, col.extant), domain = dd[,"extant"])
    
    # ratify
    r <- ratify(r()$raster)
    
    # define raster category colors and labels
    cpal <- c("red", "lightgray", "white", "lightgreen", "darkgreen")
    at <- c(-2,-1,0,1,2)
    lab <- c("forest loss (in habitat)",     # -2
             "forest (gain/loss/constant, outside habitat)",     # -1
             "no forest",                    #  0
             "forest constant (in habitat)", #  1
             "forest gain (in habitat)"      #  2
             )
    
    # drop unused color levels
    choose <- sapply(at, function(x) {any(x %in% values(r))})
    cpal <- cpal[choose]
    lab <- lab[choose]
    
    # plot
    tit <- paste0("Change in habitat [", input$period1, " - ", input$period2, 
                  ", at ", altmin(), " - ", altmax(), "m a.s.l.]")
    t <- identity(S())
    cols <- pal.ext(dd[,"extant"])
    
    P <- rasterVis::levelplot(r, col.regions = cpal, att = "ID", sub = tit,
                         maxpixels = ifelse(nb.cells, input$maxpixels, p.maxpixels), 
                         colorkey = list(col = cpal, labels = list(labels = lab))) +
        latticeExtra::layer(sp.points(s, cex = 1, pch = 23, col = col, alpha = 1), data = list(s = s()$utm, col = cols)) +
        latticeExtra::layer(sp.polygons(S, col = col, lwd = 2, alpha = 1), data = list(S = t, col = col.historic))
    
    if (verbose) cat(paste0(Sys.time(), " done.\n"))
    return(P)
  })
  output$staticplot <- renderPlot({plotInput2()}, height = psize, width = psize)
  
 
  ## Compute habitat change (in polygon and altitude range) time series
  d.habitat <- eventReactive(input$submit, {
    req(dd())
    
    if (verbose) cat(paste0(Sys.time(), " computing habitat time series for selected periods...\n"))
    
    d.res <- sapply(paste0("d.for", input$periods), FUN = function(x) {get(x)()})
    d.habitat <- get.stats(d.res, factor = fac, unit = unit, 
                           period1 = input$period1, period2 = input$period2,
                           refyear = paste0("d.for", input$period1), 
                           extent = paste0("EOO at ", altmin(), " - ", altmax(), " m"))
    
    # add IUCN category
    now <- as.numeric(unlist(strsplit(as.character(Sys.time()), split = "[-]"))[1])
    
    input <- list(period1 = 1953, period2 = 2017)
    a <- ifelse(input$period1 < now & input$period2 <= now, "A2",
                ifelse(input$period1 >= now & input$period2 > now, "A3",
                       ifelse(input$period1 < now & input$period2 > now, "A4", "A?")))
    
    d.habitat[,"IUCN Cat A1"] <- ifelse(d.habitat[,"Reduction_perc"] >= 90, "CR",
                                        ifelse(d.habitat[,"Reduction_perc"] >= 70, "EN",
                                               ifelse(d.habitat[,"Reduction_perc"] >= 50, "VU", "LC or NT")))
    d.habitat[,paste("IUCN Cat", a)] <- ifelse(d.habitat[,"Reduction_perc"] >= 80, "CR",
                                               ifelse(d.habitat[,"Reduction_perc"] >= 50, "EN",
                                                      ifelse(d.habitat[,"Reduction_perc"] >= 30, "VU", "LC or NT")))
    d.habitat[,"IUCN Code A1"] <- paste(d.habitat[,"IUCN Cat A1"], "A1c")
    d.habitat[,paste("IUCN Code", a)] <- paste0(d.habitat[,paste("IUCN Cat", a)], " ", a, "c")
    d.habitat <- data.frame(d.habitat[,1:2], "N" = nrow(dd()), d.habitat[,3:ncol(d.habitat)], check.names = FALSE)
    
    return(d.habitat)
  })
  
  
  ## Forecast habitat change (in polygon and altitude range) based on time series
  fc <- reactive({
    req(d.habitat())
    
    if (verbose) cat(paste0(Sys.time(), " forecasting habitat change...\n"))
    
    d.habitat <- d.habitat()
    d.habitat[,c("HabitatRemaining_HabitatFraction")] <- paste0(round(d.habitat$Remaining_perc, round.percent), "% (", round(d.habitat$HabitatFraction_perc, round.percent), "%)")
    d.habitat[,c("HabitatReduction_NonHabitatFraction")] <- paste0(round(d.habitat$Reduction_perc, round.percent), "% (", round(100-d.habitat$HabitatFraction_perc, round.percent), " %)")
    
    # interpolate using loess smoother
    d <- loess(formula = paste(paste0("Habitat_", unit), "~", "Year"), data = d.habitat)
    
    # inter <- seq(min(d.habitat[,"Year"], na.rm = TRUE), max(d.habitat[,"Year"], na.rm = TRUE), by = 1)
    inter <- seq(input$period1, input$period2, by = 1)
    pred <- predict(object = d, newdata = inter)
    
    # extrapolate using forecast library
    t <- ts(data = pred, start = min(inter), end = max(inter), frequency = 1)
    
    f <- forecast(t, lambda = "auto", h = input$cast)
    fc <- data.frame(Year = start(f$x)[1]:end(f$mean)[1],
                     habitat = c(f$x,f$mean),
                     upper80 = c(rep(NA, length(f$x)), f$upper[,1]),
                     upper95 = c(rep(NA, length(f$x)), f$upper[,2]),
                     lower80 = c(rep(NA, length(f$x)), f$lower[,1]),
                     lower95 = c(rep(NA, length(f$x)), f$lower[,2]),
                     Extent = "EX")
    fc$Remaining_perc <- 100*fc$habitat/fc[fc$Year == input$period1, "habitat"]
    fc$Reduction_perc <- 100-fc$Remaining_perc
    remfr <- (d.habitat$Remaining_perc/d.habitat$HabitatFraction_perc)[1]
    fc$HabitatFraction_perc <- fc$Remaining_perc/remfr
    fc[,c("HabitatRemaining_HabitatFraction")] <- paste0(round(fc$Remaining_perc, round.percent), "% (", round(fc$HabitatFraction_perc, round.percent), "%)")
    fc[,c("HabitatReduction_NonHabitatFraction")] <- paste0(round(fc$Reduction_perc, round.percent), "% (", round(100-fc$HabitatFraction_perc, round.percent), " %)")
    
    return(fc)
  })
  
  ## Habitat change and forecast plot
  plotInput3 <- eventReactive(input$submit, {
    req(d.habitat()) ; req(fc()) ; req(S())
    
    d.habitat <- d.habitat()
    d.habitat[,c("HabitatRemaining_HabitatFraction")] <- paste0(round(d.habitat$Remaining_perc, round.percent), "% (", round(d.habitat$HabitatFraction_perc, round.percent), "%)")
    d.habitat[,c("HabitatReduction_NonHabitatFraction")] <- paste0(round(d.habitat$Reduction_perc, round.percent), "% (", round(100-d.habitat$HabitatFraction_perc, round.percent), " %)")
    fc <- fc()
  
    tit <- paste0("Change in habitat [", input$period1, " - ", input$period2, " at extent ",
                  paste(round(c(xmin(S()),xmax(S()),ymin(S()),ymax(S())), 7), collapse = ", "),
                  ", at ", altmin(), " - ", altmax(), "m a.s.l.]")

    # create plot (with forecasting)
    ph <- ggplot(fc, aes_string(x = "Year", y = "habitat", group = "Extent", label = "HabitatReduction_NonHabitatFraction")) +
      geom_ribbon(aes(ymin = lower95, ymax = upper95), fill = "gray70", alpha = 0.5) +
      geom_ribbon(aes(ymin = lower80, ymax = upper80), fill = "blue", alpha = 0.5) +
      geom_line(colour = "blue") +
      geom_point(data = d.habitat, aes_string(x = "Year", y = paste0("Habitat_", unit)), colour = "#D2691E", size = 2) +
      # geom_text(data = d.habitat, aes_string(x = "Year", y = paste0("Habitat_", unit), label = "HabitatRemaining_HabitatFraction"), colour = "gray15", size = 2,
      #           nudge_y = diff(range(d.habitat[,paste0("Habitat_", unit)]))/50) +
      geom_text(data = d.habitat, aes_string(x = "Year", y = paste0("Habitat_", unit), label = "HabitatReduction_NonHabitatFraction"), colour = "gray15", size = 2,
                nudge_y = -diff(range(d.habitat[,paste0("Habitat_", unit)]))/50) +
      geom_point(data = subset(fc, Year %in% seq(round(as.numeric(input$period2), -1), max(fc$Year, na.rm = TRUE), by = 10)), 
                 aes_string(x = "Year", y = "habitat"), colour = "tomato", shape = 4, size = 2) +
      # geom_text(data = subset(fc, Year %in% seq(round(as.numeric(input$period2), -1), max(fc$Year, na.rm = TRUE), by = 10)), 
      #           aes_string(x = "Year", y = "habitat", label = "HabitatRemaining_HabitatFraction"), colour = "gray15", size = 2,
      #           nudge_y = +diff(range(d.habitat[,paste0("Habitat_", unit)]))/50) +
      geom_text(data = subset(fc, Year %in% seq(round(as.numeric(input$period2), -1), max(fc$Year, na.rm = TRUE), by = 10)),
                aes_string(x = "Year", y = "habitat", label = "HabitatReduction_NonHabitatFraction"), colour = "tomato", size = 2,
                nudge_y = -diff(range(d.habitat[,paste0("Habitat_", unit)]))/50) +
      labs(x = "Year", y = paste0("Habitat (", unit, ")")) +
      theme_minimal() +
      ggtitle(label = "", subtitle = tit)

    if (verbose) cat(paste0(Sys.time(), " done.\n"))
    return(ph)
  })
  output$forecast <- renderPlotly({
    
    suppressWarnings(
      ggplotly(plotInput3(), tooltip = c("x","y","label"), 
               width = input$psize, height = input$psize)
    )
  })
  
  
  ## Criterion A) Habitat Change Table
  TableInput1 <- eventReactive(input$submit, {
    req(dd())
    
    if (verbose) cat(paste0(Sys.time(), " computing habitat reduction table (IUCN criterion A).\n"))
    d.habitat <- d.habitat()
    rownames(d.habitat) <- NULL
    
    # round to round.percent
    d.habitat[,"HabitatFraction_perc"] <- round(d.habitat[,"HabitatFraction_perc"], round.percent)
    d.habitat[,"Remaining_perc"] <- round(d.habitat[,"Remaining_perc"], round.percent)
    d.habitat[,"Reduction_perc"] <- round(d.habitat[,"Reduction_perc"], round.percent)
  
    names(d.habitat) <- gsub(unit, paste0("(", unit, ")"), gsub("pixels", "(pixels)", gsub("perc", "(%)", gsub("_", " ", names(d.habitat)))))
    for (i in names(d.habitat)[grep(paste0(unit, "$"), names(d.habitat))]) d.habitat[,i] <- round(d.habitat[,i], round.percent)
    
    if (verbose) cat(paste0(Sys.time(), " done.\n"))
    return(d.habitat[,grep("(pixels)", names(d.habitat), invert = TRUE)])
  })
  output$table1 <- DT::renderDataTable({TableInput1()})
  
  
  ## Criterion B) AOO / EOO Change Table
  TableInput2 <- eventReactive(input$submit, {
    req(R1()) ; req(R2())
    
    if (verbose) cat(paste0(Sys.time(), " computing EOO / AOO reduction table (IUCN criterion B).\n"))
    d.eoo <- rbind(R1()$iucn.res, R2()$iucn.res)
    rownames(d.eoo) <- NULL
    colnames(d.eoo) <- gsub("#", "Nb", names(d.eoo))
    names(d.eoo)[names(d.eoo) == "Nb subpopulations"] <- "Nb sub-populations"
    names(d.eoo)[names(d.eoo) == "Occurrences in PA"] <- "Occurrences in PA (%)"
    d.eoo[,1] <- as.numeric(gsub(" km<sup>2</sup>", "", d.eoo[,1]))
    d.eoo[,2] <- as.numeric(gsub(" km<sup>2</sup>", "", d.eoo[,2]))
    d.eoo[,9] <- suppressWarnings(as.numeric(gsub("%$", "", d.eoo[,9])))
    d.eoo <- data.frame("Period" = c("Historic", "Extant"),
                        "Extent" = c("Historic EOO", "Extant EOO"),
                        "N" = c(nrow(dd()), nrow(subset(dd(), occurrenceRemarks %in% str.extant))),
                        "EOO (km2)" = d.eoo[,1], "EOO Red (%)" = round(100-(100*d.eoo[,1]/d.eoo[1,1]), round.percent),
                        "AOO (km2)" = d.eoo[,2], "AOO Red (%)" = round(100-(100*d.eoo[,2]/d.eoo[1,2]), round.percent),
                        d.eoo[c(3:6,9,7:8,10:11)], check.names = FALSE)
  
    if (verbose) cat(paste0(Sys.time(), " done.\n"))
    return(d.eoo)
  })
  output$table2 <- DT::renderDataTable({TableInput2()})
  
  
  ## Deforestation Madagascar (Vieilledent et al. 2018 benchmark dataset)
  output$ForestCover <- renderPlotly({
    p1 <- ggplot(d.vieilledent, aes(x = Year, y = ForestCover_kha, group = Source, colour = Source, label = Percent)) +
      geom_point(na.rm = T) +
      geom_line(na.rm = T) +
      labs(x = "Year", y = "Forest Cover (kha)") +
      facet_wrap(~Forest.type) +
      theme_minimal()
    
    ggplotly(p1, tooltip = c("x","y","label"),
             height = input$psize, width = input$psize)
  })
  output$Deforestation <- renderPlotly({
    p2 <- ggplot(d.vieilledent2, aes(x = Year, y = deforestation_kha.yr, group = Source, colour = Source, label = Percent1)) +
      geom_point(na.rm = T) +
      geom_line(na.rm = T) +
      labs(x = "Period", y = "Deforestation (kha per year)") +
      facet_wrap(~Forest.type) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45))
    
    ggplotly(p2, tooltip = c("x","y","label"),
             height = input$psize, width = input$psize)
  })
  output$DeforestationRate <- renderPlotly({
    p3 <- ggplot(d.vieilledent2, aes(x = Year, y = deforestation_rate_..yr, group = Source, colour = Source, label = Percent2)) +
      geom_point(na.rm = T) +
      geom_line(na.rm = T) +
      labs(x = "Period", y = "Deforestation Rate (% per year)") +
      facet_wrap(~Forest.type) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45))
    
    ggplotly(p3, tooltip = c("x","y","label"),
             height = input$psize, width = input$psize)
  })
  
  
  ## Download
  output$down1 <- downloadHandler(
    filename = function(){paste0(ifelse(print.date, gsub("-","", gsub(":", "", gsub(" ", "_", as.character(Sys.time())))), ""), 
                                 ifelse(nchar(fn()) > 0 & print.date, "_", ""), fn(), "_AB.html")},
    content = function(file) {
      htmlwidgets::saveWidget(plotInput1(), file = file)
    }
  )
  output$down2 <- downloadHandler(
    filename = function(){paste0(ifelse(print.date, gsub("-","", gsub(":", "", gsub(" ", "_", as.character(Sys.time())))), ""), 
                                 ifelse(nchar(fn()) > 0 & print.date, "_", ""), fn(), "_A_map.pdf")},
    content = function(file) {
      pdf(file, width = input$size, height = input$size)
      suppressWarnings(print(plotInput2()))
      dev.off()
    }
  )
  output$down3 <- downloadHandler(
    filename = function(){paste0(ifelse(print.date, gsub("-","", gsub(":", "", gsub(" ", "_", as.character(Sys.time())))), ""), 
                                 ifelse(nchar(fn()) > 0 & print.date, "_", ""), fn(), "_A_forecast.pdf")},
    content = function(file) {
      pdf(file, width = input$size, height = input$size)
      suppressWarnings(print(plotInput3()))
      dev.off()
    }
  )
  output$down4 <- downloadHandler(
    filename = function(){paste0(ifelse(print.date, gsub("-","", gsub(":", "", gsub(" ", "_", as.character(Sys.time())))), ""), 
                                 ifelse(nchar(fn()) > 0 & print.date, "_", ""), fn(), "_A.csv")},
    content = function(file) {
      dtab1 <- d.habitat()
      dtab1$method.range <- method.range()
      dtab1$exclude.area <- exclude.area()
      dtab1$alpha <- alpha()
      dtab1$buff.alpha <- buff.alpha()
      
      write.csv(dtab1, file = file, row.names = FALSE, quote = FALSE)
    }
  )
  output$down5 <- downloadHandler(
    filename = function(){paste0(ifelse(print.date, gsub("-","", gsub(":", "", gsub(" ", "_", as.character(Sys.time())))), ""), 
                                 ifelse(nchar(fn()) > 0 & print.date, "_", ""), fn(), "_B.html")},
    content = function(file) {
      htmlwidgets::saveWidget(plotInput5(), file = file)
    }
  )
  output$down6 <- downloadHandler(
    filename = function(){paste0(ifelse(print.date, gsub("-","", gsub(":", "", gsub(" ", "_", as.character(Sys.time())))), ""), 
                                 ifelse(nchar(fn()) > 0 & print.date, "_", ""), fn(), "_B.csv")},
    content = function(file) {
      dtab2 <- TableInput2()
      dtab2$method.range <- method.range()
      dtab2$exclude.area <- exclude.area()
      dtab2$alpha <- alpha()
      dtab2$buff.alpha <- buff.alpha()
      dtab2$Resol_sub_pop <- Resol_sub_pop()
      dtab2$method_protected_area <- method_protected_area()
      dtab2$method_locations <- method_locations
      dtab2$Cell_size_locations = ifelse(method_locations == "fixed_grid", Cell_size_locations(), 10)
      dtab2$Rel_cell_size = ifelse(method_locations == "sliding scale", Rel_cell_size(), 0.05)
      
      write.csv(dtab2, file = file, row.names = FALSE, quote = FALSE)
    }
  )
  
  
  ## Export
  output$export <- downloadHandler(
    filename = function(){paste0(ifelse(print.date, gsub("-","", gsub(":", "", gsub(" ", "_", as.character(Sys.time())))), ""), 
                                 ifelse(nchar(fn()) > 0 & print.date, "_", ""), fn(), ".rda")},
    content = function(file) {
      dd <- dd()
      S <- S()
      s <- s()
      alt <- alt()
      R1 <- R1()
      R2 <- R2()
      r <- r()
      fc <- suppressWarnings(fc())
      tab_a <- d.habitat()
      tab_b <- TableInput2()
      leafletAB <- plotInput1()
      leafletB <- plotInput5()
      params <- list()
      save(dd, S, s, alt, R1, R2, r, fc, tab_a, tab_b, leafletAB, leafletB, d.vieilledent, d.vieilledent2, file = file)
    }
  )
}

###########
### APP ###
###########
shinyApp(ui, server)
