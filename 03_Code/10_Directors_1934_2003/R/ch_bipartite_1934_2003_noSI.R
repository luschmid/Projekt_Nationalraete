library(xtable)
library(openxlsx)
library(haven)
library(dplyr)
library(igraph)
library(readstata13)
library(sf)
library(ggraph)
library(ggplot2)
library(tidygraph)  
library(magrittr)
library(tmap)
library(tidyr)
library(reshape2)
library(patchwork)
library(mclust)
library(readr)

# load map data from Fretz Parchet Robert-Nicoud
load("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/20_Highway_Firm/data/raw/maps/maps.RData")

#Path
setwd("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/14_Network_OneMode")

################################################################################
################################################################################
### Creation of bipartite network, firm-to-firm, director-to-director, mun-to-mun : all years

base_path <- "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/"

#Store all edgelists for all years
#years <- c(1997, 1998, 1999, 2000, 2001, 2002, 2003)
years <- c(1943, 1960, 1962, 1963, 1964, 1965, 1966, 1969, 1972, 1975, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003)

for (year in years) {
  # 1. Load edgelist
  path <- paste0(base_path, "edgeliste_", year, "_bipartite_noSI.dta")
  edgelist <- read_stata(path)
  assign(paste0("edgeliste_", year, "_bi"), edgelist)
  
  # 2. Create bipartite graph
  g_bi <- graph_from_data_frame(edgelist, directed = FALSE)
  V(g_bi)$type <- bipartite_mapping(g_bi)$type
  assign(paste0("g_bi_", year), g_bi)
  
  # 3. Create firm–firm and director–director projections
  projections <- bipartite_projection(g_bi, multiplicity = TRUE)
  
  g_firms <- projections$proj1
  assign(paste0("g_firms_", year), g_firms)
  
  g_directors <- projections$proj2
  assign(paste0("g_directors_", year), g_directors)
}

### Add nodes' attributes

### Firms newtork
base_attr_path <- "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/Node_attribute_comp/"

for (year in years) {
  # 1. Read the firm attribute file for the year
  attr_path <- paste0(base_attr_path, "comp_", year, "_attribute_complete_noSI.dta")
  firm_attrs <- read_stata(attr_path)
  
  # 2. Get the graph for the year (must already exist in env as g_firms_XXXX)
  graph_name <- paste0("g_firms_", year)
  g_firms <- get(graph_name)
  
  # 3. Match nodes using firm ID (V(g)$name vs firm_attrs$N_UCID)
  match_idx <- match(V(g_firms)$name, firm_attrs$N_UCID)
  
  # 4. Assign firm-level attributes 
  # Municipality
  V(g_firms)$gdenr       <- firm_attrs$gdenr[match_idx]
  V(g_firms)$gdenr_2012  <- firm_attrs$gdenr_2012[match_idx]
  V(g_firms)$gdenr_2018  <- firm_attrs$gdenr_2018[match_idx]
  #Canton
  V(g_firms)$ctn       <- firm_attrs$ctn[match_idx]
  #Name
  V(g_firms)$cname       <- firm_attrs$cname_orig[match_idx]
  #Firms info
  V(g_firms)$capital       <- firm_attrs$capital[match_idx]
  V(g_firms)$nb_owners     <- firm_attrs$nb_owners[match_idx]
  #Mun langage
  V(g_firms)$language     <- firm_attrs$language[match_idx]
  V(g_firms)$language_def     <- firm_attrs$language_temp[match_idx]
  
  # 5. Optionally reassign updated graph back into environment
  assign(graph_name, g_firms)
}

### Directors network
base_attr_path_dir <- "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/Node_attribute_pers/"

for (year in years) {
  # 1. Read the director attribute file
  attr_path <- paste0(base_attr_path_dir, "pers_", year, "_attribute_complete_noSI.dta")
  dir_attrs <- read_stata(attr_path)
  
  # 2. Get the director–director graph for the year
  graph_name <- paste0("g_directors_", year)
  g_directors <- get(graph_name)
  
  # 3. Match node names to director IDs (N_id_sug)
  match_idx <- match(V(g_directors)$name, dir_attrs$N_id_sug)
  
  # 4. Assign attributes (you can add more below)
  #gdenr & gdenr 2012 to add!
  V(g_directors)$gdenr_2018 <- dir_attrs$gdenr_2018[match_idx]
  #Characteristics
  V(g_directors)$firstname <- dir_attrs$firstname[match_idx]
  V(g_directors)$lastname <- dir_attrs$lastname[match_idx]
  V(g_directors)$male <- dir_attrs$male[match_idx]
  V(g_directors)$nCH <- dir_attrs$nCH[match_idx]
  V(g_directors)$language_w <- dir_attrs$language_w[match_idx]
  V(g_directors)$ctn <- dir_attrs$ctn[match_idx]
  
  # 5. Save the updated graph back to environment
  assign(graph_name, g_directors)
}

### Bipartite
for (year in years) {
  graph_name <- paste0("g_bi_", year)
  g_bi <- get(graph_name)
  
  # -----------------------------------
  # Director attributes (type == TRUE)
  # -----------------------------------
  dir_attrs <- read_stata(paste0(base_attr_path_dir, "pers_", year, "_attribute_complete_noSI.dta"))
  is_director <- V(g_bi)$type == TRUE
  match_idx_dir <- match(V(g_bi)$name[is_director], dir_attrs$N_id_sug)
  
  #gdenr & gdenr 2012 to add!
  V(g_bi)$gdenr_2018[is_director]   <- dir_attrs$gdenr_2018[match_idx_dir]
  #Characteristics
  V(g_directors)$firstname <- dir_attrs$firstname[match_idx]
  V(g_directors)$lastname <- dir_attrs$lastname[match_idx]
  V(g_bi)$male[is_director] <- dir_attrs$male[match_idx_dir]
  V(g_bi)$nCH[is_director] <- dir_attrs$nCH[match_idx_dir]
  V(g_bi)$language_w[is_director] <- dir_attrs$language_w[match_idx_dir]
  V(g_bi)$ctn[is_director] <- dir_attrs$ctn[match_idx_dir]
  
  # -----------------------------------
  # Firm attributes (type == FALSE)
  # -----------------------------------
  firm_attrs <- read_stata(paste0(base_attr_path, "comp_", year, "_attribute_complete_noSI.dta"))
  is_firm <- V(g_bi)$type == FALSE
  match_idx_firm <- match(V(g_bi)$name[is_firm], firm_attrs$N_UCID)
  
  # Add relevant attributes
  V(g_bi)$gdenr[is_firm]        <- firm_attrs$gdenr[match_idx_firm]
  V(g_bi)$gdenr_2012[is_firm]   <- firm_attrs$gdenr_2012[match_idx_firm]
  V(g_bi)$gdenr_2018[is_firm]   <- firm_attrs$gdenr_2018[match_idx_firm]
  #Canton
  V(g_bi)$ctn[is_firm]       <- firm_attrs$ctn[match_idx_firm]
  #Name
  V(g_bi)$cname[is_firm]       <- firm_attrs$cname_orig[match_idx_firm]
  #Firms info
  V(g_bi)$capital[is_firm]       <- firm_attrs$capital[match_idx_firm]
  V(g_bi)$nb_owners[is_firm]     <- firm_attrs$nb_owners[match_idx_firm]
  
  # Reassign to environment
  assign(graph_name, g_bi)
}

## Add agglo dummy in graphs
mun$dummy_agglo <- ifelse(mun$ID0 %in% mun_agglo$ID0, 1, 0)
agglo_lookup <- setNames(mun$dummy_agglo, mun$ID0)

for (year in years) {
  g <- get(paste0("g_firms_", year))
  
  # Add agglomeration dummy
  V(g)$dummy_agglo <- agglo_lookup[V(g)$gdenr]
  
  assign(paste0("g_firms_", year), g)
}


################################################################################
### Measuring distances to cantonal borders
# Marking internal cantonal borders        
mun_2012 <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/data/municipalities_1960_2012_1.dta")
# Merge mun_2012 with mun shape file 
mun_2012 <- mun_2012 %>%
  mutate(gdenr_2012 = as.character(gdenr_2012))

border <- st_union(mun)

canton <- mun %>% left_join(mun_2012 %>% select(gdenr_2012, ctn, ctnnr, dummy_german, dummy_french, dummy_ital), by = c("ID0" = "gdenr_2012"))

canton <- canton %>%
  mutate(
    ctn   = ifelse(is.na(ctn), "TI", ctn),
    ctnnr = ifelse(is.na(ctnnr), 21, ctnnr)
  )

sum(is.na(canton$ctn))      # should be 0
length(unique(canton$ctnnr)) # should be 26 if all cantons are there

# Dissolve municipalities into canton geometries
cantons <- canton %>%
  group_by(ctnnr) %>%
  summarise() %>%   # this automatically dissolves by canton
  ungroup()

# Extract internal canton borders
# Outer Swiss border
swiss_border <- st_union(cantons)

# Get canton boundaries (linework, not polygons)
all_canton_borders <- st_boundary(cantons)

# Internal borders = canton boundaries minus outer border
internal_borders <- st_difference(st_union(all_canton_borders), st_boundary(swiss_border))

map <- tm_shape(border)+
  tm_borders(lwd = 1)+ 
  tm_shape(mun)+
  tm_borders(lwd = 0.3, col = "black")+
  tm_shape(lake)+
  tm_fill("lightblue")+
  tm_borders("darkblue", lwd = 0.2)+
  tm_shape(mun_agglo)+
  tm_fill("grey80")+
  tm_borders(lwd = 0.3, col = "grey")+
  tm_shape(border)+ 
  tm_borders(lwd = 1, col = "black")+ 
  tm_shape(cities)+
  tm_dots(col = "black", size = 0.3)+
  tm_layout(frame = FALSE)+
  tm_shape(internal_borders)+ 
  tm_lines(lwd = 1.5, col = "red")   

print(map)

# Centroids and distance to borders (simple subsample, one-dimensional distance used)
mun <- st_transform(mun, 2056)
cantons <- st_transform(cantons, 2056)
internal_borders <- st_transform(internal_borders, 2056)

#centroids
mun_centroids <- st_centroid(mun)

# Compute distance (each row: one municipality centroid)
distances <- st_distance(mun_centroids, internal_borders)

# Take the minimum distance per municipality
mun$dist_to_canton_border <- apply(distances, 1, min)
mun$dist_to_canton_border_km <- mun$dist_to_canton_border / 1000

summary(mun$dist_to_canton_border)   # distribution of distances

median_dist <- median(mun$dist_to_canton_border, na.rm = TRUE)
print(median_dist)

# Dummy: within 10km of a cantonal border
mun$dummy_10km <- ifelse(mun$dist_to_canton_border <= 10000, 1, 0)
mun$dummy_10km <- factor(mun$dummy_10km, levels = c(0,1), labels = c("0","1"))

# Dummy: within 5km of a cantonal border
mun$dummy_5km  <- ifelse(mun$dist_to_canton_border <= 5000, 1, 0)
mun$dummy_5km <- factor(mun$dummy_5km, levels = c(0,1), labels = c("0","1"))

mun$dummy_median <- ifelse(mun$dist_to_canton_border <= median_dist, 1, 0)
mun$dummy_median <- factor(mun$dummy_median, levels = c(0,1), labels = c("0","1"))

#Subsample of those within distances (5, median, 10) & not in agglo

mun$dummy_5km_noagglo <- ifelse(mun$dummy_5km == 1 & mun$dummy_agglo == 0, 1, 0)
mun$dummy_5km_noagglo <- factor(mun$dummy_5km_noagglo, levels = c(0,1), labels = c("0","1"))

mun$dummy_10km_noagglo <- ifelse(mun$dummy_10km == 1 & mun$dummy_agglo == 0, 1, 0)
mun$dummy_10km_noagglo <- factor(mun$dummy_10km_noagglo, levels = c(0,1), labels = c("0","1"))

mun$dummy_median_noagglo <- ifelse(mun$dummy_median == 1 & mun$dummy_agglo == 0, 1, 0)
mun$dummy_median_noagglo <- factor(mun$dummy_median_noagglo, levels = c(0,1), labels = c("0","1"))

### Mun touching borders only
# st_intersects / st_touches to detect geometries that touch the border line
touches_border <- st_intersects(mun, internal_borders, sparse = FALSE)

# Logical to dummy
mun$dummy_touch_border <- ifelse(rowSums(touches_border) > 0, 1, 0)
mun$dummy_touch_border <- factor(mun$dummy_touch_border, levels = c(0,1), labels = c("0","1"))

mun$dummy_touch_border_noagglo <- ifelse(mun$dummy_touch_border == 1 & mun$dummy_agglo == 0, 1, 0)
mun$dummy_touch_border_noagglo <- factor(mun$dummy_touch_border_noagglo, levels = c(0,1), labels = c("0","1"))

### Agglomeration mun

mun$dummy_agglo <- ifelse(mun$dummy_agglo == 1, 1, 0)
mun$dummy_agglo <- factor(mun$dummy_agglo, levels = c(0,1), labels = c("0","1"))

mun$dummy_nonborder_agglo <- ifelse(mun$dummy_touch_border == 0 & mun$dummy_agglo == 1, 1, 0)
mun$dummy_nonborder_agglo <- factor(mun$dummy_nonborder_agglo, levels = c(0,1), labels = c("0","1"))

mun$dummy_median_agglo <- ifelse(mun$dist_to_canton_border <= median_dist & mun$dummy_agglo == 1, 1, 0)
mun$dummy_median_agglo <- factor(mun$dummy_median_agglo, levels = c(0,1), labels = c("0","1"))

mun$dummy_border_agglo <- ifelse(mun$dummy_touch_border == 1 & mun$dummy_agglo == 1, 1, 0)
mun$dummy_border_agglo <- factor(mun$dummy_border_agglo, levels = c(0,1), labels = c("0","1"))

### Add all dummies in graph
#dist_to_canton_border
dummy_lookup <- setNames(mun$dist_to_canton_border, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dist_to_canton_border <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dist_to_canton_border_km
dummy_lookup <- setNames(mun$dist_to_canton_border_km, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dist_to_canton_border_km <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dummy_median
dummy_lookup <- setNames(mun$dummy_median, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dummy_median <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dummy_median_noagglo
dummy_lookup <- setNames(mun$dummy_median_noagglo, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dummy_median_noagglo <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dummy_touch_border
dummy_lookup <- setNames(mun$dummy_touch_border, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dummy_touch_border <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dummy_touch_border_noagglo
dummy_lookup <- setNames(mun$dummy_touch_border_noagglo, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dummy_touch_border_noagglo <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dummy_nonborder_agglo
dummy_lookup <- setNames(mun$dummy_nonborder_agglo, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dummy_nonborder_agglo <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dummy_median_agglo
dummy_lookup <- setNames(mun$dummy_median_agglo, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dummy_median_agglo <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dummy_border_agglo
dummy_lookup <- setNames(mun$dummy_border_agglo, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dummy_border_agglo <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

################################################################################

### Measuring distances to language border

language <- mun %>% left_join(mun_2012 %>% select(gdenr_2012, gdenr_2018, ctn, ctnnr, dummy_german, dummy_french, dummy_ital), by = c("ID0" = "gdenr_2012"))

language <- language %>%
  mutate(
    dummy_german   = ifelse(is.na(dummy_german), 0, dummy_german),
    dummy_french   = ifelse(is.na(dummy_french), 0, dummy_french),
    dummy_ital   = ifelse(is.na(dummy_ital), 0, dummy_ital)
  )

# 3. Create a single language region variable
language$language_region <- case_when(
  language$dummy_german == 1 ~ "German",
  language$dummy_french == 1 ~ "French",
  language$dummy_ital   == 1 ~ "Italian",
  TRUE ~ "Other"
)

# 4. Dissolve municipalities into language regions
language_regions <- language %>%
  group_by(language_region) %>%
  summarise(.groups = "drop")  # st_union is automatic

# 5. Extract internal borders
# Union of all regions = Switzerland outline
swiss_border <- st_union(language_regions)

# All language boundaries
all_language_borders <- st_union(st_boundary(language_regions))

# Internal borders = all borders minus the outer Swiss border
internal_language_borders <- st_sym_difference(
  all_language_borders,
  st_boundary(swiss_border)
)

# 6. Plot
map <- tm_shape(border)+
  tm_borders(lwd = 1)+ 
  tm_shape(mun)+
  tm_borders(lwd = 0.3, col = "black")+
  tm_shape(lake)+
  tm_fill("lightblue")+
  tm_borders("darkblue", lwd = 0.2)+
  tm_shape(mun_agglo)+
  tm_fill("grey80")+
  tm_borders(lwd = 0.3, col = "grey")+
  tm_shape(border)+ 
  tm_borders(lwd = 1, col = "black")+ 
  tm_shape(cities)+
  tm_dots(col = "black", size = 0.3)+
  tm_layout(frame = FALSE)+
  tm_shape(internal_language_borders)+ 
  tm_lines(lwd = 1.5, col = "red")   

print(map)

### Focusing on Röstigraben
bfv <- language %>% 
  filter(ctn %in% c("BE", "FR", "VS"))

bfv_regions <- bfv %>%
  group_by(language_region) %>%
  summarise(.groups = "drop")

bfv_canton_lang <- bfv %>%
  group_by(ctn, language_region) %>%
  summarise(.groups = "drop")

# Get German and French polygons
german_poly <- bfv_regions %>% filter(language_region == "German")
french_poly <- bfv_regions %>% filter(language_region == "French")

# Intersection of their boundaries = shared border
rostigraben_raw <- st_intersection(st_boundary(german_poly), st_boundary(french_poly))

# Outline of BE/FR/VS
bfv_outline <- st_union(bfv)

# Clip border to this outline
rostigraben_clipped <- st_intersection(rostigraben_raw, bfv_outline)

# Explode into LINESTRING pieces
rostigraben_parts <- st_cast(rostigraben_clipped, "LINESTRING")

map <- tm_shape(bfv) +
  tm_fill("language_region", palette = c("darksalmon", "orange")) +
  tm_borders(lwd = 0.3, col = "grey40") +
  tm_shape(rostigraben_parts) +
  tm_lines(lwd = 2, col = "red") +
  tm_shape(border)+
  tm_borders(lwd = 1)+ 
  tm_shape(mun)+
  tm_borders(lwd = 0.3, col = "black")+
  tm_shape(lake)+
  tm_fill("lightblue")+
  tm_borders("darkblue", lwd = 0.2)+ 
  tm_shape(internal_borders)+ 
  tm_lines(lwd = 2, col = "black")+
  tm_layout(frame = FALSE)

print(map)
tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/Rostigraben.pdf", width = 11.69, height = 11.69)


#For each canton, extract German-French borders inside it
rostigraben_list <- list()

for (canton in c("BE", "FR", "VS")) {
  canton_regions <- bfv_canton_lang %>% filter(ctn == canton)
  
  german_poly <- canton_regions %>% filter(language_region == "German")
  french_poly <- canton_regions %>% filter(language_region == "French")
  
  if (nrow(german_poly) > 0 && nrow(french_poly) > 0) {
    border <- st_intersection(st_boundary(german_poly), st_boundary(french_poly))
    rostigraben_list[[canton]] <- border
  }
}

rostigraben_within_canton <- do.call(rbind, rostigraben_list)

# Centroids and distance to borders (simple subsample, one-dimensional distance used)
rostigraben_parts <- st_transform(rostigraben_parts, 2056)
rostigraben_within_canton <- st_transform(rostigraben_within_canton, 2056)

#centroids
mun_centroids <- st_centroid(mun)

# Compute distance (each row: one municipality centroid)
distances <- st_distance(mun_centroids, rostigraben_within_canton)

# Take the minimum distance per municipality
mun$dist_to_rosti_border <- apply(distances, 1, min)
mun$dist_to_rosti_border_km <- mun$dist_to_rosti_border / 1000

summary(mun$dist_to_rosti_border_km)   # distribution of distances

# Adding ctn & gdenr_2018 to mun 
mun <- mun %>%
  left_join(mun_2012 %>% select(gdenr_2012, gdenr_2018, ctn, ctnnr),
            by = c("ID0" = "gdenr_2012"))

# Dummy: within 15km of a rostigraben
mun$dummy_15km_rosti <- ifelse(mun$ctn %in% c("BE", "FR", "VS") & mun$dist_to_rosti_border <= 15000, 1, 0)
mun$dummy_15km_rosti <- factor(mun$dummy_15km_rosti, levels = c(0,1), labels = c("0","1"))

# Dummy: within 10km of a rostigraben
mun$dummy_10km_rosti <- ifelse(mun$ctn %in% c("BE", "FR", "VS") & mun$dist_to_rosti_border <= 10000, 1, 0)
mun$dummy_10km_rosti <- factor(mun$dummy_10km_rosti, levels = c(0,1), labels = c("0","1"))

# Dummy: within 5km of a rostigraben
mun$dummy_5km_rosti <- ifelse(mun$ctn %in% c("BE", "FR", "VS") & mun$dist_to_rosti_border <= 5000, 1, 0)
mun$dummy_5km_rosti <- factor(mun$dummy_5km_rosti, levels = c(0,1), labels = c("0","1"))

# Dummy: within 2.5km of a rostigraben
mun$dummy_25km_rosti <- ifelse(mun$ctn %in% c("BE", "FR", "VS") & mun$dist_to_rosti_border <= 2500, 1, 0)
mun$dummy_25km_rosti <- factor(mun$dummy_5km_rosti, levels = c(0,1), labels = c("0","1"))

### Mun touching borders only
# st_intersects / st_touches to detect geometries that touch the border line
touches_rosti_border <- st_intersects(mun, rostigraben_within_canton, sparse = FALSE)

# Logical to dummy
mun$dummy_touch_rosti_border <- ifelse(mun$ctn %in% c("BE", "FR", "VS") & rowSums(touches_rosti_border) > 0, 1, 0)
mun$dummy_touch_rosti_border <- factor(mun$dummy_touch_rosti_border, levels = c(0,1), labels = c("0","1"))

### Add all dummies in graph
#dist_to_rosti_border
dummy_lookup <- setNames(mun$dist_to_rosti_border, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dist_to_rosti_border <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dist_to_rosti_border_km
dummy_lookup <- setNames(mun$dist_to_rosti_border_km, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dist_to_rosti_border_km <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dummy_15km_rosti
dummy_lookup <- setNames(mun$dummy_15km_rosti, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dummy_15km_rosti <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dummy_10km_rosti
dummy_lookup <- setNames(mun$dummy_10km_rosti, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dummy_10km_rosti <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dummy_5km_rosti
dummy_lookup <- setNames(mun$dummy_5km_rosti, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dummy_5km_rosti <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dummy_2.5km_rosti
dummy_lookup <- setNames(mun$dummy_25km_rosti, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dummy_25km_rosti <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

#dummy_touch_rosti_border
dummy_lookup <- setNames(mun$dummy_touch_rosti_border, mun$ID0)

for (year in years) {
  g <- get(paste0("g_bi_", year))
  
  # Add agglomeration dummy
  V(g)$dummy_touch_rosti_border <- dummy_lookup[V(g)$gdenr]
  
  assign(paste0("g_bi_", year), g)
}

################################################################################
### Largest Connected Component (LCC)
for (year in years) {
  # Get the bipartite network for the year
  g_firm <- get(paste0("g_bi_", year))
  
  # Identify the largest component
  comp <- components(g_firm)
  giant_id <- which.max(comp$csize)
  lcc_nodes <- which(comp$membership == giant_id)
  
  # Induce subgraph (largest connected component)
  g_lcc <- induced_subgraph(g_firm, vids = lcc_nodes)
  
  # Assign to environment with name like g_bi_lcc_1943
  assign(paste0("g_bi_lcc_", year), g_lcc)
}

### Municipality Coordinates
muni_coords_lv95 <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/Gde_GdeNr2018_Centroid_1931-2018_tmp.dta" )

# Create an sf object using LV95 CRS (EPSG:2056)
muni_sf <- st_as_sf(
  muni_coords_lv95,
  coords = c("E_CNTR", "N_CNTR"),
  crs = 2056  # LV95
)

# Transform to WGS84 (EPSG:4326)
muni_sf_wgs84 <- st_transform(muni_sf, crs = 4326)

# Extract coordinates back into columns
coords_wgs84 <- st_coordinates(muni_sf_wgs84)

# Add lon/lat columns
muni_coords_lv95$lon_wgs84 <- coords_wgs84[, 1]
muni_coords_lv95$lat_wgs84 <- coords_wgs84[, 2]

################################################################################
################################################################################
### Full Bimodale network : general information before going into LCC details

bi_network_summaries <- list()

# Loop over years
for (year in years) {
  
  # Get the firm–firm graph object
  g_bi <- get(paste0("g_bi_", year))
  
  # Separate vertices by type
  type_vec <- V(g_bi)$type
  v_firms <- V(g_bi)[type_vec == FALSE]  
  v_directors <- V(g_bi)[type_vec == TRUE]
  
  # Degrees
  deg <- degree(g_bi)
  
  # Compute macro-level stats (unweighted)
  net_sum <- list(
    year = year,
    
    # General
    number_of_edges = gsize(g_bi),
    number_of_nodes_total = gorder(g_bi),
    
    # Firm-type nodes
    number_of_firm_nodes = length(v_firms),
    avg_degree_firm = mean(deg[v_firms]),
    max_degree_firm = max(deg[v_firms]),
    
    # Other-type nodes
    number_of_dir_nodes = length(v_directors),
    avg_degree_dir = mean(deg[v_directors]),
    max_degree_dir = max(deg[v_directors]),
    
    # Connectivity
    density = edge_density(g_bi),
    is_connected = is_connected(g_bi),
    number_of_components = components(g_bi)$no,
    share_isolates = sum(deg == 0) / gorder(g_bi),
    size_LCC = max(components(g_bi)$csize),
    share_LCC = max(components(g_bi)$csize) / gorder(g_bi)
  )
  
  # Store results
  bi_network_summaries[[as.character(year)]] <- net_sum
}

# Convert to data.frame for summary or plotting
bi_network_summary_df <- bind_rows(bi_network_summaries)

# View
print(bi_network_summary_df) 

# Cleaning of table
bi_network_summary_df_clean <- bi_network_summary_df %>%
  mutate(
    year = as.integer(year),
    density = round(density, 3),
    number_of_components = as.integer(number_of_components),
    number_of_edges = as.integer(number_of_edges),
    number_of_nodes_total = as.integer(number_of_nodes_total),
    number_of_firm_nodes = as.integer(number_of_firm_nodes),
    number_of_dir_nodes = as.integer(number_of_dir_nodes),
    max_degree_dir = as.integer(max_degree_dir),
    max_degree_firm = as.integer(max_degree_firm),
    size_LCC = as.integer(size_LCC)
  )

# Convert to LaTeX
latex_table <- xtable(bi_network_summary_df_clean,
                      caption = "Macro-level measures of the bipartite network",
                      label = "tab:bipartite_net_summary")

# Print to .tex file
print(latex_table,
      include.rownames = FALSE,  
      file = "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/bipartite/bipartite_net_summary.tex",
      booktabs = TRUE)   

######
### bipartite network analysis (in LCC)

# Create empty list to store results
bi_network_centralization <- list()

# Loop over years
for (year in years) {
  # Get the firm–firm graph object
  g_bi_lcc <- get(paste0("g_bi_lcc_", year))
  
  # Separate vertices by type
  type_vec <- V(g_bi)$type
  v_firms <- V(g_bi)[type_vec == FALSE]  
  v_directors <- V(g_bi)[type_vec == TRUE]
  
  # Degrees
  deg <- degree(g_bi_lcc)
  
  # Compute macro-level stats (unweighted)
  net_centr <- list(
    year = year,
    
    # Global network-level metrics (do not split)
    number_of_edges = gsize(g_bi_lcc),
    number_of_nodes = gorder(g_bi_lcc),
    density = edge_density(g_bi_lcc),
    transitivity = transitivity(g_bi_lcc, type = "global"),
    assortativity_degree = assortativity_degree(g_bi_lcc),
    degree_centralization = centr_degree(g_bi_lcc)$centralization,
    clo_centralization = centr_clo(g_bi_lcc)$centralization,
    betw_centralization = centr_betw(g_bi_lcc)$centralization,
    eigen_centralization = centr_eigen(g_bi_lcc)$centralization,
    
    # Separate by node type
    number_of_firm_nodes = sum(V(g_bi_lcc)$type == FALSE),
    number_of_dir_nodes = sum(V(g_bi_lcc)$type == TRUE),
    
    # Degree
    avg_degree_firm = mean(degree(g_bi_lcc)[V(g_bi_lcc)$type == FALSE]),
    max_degree_firm = max(degree(g_bi_lcc)[V(g_bi_lcc)$type == FALSE]),
    avg_degree_dir = mean(degree(g_bi_lcc)[V(g_bi_lcc)$type == TRUE]),
    max_degree_dir = max(degree(g_bi_lcc)[V(g_bi_lcc)$type == TRUE]),
    
    # Strength (weighted degree)
    avg_strength_firm = mean(strength(g_bi_lcc, weights = E(g_bi_lcc)$weight)[V(g_bi_lcc)$type == FALSE]),
    max_strength_firm = max(strength(g_bi_lcc, weights = E(g_bi_lcc)$weight)[V(g_bi_lcc)$type == FALSE]),
    avg_strength_dir = mean(strength(g_bi_lcc, weights = E(g_bi_lcc)$weight)[V(g_bi_lcc)$type == TRUE]),
    max_strength_dir = max(strength(g_bi_lcc, weights = E(g_bi_lcc)$weight)[V(g_bi_lcc)$type == TRUE]),
    
    # Optional: betweenness and closeness can also be split by type
    avg_betw_firm = mean(betweenness(g_bi_lcc, normalized = TRUE)[V(g_bi_lcc)$type == FALSE]),
    avg_betw_dir = mean(betweenness(g_bi_lcc, normalized = TRUE)[V(g_bi_lcc)$type == TRUE]),
    
    avg_clo_firm = mean(closeness(g_bi_lcc, normalized = TRUE)[V(g_bi_lcc)$type == FALSE]),
    avg_clo_dir = mean(closeness(g_bi_lcc, normalized = TRUE)[V(g_bi_lcc)$type == TRUE])
  )
  # Store
  bi_network_centralization[[as.character(year)]] <- net_centr
}

# Convert to data.frame for summary or plotting
bi_network_centralization_df <- bind_rows(bi_network_centralization)

# View
print(bi_network_centralization_df) 

# Cleaning of table
bi_network_centralization_df_clean <- bi_network_centralization_df %>%
  mutate(
    year = as.integer(year),
    
    # Counts
    number_of_edges = as.integer(number_of_edges),
    number_of_nodes = as.integer(number_of_nodes),
    number_of_firm_nodes = as.integer(number_of_firm_nodes),
    number_of_dir_nodes = as.integer(number_of_dir_nodes),
    
    # Global metrics
    density = round(as.numeric(density), 5),
    transitivity = round(as.numeric(transitivity), 5),
    assortativity_degree = round(as.numeric(assortativity_degree), 5),
    degree_centralization = round(as.numeric(degree_centralization), 5),
    clo_centralization = round(as.numeric(clo_centralization), 5),
    betw_centralization = round(as.numeric(betw_centralization), 5),
    eigen_centralization = round(as.numeric(eigen_centralization), 5),
    
    # Degree stats
    avg_degree_firm = round(as.numeric(avg_degree_firm), 3),
    max_degree_firm = as.integer(max_degree_firm),
    avg_degree_dir = round(as.numeric(avg_degree_dir), 3),
    max_degree_dir = as.integer(max_degree_dir),
    
    # Strength stats
    avg_strength_firm = round(as.numeric(avg_strength_firm), 3),
    max_strength_firm = round(as.numeric(max_strength_firm), 0),
    avg_strength_dir = round(as.numeric(avg_strength_dir), 3),
    max_strength_dir = round(as.numeric(max_strength_dir), 0),
    
    # Centrality averages
    avg_betw_firm = round(as.numeric(avg_betw_firm), 5),
    avg_betw_dir = round(as.numeric(avg_betw_dir), 5),
    avg_clo_firm = round(as.numeric(avg_clo_firm), 5),
    avg_clo_dir = round(as.numeric(avg_clo_dir), 5)
  )

# Convert to LaTeX
latex_table <- xtable(bi_network_centralization_df_clean,
                      caption = "Macro-level measures of the Bipartite network",
                      label = "tab:bipartite_net_summary")

# Print to .tex file
print(latex_table,
      include.rownames = FALSE,  
      file = "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/bipartite/bi_LCC_summary.tex",
      booktabs = TRUE)

################################################################################
### Node-level centrality measure
# Storage
all_full_ranks <- list()
all_top25      <- list()
all_avg_rank   <- list()

# Output directory
outdir <- "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/bipartite/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Loop over years
for (year in years) {
  g_bi_lcc <- get(paste0("g_bi_lcc_", year))
  
  # Centralities
  V(g_bi_lcc)$degree      <- degree(g_bi_lcc)
  V(g_bi_lcc)$betweenness <- betweenness(g_bi_lcc, normalized = TRUE)
  V(g_bi_lcc)$closeness   <- closeness(g_bi_lcc, normalized = TRUE)
  V(g_bi_lcc)$eigenvector <- eigen_centrality(g_bi_lcc)$vector
  V(g_bi_lcc)$pagerank    <- page_rank(g_bi_lcc)$vector
  
  # Node type (firms vs directors)
  node_type <- ifelse(V(g_bi_lcc)$type == FALSE, "firm", "dir")
  
  centrality_df <- data.frame(
    node_id    = V(g_bi_lcc)$name,
    node_label = V(g_bi_lcc)$cname,   # works for firms; for directors may be NA if cname missing
    node_type  = node_type,
    degree      = V(g_bi_lcc)$degree,
    betweenness = V(g_bi_lcc)$betweenness,
    closeness   = V(g_bi_lcc)$closeness,
    eigenvector = V(g_bi_lcc)$eigenvector,
    pagerank    = V(g_bi_lcc)$pagerank,
    stringsAsFactors = FALSE
  )
  
  # Function: full ranking for one measure
  make_full_rank <- function(var) {
    centrality_df %>%
      arrange(desc(.data[[var]])) %>%
      mutate(rank = row_number(),
             centrality = var,
             year = year)
  }
  
  # Stack all measures
  full_rank_df <- bind_rows(
    make_full_rank("degree"),
    make_full_rank("betweenness"),
    make_full_rank("closeness"),
    make_full_rank("eigenvector"),
    make_full_rank("pagerank")
  )
  
  # Top 25 per measure, *by node type*
  top25_df <- full_rank_df %>%
    group_by(year, centrality, node_type) %>%
    slice_head(n = 25) %>%
    ungroup()
  
  # Average rank across measures, by node type
  avg_rank_df <- full_rank_df %>%
    group_by(year, node_id, node_label, node_type) %>%
    summarise(
      avg_rank = mean(rank, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(year, node_type, avg_rank)
  
  top25_avg_rank <- avg_rank_df %>%
    group_by(year, node_type) %>%
    slice_head(n = 25) %>%
    ungroup()
  
  # Store
  all_full_ranks[[as.character(year)]] <- full_rank_df
  all_top25[[as.character(year)]]      <- top25_df
  all_avg_rank[[as.character(year)]]   <- avg_rank_df
  
  # Export LaTeX for top25 per measure & type
  for (ntype in c("firm", "dir")) {
    top25_df_ntype <- top25_df %>% filter(node_type == ntype)
    
    print(xtable(top25_df_ntype,
                 caption = paste("Top 25", ntype, "nodes by centrality measures,", year),
                 label = paste0("tab:top25_", ntype, "_", year)),
          include.rownames = FALSE,
          file = file.path(outdir, paste0("top25_", ntype, "_", year, ".tex")),
          booktabs = TRUE)
    
    top25_avg_rank_ntype <- top25_avg_rank %>% filter(node_type == ntype)
    
    print(xtable(top25_avg_rank_ntype,
                 caption = paste("Top 25", ntype, "nodes by average centrality rank,", year),
                 label = paste0("tab:top25_avg_", ntype, "_", year)),
          include.rownames = FALSE,
          file = file.path(outdir, paste0("top25_avg_", ntype, "_", year, ".tex")),
          booktabs = TRUE)
  }
}

# Combine across years
all_full_ranks_df <- bind_rows(all_full_ranks)
all_top25_df      <- bind_rows(all_top25)
all_avg_rank_df   <- bind_rows(all_avg_rank)


################################################################################
################################################################################            

### Direct edge counts  (share within/across cantons)

# Empty list to hold results
cross_summary <- list()

for (year in years) {
  g <- get(paste0("g_bi_lcc_", year))
  
  if (ecount(g) == 0 || vcount(g) == 0) {
    if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
    cross_summary[[as.character(year)]] <- data.frame(
      year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
    )
    assign(paste0("g_bi_lcc_", year), g)
    next
  }
  
  # Vertex canton labels
  ctn <- V(g)$ctn
  
  # Endpoints as numeric vertex IDs : allow comparison of vertex attributes
  ep <- ends(g, E(g), names = FALSE)
  
  # Mark bad/missing labels
  bad  <- is.na(ctn[ep[,1]]) | is.na(ctn[ep[,2]]) | ctn[ep[,1]] == "" | ctn[ep[,2]] == ""
  diff <- ctn[ep[,1]] != ctn[ep[,2]]
  
  # Edge attribute: 1 = cross-canton, 0 = within; NA if a label missing
  E(g)$cross_canton <- ifelse(bad, NA_integer_, as.integer(diff))
  
  # Save graph back with new edge attribute
  assign(paste0("g_bi_lcc_", year), g)
  
  # Count edges within/cross cantons
  within_cnt <- sum(E(g)$cross_canton == 0L, na.rm = TRUE)
  cross_cnt  <- sum(E(g)$cross_canton == 1L, na.rm = TRUE)
  tot_defined <- within_cnt + cross_cnt
  cross_summary[[as.character(year)]] <- data.frame(
    year = year,
    within = within_cnt,
    cross  = cross_cnt,
    share_within = if (tot_defined > 0) within_cnt / tot_defined else NA_real_,
    share_cross  = if (tot_defined > 0) cross_cnt  / tot_defined else NA_real_
  )
}

# Bind the summary into a data.frame for plotting
cross_summary_df <- dplyr::bind_rows(cross_summary)
print(cross_summary_df)

# Plot that 
share_long <- cross_summary_df %>%
  select(year, share_within, share_cross) %>%
  pivot_longer(-year, names_to = "measure", values_to = "value")

p_both <- ggplot(share_long, aes(x = year, y = value, linetype = measure)) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE) +
  scale_linetype_manual(values = c("solid","dashed"),
                        labels = c("Share cross-canton", "Share within canton")) +
  labs(x = "Year", y = "Share of edges",
       title = "Within- vs cross-canton edges over time",
       linetype = NULL) +
  theme_minimal(base_size = 14)

print(p_both)
ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/bipartite/fig_within_vs_cross.pdf", p_both, width = 11.69, height = 11.69)

### Direct edge counts  (share within/across languages)

# Storage for summary per year
cross_summary_lang <- list()

for (year in years) {
  g <- get(paste0("g_bi_lcc_", year))   # your LCC graph for that year
  
  if (ecount(g) == 0 || vcount(g) == 0) {
    # no edges or no nodes: set empty attribute + summary
    if (ecount(g) == 0) E(g)$cross_lang <- integer(0)
    cross_summary_lang[[as.character(year)]] <- data.frame(
      year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
    )
    assign(paste0("g_bi_lcc_", year), g)
    next
  }
  
  # Vertex canton labels (character is fine)
  lang <- V(g)$language
  
  # Endpoints as numeric vertex IDs (important!)
  ep <- ends(g, E(g), names = FALSE)
  
  # Mark bad/missing labels
  bad  <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]] == "" | lang[ep[,2]] == ""
  diff <- lang[ep[,1]] != lang[ep[,2]]
  
  # Edge attribute: 1 = cross-canton, 0 = within; NA if a label missing
  E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
  
  # Save graph back with new edge attribute
  assign(paste0("g_firms_lcc_", year), g)
  
  # Yearly summary
  within_lang <- sum(E(g)$cross_lang == 0L, na.rm = TRUE)
  cross_lang  <- sum(E(g)$cross_lang == 1L, na.rm = TRUE)
  tot_defined <- within_lang + cross_lang
  cross_summary_lang[[as.character(year)]] <- data.frame(
    year = year,
    within = within_lang,
    cross  = cross_lang,
    share_within_lang = if (tot_defined > 0) within_lang / tot_defined else NA_real_,
    share_cross_lang  = if (tot_defined > 0) cross_lang  / tot_defined else NA_real_
  )
}

# Bind the optional summary into a data.frame (handy for plotting)
cross_summary_lang_df <- dplyr::bind_rows(cross_summary_lang)
print(cross_summary_lang_df)

# Plot that 
share_long_lang <- cross_summary_lang_df %>%
  select(year, share_within_lang, share_cross_lang) %>%
  pivot_longer(-year, names_to = "measure", values_to = "value")

p_both_lang <- ggplot(share_long_lang, aes(x = year, y = value, linetype = measure)) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE) +
  scale_linetype_manual(values = c("solid","dashed"),
                        labels = c("Share across", "Share within language region")) +
  labs(x = "Year", y = "Share of edges",
       title = "Within- vs across language regions edges over time",
       linetype = NULL) +
  theme_minimal(base_size = 14)

print(p_both_lang)
ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/bipartite/fig_within_vs_cross_lang.pdf", p_both_lang, width = 11.69, height = 11.69)

### Comparing within/across cantons & within/across language regions 
# Unweighted
unw_combined <- bind_rows(
  share_long %>%
    mutate(family = "canton"),
  share_long_lang %>%
    mutate(family = "language")
) %>%
  mutate(
    kind = ifelse(grepl("cross", measure), "cross", "within"),
    color_key = ifelse(family == "canton", "canton", "language"),
    lty_key   = ifelse(kind == "cross", "solid", "dashed")
  )

p_unw <- ggplot(unw_combined, aes(x = year, y = value,
                                  group = interaction(family, kind))) +
  geom_line(aes(color = color_key, linetype = lty_key), na.rm = TRUE) +
  geom_point(aes(color = color_key), na.rm = TRUE, size = 1.8) +
  scale_color_manual(values = c(canton = "red", language = "black"),
                     labels = c("Canton", "Language region"), name = NULL) +
  scale_linetype_identity(guide = "none") +
  labs(x = "Year", y = "Share of edges",
       title = "Within vs. cross edges — Canton (red) vs Language (black)",
       subtitle = "Solid = cross; Dashed = within") +
  theme_minimal(base_size = 13)

print(p_unw)
ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_within_cross_canton_language_unweighted.pdf", p_unw, width = 11.69, height = 11.69)
