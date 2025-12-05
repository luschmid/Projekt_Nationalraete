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

  # Dummy: within 15km of a cantonal border
  mun$dummy_15km  <- ifelse(mun$dist_to_canton_border <= 15000, 1, 0)
  mun$dummy_15km <- factor(mun$dummy_15km, levels = c(0,1), labels = c("0","1"))  
  
# Dummy: within 10km of a cantonal border
  mun$dummy_10km <- ifelse(mun$dist_to_canton_border <= 10000, 1, 0)
  mun$dummy_10km <- factor(mun$dummy_10km, levels = c(0,1), labels = c("0","1"))

# Dummy: within 5km of a cantonal border
  mun$dummy_5km  <- ifelse(mun$dist_to_canton_border <= 5000, 1, 0)
  mun$dummy_5km <- factor(mun$dummy_5km, levels = c(0,1), labels = c("0","1"))
  
# Dummy: Median to cantonal border  
  mun$dummy_median <- ifelse(mun$dist_to_canton_border <= median_dist, 1, 0)
  mun$dummy_median <- factor(mun$dummy_median, levels = c(0,1), labels = c("0","1"))

#Subsample of those within distances (5, median, 10) & not in agglo
  mun$dummy_5km_noagglo <- ifelse(mun$dummy_5km == 1 & mun$dummy_agglo == 0, 1, 0)
  mun$dummy_5km_noagglo <- factor(mun$dummy_5km_noagglo, levels = c(0,1), labels = c("0","1"))
  
  mun$dummy_10km_noagglo <- ifelse(mun$dummy_10km == 1 & mun$dummy_agglo == 0, 1, 0)
  mun$dummy_10km_noagglo <- factor(mun$dummy_10km_noagglo, levels = c(0,1), labels = c("0","1"))
  
  mun$dummy_median_noagglo <- ifelse(mun$dummy_median == 1 & mun$dummy_agglo == 0, 1, 0)
  mun$dummy_median_noagglo <- factor(mun$dummy_median_noagglo, levels = c(0,1), labels = c("0","1"))

#Subsample of bins (0-5, 5-10, 10-15)
# 0–5 km
  mun$dummy_0_5km <- ifelse(mun$dist_to_canton_border <= 5000, 1, 0)
  mun$dummy_0_5km <- factor(mun$dummy_0_5km, levels = c(0,1), labels = c("0","1"))

# 5–10 km
  mun$dummy_5_10km <- ifelse(mun$dist_to_canton_border > 5000 & mun$dist_to_canton_border <= 10000, 1, 0)
  mun$dummy_5_10km <- factor(mun$dummy_5_10km, levels = c(0,1), labels = c("0","1"))

# 10–15 km
  mun$dummy_10_15km <- ifelse(mun$dist_to_canton_border > 10000 & mun$dist_to_canton_border <= 15000, 1, 0)
  mun$dummy_10_15km <- factor(mun$dummy_10_15km, levels = c(0,1), labels = c("0","1"))

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
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dist_to_canton_border <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dist_to_canton_border_km
  dummy_lookup <- setNames(mun$dist_to_canton_border_km, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dist_to_canton_border_km <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dummy_median
  dummy_lookup <- setNames(mun$dummy_median, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_median <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dummy_median_noagglo
  dummy_lookup <- setNames(mun$dummy_median_noagglo, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_median_noagglo <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dummy_touch_border
  dummy_lookup <- setNames(mun$dummy_touch_border, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_touch_border <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dummy_touch_border_noagglo
  dummy_lookup <- setNames(mun$dummy_touch_border_noagglo, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_touch_border_noagglo <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dummy_nonborder_agglo
  dummy_lookup <- setNames(mun$dummy_nonborder_agglo, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_nonborder_agglo <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dummy_median_agglo
  dummy_lookup <- setNames(mun$dummy_median_agglo, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_median_agglo <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dummy_border_agglo
  dummy_lookup <- setNames(mun$dummy_border_agglo, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_border_agglo <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

# 0–5 km
  dummy_lookup <- setNames(mun$dummy_0_5km, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_0_5km <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

# 5–10 km
  dummy_lookup <- setNames(mun$dummy_5_10km, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_5_10km <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

# 10-15 km
  dummy_lookup <- setNames(mun$dummy_10_15km, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_10_15km <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }
  
# 0-15 km
  dummy_lookup <- setNames(mun$dummy_15km, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_15km <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

################################################################################

### Measuring distances to language border

  language <- mun %>% left_join(mun_2012 %>% select(gdenr_2012, ctn, ctnnr, dummy_german, dummy_french, dummy_ital), by = c("ID0" = "gdenr_2012"))
  
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

# Adding ctn to mun 
  mun <- mun %>%
    left_join(mun_2012 %>% select(gdenr_2012, ctn, ctnnr),
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

#Subsample of bins (0-5, 5-10, 10-15)
# 0–5 km
  mun$dummy_0_5km_rosti <- ifelse(mun$dist_to_rosti_border <= 5000, 1, 0)
  mun$dummy_0_5km_rosti <- factor(mun$dummy_0_5km_rosti, levels = c(0,1), labels = c("0","1"))

# 5–10 km
  mun$dummy_5_10km_rosti <- ifelse(mun$dist_to_rosti_border > 5000 & mun$dist_to_rosti_border <= 10000, 1, 0)
  mun$dummy_5_10km_rosti <- factor(mun$dummy_5_10km_rosti, levels = c(0,1), labels = c("0","1"))

# 10–15 km
  mun$dummy_10_15km_rosti <- ifelse(mun$dist_to_rosti_border > 10000 & mun$dist_to_rosti_border <= 15000, 1, 0)
  mun$dummy_10_15km_rosti <- factor(mun$dummy_10_15km_rosti, levels = c(0,1), labels = c("0","1"))

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
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dist_to_rosti_border <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dist_to_rosti_border_km
  dummy_lookup <- setNames(mun$dist_to_rosti_border_km, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dist_to_rosti_border_km <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dummy_15km_rosti
  dummy_lookup <- setNames(mun$dummy_15km_rosti, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_15km_rosti <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dummy_10km_rosti
  dummy_lookup <- setNames(mun$dummy_10km_rosti, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_10km_rosti <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dummy_5km_rosti
  dummy_lookup <- setNames(mun$dummy_5km_rosti, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_5km_rosti <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dummy_25km_rosti
  dummy_lookup <- setNames(mun$dummy_25km_rosti, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_25km_rosti <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

#dummy_touch_rosti_border
  dummy_lookup <- setNames(mun$dummy_touch_rosti_border, mun$ID0)
  
  for (year in years) {
    g <- get(paste0("g_firms_", year))
    
    # Add agglomeration dummy
    V(g)$dummy_touch_rosti_border <- dummy_lookup[V(g)$gdenr]
    
    assign(paste0("g_firms_", year), g)
  }

################################################################################
### Largest Connected Component (LCC)
  for (year in years) {
    # Get the firm network for the year
    g_firm <- get(paste0("g_firms_", year))
    
    # Identify the largest component
    comp <- components(g_firm)
    giant_id <- which.max(comp$csize)
    lcc_nodes <- which(comp$membership == giant_id)
    
    # Induce subgraph (largest connected component)
    g_lcc <- induced_subgraph(g_firm, vids = lcc_nodes)
    
    # Assign to environment with name like g_firms_lcc_1943
    assign(paste0("g_firms_lcc_", year), g_lcc)
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
### Full firm-to-firm network : general information before going into LCC details

  firm_network_summaries <- list()

  # Loop over years
    for (year in years) {
      
      # Get the firm–firm graph object
      g_firm <- get(paste0("g_firms_", year))
      
      # Precompute values
      comp <- components(g_firm)
      size_LCC <- max(comp$csize)
      number_of_nodes <- gorder(g_firm)
      number_of_isolates <- sum(degree(g_firm) == 0)
      
      # Compute macro-level stats (unweighted)
      net_sum <- list(
        year = year,
        number_of_edges = gsize(g_firm),
        number_of_nodes = number_of_nodes,
        size_LCC = size_LCC,
        share_LCC = size_LCC / number_of_nodes,
        size_nonLCC = number_of_nodes - size_LCC,
        number_of_isolates = number_of_isolates,
        nodes_small_comps = number_of_nodes - size_LCC - number_of_isolates,
        share_isolates = number_of_isolates / number_of_nodes,
        density = edge_density(g_firm),
        average_degree = mean(degree(g_firm)),
        max_degree = max(degree(g_firm)),
        is_connected = is_connected(g_firm),
        number_of_components = comp$no
      )
      
      # Store results
      firm_network_summaries[[as.character(year)]] <- net_sum
    }
  
  # Convert to data.frame for summary or plotting
    firm_network_summary_df <- bind_rows(firm_network_summaries)
  
  # View
    print(firm_network_summary_df) 

# Cleaning of table
  firm_network_summary_df_clean <- firm_network_summary_df %>%
    mutate(
      year = as.integer(year),
      number_of_components = as.integer(number_of_components),
      number_of_edges = as.integer(number_of_edges),
      number_of_nodes = as.integer(number_of_nodes),
      max_degree = as.integer(max_degree),
      size_LCC = as.integer(size_LCC)
    )

# Convert to LaTeX
  latex_table <- xtable(firm_network_summary_df_clean,
                        caption = "Macro-level measures of the firm network",
                        label = "tab:firm_net_summary")

# Print to .tex file
  print(latex_table,
        include.rownames = FALSE,  
        file = "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/firm_net_summary.tex",
        booktabs = TRUE)   

  ### Plot: fragmentation of the firm-to-firm network
  p_frag <- ggplot(firm_network_summary_df, aes(x = year)) +

    # Lines
    geom_line(aes(y = number_of_nodes, color = "Total Number of Nodes"), linewidth = 1.2) +
    geom_line(aes(y = size_LCC, color = "Nodes in LCC"), linewidth = 1.2) +
    geom_line(aes(y = nodes_small_comps, color = "Nodes not in LCC or Isolated"), linewidth = 1.2) +
    geom_line(aes(y = number_of_isolates, color = "Isolated Nodes"), linewidth = 1.2) +
    
    # Labels & colors
    scale_color_manual(values = c("Total Number of Nodes" = "black",
                                  "Nodes in LCC" = "blue",
                                  "Nodes not in LCC or Isolated" = "#8856a7",
                                  "Isolated Nodes" = "#7fcdbb")) +
    
    labs(x = "Year",
         y = "Number of Firms",
         title = "Evolution and Fragmentation of the Swiss Firm Network",
         color = NULL) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 18),      # legend labels
      axis.title = element_text(size = 18),       # axis titles
      axis.text  = element_text(size = 18),       # axis tick labels
      plot.title = element_text(size = 18, face = "bold"), # title
    )
  
  p_frag <- p_frag +
    coord_cartesian(xlim = c(1960, 2003))
  
  print(p_frag)
  
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_1990s_60to2000.pdf",
         p_frag, width = 11.69, height = 8.27)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_1990s_60to2000.png",
         p_frag, width = 11.69, height = 8.27, dpi = 300)
  
  p_frag_2 <- p_frag +
    coord_cartesian(xlim = c(1970, 2003))
  
  print(p_frag_2)
  
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_1990s_70to2000.pdf",
         p_frag_2, width = 11.69, height = 8.27)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_1990s_70to2000.png",
         p_frag_2, width = 11.69, height = 8.27, dpi = 300)
  
  p_frag_3 <- p_frag +
    coord_cartesian(xlim = c(1980, 2003))
  
  print(p_frag_3)
  
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_1990s_80to2000.pdf",
         p_frag_3, width = 11.69, height = 8.27)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_1990s_80to2000.png",
         p_frag_3, width = 11.69, height = 8.27, dpi = 300)
  
######
### Firm to firm network analysis (in LCC)

# Create empty list to store results
  firm_network_centralization <- list()

# Loop over years
  for (year in years) {
    # Get the firm–firm graph object
    g_firm_lcc <- get(paste0("g_firms_lcc_", year))
    
    # Compute macro-level stats (unweighted)
    net_centr <- list(
      year = year,
      number_of_edges = gsize(g_firm_lcc),
      number_of_nodes = gorder(g_firm_lcc),
      density = edge_density(g_firm_lcc),
      average_degree = mean(degree(g_firm_lcc)),
      max_degree = max(degree(g_firm_lcc)),
      transitivity = transitivity(g_firm_lcc, type = "global"),
      assortativity_degree = assortativity_degree(g_firm_lcc),
      degree_centralization = centr_degree(g_firm_lcc)$centralization,
      clo_centralization = centr_clo(g_firm_lcc)$centralization,
      betw_centralization = centr_betw(g_firm_lcc)$centralization,
      eigen_centralization = centr_betw(g_firm_lcc)$centralization,
      average_strength = mean(strength(g_firm_lcc, weights = E(g_firm_lcc)$weight)),
      max_strength = max(strength(g_firm_lcc, weights = E(g_firm_lcc)$weight))
    )
    # Store
    firm_network_centralization[[as.character(year)]] <- net_centr
    }

# Convert to data.frame for summary or plotting
  firm_network_centralization_df <- bind_rows(firm_network_centralization)

# View
  print(firm_network_centralization_df) 

# Cleaning of table
  firm_network_centralization_df_clean <- firm_network_centralization_df %>%
    mutate(
      year = as.integer(year),
      density = round(density, 5),
      number_of_edges = as.integer(number_of_edges),
      number_of_nodes = as.integer(number_of_nodes),
      max_degree = as.integer(max_degree)
     )

# Convert to LaTeX
  latex_table <- xtable(firm_network_summary_df_clean,
                        caption = "Macro-level measures of the firm network",
                        label = "tab:firm_net_summary")

# Print to .tex file
  print(latex_table,
        include.rownames = FALSE,  
        file = "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/firm_LCC_summary.tex",
        booktabs = TRUE)

################################################################################
################################################################################            

### Direct edge counts  (share within/across cantons, with 95% CI)

# Empty list to hold results
  cross_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))   # call graph
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
      cross_summary[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, 
        share_within = NA_real_, share_cross = NA_real_,
        ci95_low_within = NA_real_, ci95_high_within = NA_real_,
        ci95_low_cross = NA_real_, ci95_high_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_", year), g)
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
  assign(paste0("g_firms_lcc_", year), g)
  
# Count edges within/cross cantons
  within_cnt <- sum(E(g)$cross_canton == 0L, na.rm = TRUE)
  cross_cnt  <- sum(E(g)$cross_canton == 1L, na.rm = TRUE)
  tot_defined <- within_cnt + cross_cnt
  
# Confidence intervals (only if edges exist)
  if (tot_defined > 0) {
    mean_cross  <- cross_cnt / tot_defined
    mean_within <- within_cnt / tot_defined
    
    se_cross  <- sqrt(mean_cross  * (1 - mean_cross)  / tot_defined)
    se_within <- sqrt(mean_within * (1 - mean_within) / tot_defined)
    
    ci95_low_cross  <- max(0, mean_cross  - 1.96 * se_cross)
    ci95_high_cross <- min(1, mean_cross  + 1.96 * se_cross)
    ci95_low_within <- max(0, mean_within - 1.96 * se_within)
    ci95_high_within<- min(1, mean_within + 1.96 * se_within)
  } else {
    mean_cross <- mean_within <- NA
    ci95_low_cross <- ci95_high_cross <- NA
    ci95_low_within <- ci95_high_within <- NA
  }
  
  cross_summary[[as.character(year)]] <- data.frame(
    year = year,
    cross  = cross_cnt,
    share_cross  = mean_cross,
    ci95_low_cross   = ci95_low_cross,
    ci95_high_cross  = ci95_high_cross
  )
}

# Bind the summary into a data.frame for plotting
  cross_summary_df <- dplyr::bind_rows(cross_summary)
  print(cross_summary_df)

# Reshape for plotting
  share_long <- cross_summary_df %>%
    select(year, share_cross, ci95_low_cross, ci95_high_cross) %>%
    pivot_longer(cols = c(share_cross),
                 names_to = "measure", values_to = "value") %>%
    mutate(
      ci_low  = ifelse(measure == "share_cross", ci95_low_cross),
      ci_high = ifelse(measure == "share_cross", ci95_high_cross)
    )

# Plot with 95% CI ribbon
  p_both <- ggplot(share_long, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = measure),
                alpha = 0.2, colour = NA) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Share cross-canton")) +
    scale_fill_manual(values = c("grey70","grey40"), guide = "none") +
    labs(x = "Year", y = "Share of edges",
         title = "Cross-canton edges over time (95% CI)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_within_vs_cross.pdf", p_both, width = 11.69, height = 11.69)

### Direct edge counts (share within/across languages, with 95% CI)

# Storage for summary per year
  cross_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))   # your LCC graph for that year
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      if (ecount(g) == 0) E(g)$cross_lang <- integer(0)
      cross_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, 
        share_within_lang = NA_real_, share_cross_lang = NA_real_,
        ci95_low_within_lang = NA_real_, ci95_high_within_lang = NA_real_,
        ci95_low_cross_lang  = NA_real_, ci95_high_cross_lang  = NA_real_
      )
      assign(paste0("g_firms_lcc_", year), g)
      next
    }
  
# Vertex language labels
  lang <- V(g)$language
  
# Endpoints as numeric vertex IDs
  ep <- ends(g, E(g), names = FALSE)
  
# Mark bad/missing labels
  bad  <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]] == "" | lang[ep[,2]] == ""
  diff <- lang[ep[,1]] != lang[ep[,2]]
  
# Edge attribute: 1 = cross-language, 0 = within
  E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
  
# Save graph back with new edge attribute
  assign(paste0("g_firms_lcc_", year), g)
  
# Yearly counts
  within_lang <- sum(E(g)$cross_lang == 0L, na.rm = TRUE)
  cross_lang  <- sum(E(g)$cross_lang == 1L, na.rm = TRUE)
  tot_defined <- within_lang + cross_lang
  
# Confidence intervals
  if (tot_defined > 0) {
    mean_cross  <- cross_lang / tot_defined
    mean_within <- within_lang / tot_defined
    
    se_cross  <- sqrt(mean_cross  * (1 - mean_cross)  / tot_defined)
    se_within <- sqrt(mean_within * (1 - mean_within) / tot_defined)
    
    ci95_low_cross  <- max(0, mean_cross  - 1.96 * se_cross)
    ci95_high_cross <- min(1, mean_cross  + 1.96 * se_cross)
    ci95_low_within <- max(0, mean_within - 1.96 * se_within)
    ci95_high_within<- min(1, mean_within + 1.96 * se_within)
  } else {
    mean_cross <- mean_within <- NA
    ci95_low_cross <- ci95_high_cross <- NA
    ci95_low_within <- ci95_high_within <- NA
  }
  
  cross_summary_lang[[as.character(year)]] <- data.frame(
    year = year,
    cross  = cross_lang,
    share_cross_lang  = mean_cross,
    ci95_low_cross_lang   = ci95_low_cross,
    ci95_high_cross_lang  = ci95_high_cross
  )
}

# Bind into a data.frame
  cross_summary_lang_df <- dplyr::bind_rows(cross_summary_lang)
  print(cross_summary_lang_df)

# Reshape for plotting
  share_long_lang <- cross_summary_lang_df %>%
    select(year, share_cross_lang, ci95_low_cross_lang, ci95_high_cross_lang) %>%
    pivot_longer(cols = c(share_cross_lang),
                 names_to = "measure", values_to = "value") %>%
    mutate(
      ci_low  = ifelse(measure == "share_cross_lang", ci95_low_cross_lang),
      ci_high = ifelse(measure == "share_cross_lang", ci95_high_cross_lang)
    )

# Plot with 95% CI ribbon
  p_both_lang <- ggplot(share_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = measure),
                alpha = 0.2, colour = NA) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Share across language regions")) +
    scale_fill_manual(values = c("grey70","grey40"), guide = "none") +
    labs(x = "Year", y = "Share of edges",
         title = "Across language regions edges over time (95% CI)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_lang)

  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_within_vs_cross_lang.pdf",
       p_both_lang, width = 11.69, height = 11.69)

  
### Add weighted cases (for both langage and cantonal borders)  
  
################################################################################
################################################################################            
### E-I (External-Internal) Index   
  
  ei_results <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    
    cc <- E(g)$cross_canton
    keep <- !is.na(cc)
    
    if (!any(keep)) {
      ei_results[[as.character(year)]] <- data.frame(year = year, ei_index = NA_real_)
      next
    }
    
    E <- sum(cc[keep] == 1L)
    I <- sum(cc[keep] == 0L)
    ei_index <- (E - I) / (E + I)
    
    ei_results[[as.character(year)]] <- data.frame(
      year = year,
      ei_index = ei_index,
      within_edges = I,
      cross_edges = E
    )
  }
  
  ei_df <- dplyr::bind_rows(ei_results)
  
### E-I Index weighted
  
  eiw_results <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    if (ecount(g) == 0) { eiw_results[[as.character(year)]] <- data.frame(year=year, ei_w=NA_real_); next }
    
  # ensure cross_canton exists; otherwise create it as before
    if (is.null(E(g)$cross_canton)) {
      ep  <- ends(g, E(g), names = FALSE); ctn <- V(g)$ctn
      bad <- is.na(ctn[ep[,1]]) | is.na(ctn[ep[,2]]) | ctn[ep[,1]]=="" | ctn[ep[,2]]==""
      diff <- ctn[ep[,1]] != ctn[ep[,2]]
      E(g)$cross_canton <- ifelse(bad, NA_integer_, as.integer(diff))
    }
    
    w <- if (!is.null(E(g)$weight)) E(g)$weight else
      if (!is.null(E(g)$n_shared_directors)) E(g)$n_shared_directors else rep(1, ecount(g))
    
    keep <- !is.na(E(g)$cross_canton) & !is.na(w)
    if (!any(keep)) { eiw_results[[as.character(year)]] <- data.frame(year=year, ei_w=NA_real_); next }
    
    Iw <- sum(w[keep & E(g)$cross_canton == 0L])
    Ew <- sum(w[keep & E(g)$cross_canton == 1L])
    ei_w <- (Ew - Iw) / (Ew + Iw)
    
    eiw_results[[as.character(year)]] <- data.frame(year = year, ei_w = ei_w)
  }
  
  ei_weighted_df <- dplyr::bind_rows(eiw_results)              
  
### Plot that
# Combine the two series
  ei_both <- ei_df %>%
    select(year, ei_unweighted = ei_index) %>%
    full_join(ei_weighted_df %>% select(year, ei_weighted = ei_w), by = "year") %>%
    arrange(year)
  
  print(ei_both)
  
# Long format plot
  ei_long <- ei_both %>%
    pivot_longer(cols = c(ei_unweighted, ei_weighted),
                 names_to = "type", values_to = "value") %>%
    mutate(type = recode(type,
                         ei_unweighted = "Unweighted E–I",
                         ei_weighted   = "Weighted E–I"))
  
# Plot
  p_ei <- ggplot(ei_long, aes(x = year, y = value, linetype = type)) +
    geom_hline(yintercept = 0, color = "gray60", linetype = "dashed") +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
    labs(x = "Year",
         y = "E–I index",
         title = "Evolution of E–I index over time (unweighted vs. weighted)") +
    theme_minimal(base_size = 14)
  
  print(p_ei)
  
# Save PDF for Overleaf
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_ei_unw_vs_w.pdf", p_ei, width = 11.69, height = 11.69)
  
  # Notes
  #EI=−1: all ties are internal (only within cantons).
  #EI=+1: all ties are external (only across cantons).
  #EI=0: balance of internal and external ties (in relative terms).  
  
### E-I (External-Internal) Index for Language   
  ei_results_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    
    cl <- E(g)$cross_lang
    keep <- !is.na(cl)
    
    if (!any(keep)) {
      ei_results_lang[[as.character(year)]] <- data.frame(year = year, ei_index = NA_real_)
      next
    }
    
    E <- sum(cl[keep] == 1L)
    I <- sum(cl[keep] == 0L)
    ei_index <- (E - I) / (E + I)
    
    ei_results_lang[[as.character(year)]] <- data.frame(
      year = year,
      ei_index = ei_index,
      within_edges = I,
      cross_edges = E
    )
  }
  
  ei_lang_df <- dplyr::bind_rows(ei_results_lang)        
  
### E-I Index for Language weighted
  
  eiw_results_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    if (ecount(g) == 0) { eiw_results_lang[[as.character(year)]] <- data.frame(year=year, ei_w=NA_real_); next }
    
    # ensure cross_lang exists; otherwise create it as before
    if (is.null(E(g)$cross_lang)) {
      ep  <- ends(g, E(g), names = FALSE); lang <- V(g)$language
      bad <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]]=="" | lang[ep[,2]]==""
      diff <- lang[ep[,1]] != lang[ep[,2]]
      E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
    }
    
    w <- if (!is.null(E(g)$weight)) E(g)$weight else
      if (!is.null(E(g)$n_shared_directors)) E(g)$n_shared_directors else rep(1, ecount(g))
    
    keep <- !is.na(E(g)$cross_lang) & !is.na(w)
    if (!any(keep)) { eiw_results_lang[[as.character(year)]] <- data.frame(year=year, ei_w=NA_real_); next }
    
    Iw <- sum(w[keep & E(g)$cross_lang == 0L])
    Ew <- sum(w[keep & E(g)$cross_lang == 1L])
    ei_w <- (Ew - Iw) / (Ew + Iw)
    
    eiw_results_lang[[as.character(year)]] <- data.frame(year = year, ei_w = ei_w)
  }
  
  ei_weighted_lang_df <- dplyr::bind_rows(eiw_results_lang)                      
  
### Plot that
  # Combine the two series
    ei_lang_both <- ei_lang_df %>%
      select(year, ei_unweighted = ei_index) %>%
      full_join(ei_weighted_lang_df %>% select(year, ei_weighted = ei_w), by = "year") %>%
      arrange(year)
    
    print(ei_lang_both)
    
  # Long format plot
    ei_lang_long <- ei_lang_both %>%
      pivot_longer(cols = c(ei_unweighted, ei_weighted),
                   names_to = "type", values_to = "value") %>%
      mutate(type = recode(type,
                           ei_unweighted = "Unweighted E–I",
                           ei_weighted   = "Weighted E–I"))
    
  # Plot
    p_ei <- ggplot(ei_lang_long, aes(x = year, y = value, linetype = type)) +
      geom_hline(yintercept = 0, color = "gray60", linetype = "dashed") +
      geom_line(na.rm = TRUE) +
      geom_point(na.rm = TRUE) +
      scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
      labs(x = "Year",
           y = "E–I index",
           title = "Evolution of E–I index over time (unweighted vs. weighted)") +
      theme_minimal(base_size = 14)
    
    print(p_ei)   
    
  # Save PDF for Overleaf
    ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_ei_unw_vs_w_lang.pdf", p_ei, width = 11.69, height = 11.69)


################################################################################
################################################################################            
  
### Matrices (26x26 at canton level, diag is within, rest is across)       
  
  make_canton_matrix_rowshare <- function(g, weighted = FALSE) {
    ep <- ends(g, E(g), names = FALSE)
    c1 <- V(g)$ctn[ep[,1]]
    c2 <- V(g)$ctn[ep[,2]]
    
  # drop missing/empty labels
    keep <- !is.na(c1) & !is.na(c2) & c1 != "" & c2 != ""
    c1 <- c1[keep]; c2 <- c2[keep]
    
    cantons <- sort(unique(V(g)$ctn[V(g)$ctn != "" & !is.na(V(g)$ctn)]))
    c1 <- factor(c1, levels = cantons)
    c2 <- factor(c2, levels = cantons)
    
    if (weighted) {
      w <- if (!is.null(E(g)$weight)) E(g)$weight else
        if (!is.null(E(g)$n_shared_directors)) E(g)$n_shared_directors else
          rep(1, ecount(g))
      w <- w[keep]
      M_dir <- xtabs(w ~ c1 + c2)
    } else {
      M_dir <- xtabs(~ c1 + c2)
    }
    
  # Symmetrize (undirected graph)
    M_sym <- M_dir + t(M_dir)
    diag(M_sym) <- diag(M_dir) + diag(M_dir)
    
  # Row-normalize: each row sums to 1
    row_sums <- rowSums(M_sym)
    M_rowshare <- sweep(M_sym, 1, row_sums, FUN = "/")
    M_rowshare
  }
  
### All years
  # Symmetric canton–canton matrix
  make_canton_matrix <- function(g, weight_attr = NULL) {
    if (ecount(g) == 0 || vcount(g) == 0) return(matrix(0, 0, 0))
    
    ep <- ends(g, E(g), names = FALSE)
    c1 <- V(g)$ctn[ep[,1]]
    c2 <- V(g)$ctn[ep[,2]]
    
    keep <- !is.na(c1) & !is.na(c2) & c1 != "" & c2 != ""
    c1 <- c1[keep]; c2 <- c2[keep]
    
    if (length(c1) == 0) return(matrix(0, 0, 0))
    
    cantons <- sort(unique(V(g)$ctn[!is.na(V(g)$ctn) & V(g)$ctn != ""]))
    a <- factor(pmin(c1, c2), levels = cantons)
    b <- factor(pmax(c1, c2), levels = cantons)
    
  # weights
    if (is.null(weight_attr)) {
      w <- rep(1, length(c1))
    } else {
      w_all <- E(g)[[weight_attr]]
      if (is.null(w_all)) w_all <- rep(1, ecount(g))
      w <- w_all[keep]
    }
    
    U <- xtabs(w ~ a + b)  # upper triangle
    
    M <- U + t(U)
    diag(M) <- diag(U)  # keep diagonal only once
    M
  }
  
  # Row-normalize
    row_normalize <- function(M) {
      if (length(M) == 0) return(M)
      rs <- rowSums(M)
      M_row <- sweep(M, 1, rs, FUN = "/")
      M_row[!is.finite(M_row)] <- 0
      M_row
    }
  
  # Heatmap function
    plot_rowshare <- function(M, year) {
      df <- melt(as.matrix(M), varnames = c("RowCanton","ColCanton"), value.name = "share")
      ggplot(df, aes(x = ColCanton, y = RowCanton, fill = share)) +
        geom_tile(color = "white") +
        geom_text(aes(label = sprintf("%.2f", share)),
                  size = 2, color = "white") +
        scale_fill_viridis_c(
          option = "cividis",
          limits = c(0, 0.8),
          trans = "sqrt"   # or "log"
        ) +
        coord_equal() +
        labs(title = paste("Row-normalized canton–canton links (", year, ")", sep=""),
             subtitle = "Each row sums to 1 (share of a canton’s edges)",
             fill = "Share") +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  
    plot_rowshare_clean <- function(M, year) {
      df <- reshape2::melt(as.matrix(M), varnames = c("RowCanton","ColCanton"), value.name = "share")
      ggplot(df, aes(x = ColCanton, y = RowCanton, fill = share)) +
        geom_tile(color = "white") +
        scale_fill_viridis_c(
          option = "cividis",
          limits = c(0, 0.8),
          trans = "sqrt"   # or "log"
        ) +
        coord_equal() +
        labs(title = paste0(year), fill = "Share") +
        theme_minimal(base_size = 9) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "right"
        )
    }
  
  # Loop over all years
    outdir <- "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/canton_heatmaps"
    dir.create(outdir, showWarnings = FALSE)
    
    for (year in years) {
      g <- get(paste0("g_firms_lcc_", year))
      M  <- make_canton_matrix(g)       # unweighted; add "weight" if you prefer weighted
      M_row <- row_normalize(M)
      
      if (length(M_row) > 0) {
        p <- plot_rowshare(M_row, year)
        ggsave(file.path(outdir, paste0("canton_heatmap_", year, ".pdf")),
               p, width = 8, height = 7)
      }
    }
  
  # Plot several years to illustrate some changes
    sel_years <- c(1943, 1960, 1980, 1990, 2000, 2003)
    
    plots <- list()
    for (year in sel_years) {
      g <- get(paste0("g_firms_lcc_", year))
      M <- make_canton_matrix(g)
      M_row <- row_normalize(M)
      if (length(M_row) > 0) {
        plots[[as.character(year)]] <- plot_rowshare_clean(M_row, year)
      }
    }
    
  # Arrange them in 2x3 grid
    multi_plot <- (plots[["1943"]] | plots[["1960"]] | plots[["1980"]]) /
      (plots[["1990"]] | plots[["2000"]] | plots[["2003"]])
    
  # Shared legend + one big title
    multi_plot <- multi_plot + plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
    
    final_plot <- wrap_plots(multi_plot) +
      plot_annotation(title = "Row-normalized canton–canton connectivity, selected years")
  
    ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/canton_heatmaps/canton_heatmaps_selected_years.pdf", final_plot,
           width = 11.7, height = 11.7)   # A4 landscape 
    
################################################################################
################################################################################            
    
### Assortativity by canton & language
    
### Assortativity with canton
    
    assortativity_results <- list()
    
    for (year in years) {
      g <- get(paste0("g_firms_lcc_", year))   # largest connected component graph
      
      if (vcount(g) == 0 || ecount(g) == 0) {
        assortativity_results[[as.character(year)]] <- data.frame(
          year = year,
          assortativity_canton = NA_real_
        )
        next
      }
      
      # Convert canton to numeric factor for assortativity_nominal
      ctn <- V(g)$ctn
      ctn_num <- as.numeric(as.factor(ctn))
      
      assortativity_results[[as.character(year)]] <- data.frame(
        year = year,
        assortativity_canton = assortativity_nominal(g, ctn_num, directed = FALSE)
      )
    }
    
    assortativity_df <- dplyr::bind_rows(assortativity_results)
    
    print(assortativity_df)
    
### Assortativity, language
    assortativity_lang_results <- list()
    
    for (year in years) {
      g <- get(paste0("g_firms_lcc_", year))
      
      if (vcount(g) == 0 || ecount(g) == 0) {
        assortativity_lang_results[[as.character(year)]] <- data.frame(
          year = year,
          assortativity_lang = NA_real_
        )
        next
      }
      
      lang <- V(g)$language
      lang[is.na(lang) | lang == ""] <- "Unknown"
      lang_num <- as.numeric(as.factor(lang))
      
      assortativity_lang_results[[as.character(year)]] <- data.frame(
        year = year,
        assortativity_lang = assortativity_nominal(g, lang_num, directed = FALSE)
      )
    }
    
    assortativity_lang_df <- dplyr::bind_rows(assortativity_lang_results)
    
    print(assortativity_lang_df)
    
    ### Combined plot
    assortativity_all <- assortativity_df %>%
      left_join(assortativity_lang_df, by = "year") %>%
      tidyr::pivot_longer(cols = c(assortativity_canton, assortativity_lang),
                          names_to = "measure", values_to = "value") %>%
      mutate(measure = dplyr::recode(measure,
                                     assortativity_canton = "Canton",
                                     assortativity_lang   = "Language"))
    
    assortativity_graph <- ggplot(assortativity_all, aes(x = year, y = value, color = measure)) +
      geom_line(linewidth = 1, na.rm = TRUE) +
      geom_point(size = 2, na.rm = TRUE) +
      labs(x = "Year", y = "Assortativity",
           title = "Assortativity of firm-to-firm networks") +
      theme_minimal(base_size = 14)
    
    print(assortativity_graph)
    
    ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/assortativity_graph.pdf", assortativity_graph, width = 11.69, height = 11.69)
    
### Weighted Assortativity (Canton)
    assortativity_results_w <- list()
    
    for (year in years) {
      g <- get(paste0("g_firms_lcc_", year))   # LCC graph for that year
      
      if (vcount(g) == 0 || ecount(g) == 0) {
        assortativity_results_w[[as.character(year)]] <- data.frame(
          year = year,
          assortativity_canton_w = NA_real_
        )
        next
      }
      
      # Extract attributes
      ctn <- V(g)$ctn
      ctn_num <- as.numeric(as.factor(ctn))
      
      # Edge list + weights
      ep <- ends(g, E(g), names = FALSE)
      if (!"weight" %in% edge_attr_names(g)) {
        w <- rep(1, ecount(g))
      } else {
        w <- E(g)$weight
      }
      
      
      # Weighted observed same-attribute fraction
      same <- (ctn_num[ep[,1]] == ctn_num[ep[,2]])
      obs  <- sum(w * same, na.rm = TRUE) / sum(w, na.rm = TRUE)
      
      # Weighted expected fraction (by marginal distribution)
      s <- strength(g, weights = w)
      total_w <- sum(w)
      p_attr <- tapply(s, ctn_num, sum, na.rm = TRUE) / (2 * total_w)
      exp <- sum(p_attr^2)
      
      # Weighted assortativity
      r_w <- obs - exp
      
      assortativity_results_w[[as.character(year)]] <- data.frame(
        year = year,
        assortativity_canton_w = r_w
      )
    }
    
    assortativity_df_w <- dplyr::bind_rows(assortativity_results_w)
    print(assortativity_df_w)

### Weighted Assortativity (Language)
    assortativity_lang_results_w <- list()
    
    for (year in years) {
      g <- get(paste0("g_firms_lcc_", year))
      
      if (vcount(g) == 0 || ecount(g) == 0) {
        assortativity_lang_results_w[[as.character(year)]] <- data.frame(
          year = year,
          assortativity_lang_w = NA_real_
        )
        next
      }
      
      lang <- V(g)$language
      lang[is.na(lang) | lang == ""] <- "Unknown"
      lang_num <- as.numeric(as.factor(lang))
      
      ep <- ends(g, E(g), names = FALSE)
      if (!"weight" %in% edge_attr_names(g)) {
        w <- rep(1, ecount(g))
      } else {
        w <- E(g)$weight
      }
      
      same <- (lang_num[ep[,1]] == lang_num[ep[,2]])
      obs  <- sum(w * same, na.rm = TRUE) / sum(w, na.rm = TRUE)
      
      s <- strength(g, weights = w)
      total_w <- sum(w)
      p_attr <- tapply(s, lang_num, sum, na.rm = TRUE) / (2 * total_w)
      exp <- sum(p_attr^2)
      
      r_w <- obs - exp
      
      assortativity_lang_results_w[[as.character(year)]] <- data.frame(
        year = year,
        assortativity_lang_w = r_w
      )
    }
    
    assortativity_lang_df_w <- dplyr::bind_rows(assortativity_lang_results_w)
    print(assortativity_lang_df_w)
    
### Combined graph
    
    assortativity_all_w <- assortativity_df_w %>%
      left_join(assortativity_lang_df_w, by = "year") %>%
      tidyr::pivot_longer(cols = c(assortativity_canton_w, assortativity_lang_w),
                          names_to = "measure", values_to = "value") %>%
      mutate(measure = dplyr::recode(measure,
                                     assortativity_canton_w = "Canton (weighted)",
                                     assortativity_lang_w   = "Language (weighted)"))
    
    assortativity_graph_w <- ggplot(assortativity_all_w, aes(x = year, y = value, color = measure)) +
      geom_line(linewidth = 1, na.rm = TRUE) +
      geom_point(size = 2, na.rm = TRUE) +
      labs(x = "Year", y = "Weighted assortativity",
           title = "Weighted assortativity of firm-to-firm networks") +
      theme_minimal(base_size = 14)
    
    print(assortativity_graph_w)
    
    ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/assortativity_graph_w.pdf", assortativity_graph_w, width = 11.69, height = 11.69)
    
################################################################################
################################################################################            
    
### Node-level homophily      
    #Assortativity (global): one single correlation coefficient across all edges (a global summary).
    #Node-level homophily: computes local proportions of same-type neighbors for each node, then averages them.
    # Would allow heterogeneity analysis more easily : take nodes in agglo and compute, take nodes in bordering regions, ...
    
# Homophily by canton
    homophily_results <- list()
    
    for (year in years) {
      g <- get(paste0("g_firms_lcc_", year))
      
      if (vcount(g) == 0 || ecount(g) == 0) {
        homophily_results[[as.character(year)]] <- data.frame(
          year = year,
          mean_homophily_canton = NA_real_
        )
        next
      }
      
      ctn <- V(g)$ctn
      ctn[is.na(ctn) | ctn == ""] <- "Unknown"
      
      neighs <- igraph::adjacent_vertices(g, V(g))
      
      node_hom <- sapply(seq_along(V(g)), function(i) {
        if (length(neighs[[i]]) == 0) return(NA_real_)  # isolates
        neigh_ctn <- ctn[neighs[[i]]]
        mean(neigh_ctn == ctn[i], na.rm = TRUE)
      })
      
      homophily_results[[as.character(year)]] <- data.frame(
        year = year,
        mean_homophily_canton = mean(node_hom, na.rm = TRUE)
      )
    }
    
    homophily_canton_df <- dplyr::bind_rows(homophily_results)
    print(homophily_canton_df)
    
# Homophily by language
    homophily_lang_results <- list()
    
    for (year in years) {
      g <- get(paste0("g_firms_lcc_", year))
      
      if (vcount(g) == 0 || ecount(g) == 0) {
        homophily_lang_results[[as.character(year)]] <- data.frame(
          year = year,
          mean_homophily_lang = NA_real_
        )
        next
      }
      
      lang <- V(g)$language
      lang[is.na(lang) | lang == ""] <- "Unknown"
      
      neighs <- igraph::adjacent_vertices(g, V(g))
      
      node_hom <- sapply(seq_along(V(g)), function(i) {
        if (length(neighs[[i]]) == 0) return(NA_real_)
        neigh_lang <- lang[neighs[[i]]]
        mean(neigh_lang == lang[i], na.rm = TRUE)
      })
      
      homophily_lang_results[[as.character(year)]] <- data.frame(
        year = year,
        mean_homophily_lang = mean(node_hom, na.rm = TRUE)
      )
    }
    
    homophily_lang_df <- dplyr::bind_rows(homophily_lang_results)
    print(homophily_lang_df)
    
# Plot that
    homophily_all <- homophily_canton_df %>%
      left_join(homophily_lang_df, by = "year") %>%
      tidyr::pivot_longer(cols = c(mean_homophily_canton, mean_homophily_lang),
                          names_to = "measure", values_to = "value") %>%
      mutate(measure = dplyr::recode(measure,
                                     mean_homophily_canton = "Canton",
                                     mean_homophily_lang   = "Language"))
    
    homophily_graph <- ggplot(homophily_all, aes(x = year, y = value, color = measure)) +
      geom_line(linewidth = 1, na.rm = TRUE) +
      geom_point(size = 2, na.rm = TRUE) +
      scale_color_manual(values = c("Canton" = "red", "Language" = "black")) +
      labs(x = "Year", y = "Mean node-level homophily",
           title = "Node-level homophily in firm-to-firm networks",
           subtitle = "Average share of neighbors in the same group") +
      theme_minimal(base_size = 14) +
      theme(legend.title = element_blank())
    
    print(homophily_graph)
    
    ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/homophily_graph.pdf", homophily_graph, width = 11.69, height = 11.69)

################################################################################
    
### Louvain community detection on LCC of firm-to-firm networks    
    
    # Same w/ LCC
    for (year in years) {
      g_firm <- get(paste0("g_firms_lcc_", year))
      if (gorder(g_firm) > 0 && gsize(g_firm) > 0) {
        comm_firm <- cluster_louvain(g_firm, weights = E(g_firm)$weight)
        V(g_firm)$community <- membership(comm_firm)
        assign(paste0("g_firms_lcc_", year), g_firm)
        assign(paste0("comm_firms_", year), comm_firm)
      }
    }
    
  # Initialize
    firm_community_evolution <- data.frame(
      year = integer(),
      num_firms = integer(),
      num_communities = integer(),
      communities_per_1000_firms = numeric()
    )
    
    for (year in years) {
      comm_name <- paste0("comm_firms_", year)
      g_name <- paste0("g_firms_lcc_", year)
      
      if (exists(comm_name) && exists(g_name)) {
        comm <- get(comm_name)
        g <- get(g_name)
        
        num_comm <- length(unique(membership(comm)))
        num_nodes <- gorder(g)
        comm_per_1000 <- num_comm / num_nodes * 1000
        
        firm_community_evolution <- rbind(
          firm_community_evolution,
          data.frame(
            year = year,
            num_firms = num_nodes,
            num_communities = num_comm,
            communities_per_1000_firms = comm_per_1000
          )
        )
      }
    }
    
  # View results
    print(firm_community_evolution)
    
  # Optional plot
    plot(
      firm_community_evolution$year,
      firm_community_evolution$communities_per_1000_firms,
      type = "b", pch = 16,
      xlab = "Year", ylab = "Communities per 1000 firms",
      main = "Normalized Firm Community Count (Louvain)"
    )
    
### Comparing communities and cantonal borders in LCC
    
    for (year in years) {
      g <- get(paste0("g_firms_lcc_", year))
      
      # Only if both attributes are present
      if (!is.null(V(g)$community) && !is.null(V(g)$ctn)) {
        tab_firm_canton <- table(V(g)$community, V(g)$ctn)
        assign(paste0("tab_firm_canton_", year), tab_firm_canton)
      }
    }
    
  # Adjusted Rand Index
    if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
    library(mclust)  # For adjustedRandIndex()
    
    ari_results <- data.frame(
      year = integer(),
      adjusted_rand_index = numeric()
    )
    
    for (year in years) {
      g_name <- paste0("g_firms_lcc_", year)
      g <- get(g_name)
      
  # Ensure both attributes are present
      if (!is.null(V(g)$community) && !is.null(V(g)$ctn)) {
        # Remove nodes with missing canton
        valid <- !is.na(V(g)$community) & !is.na(V(g)$ctn)
        ari_lcc <- adjustedRandIndex(V(g)$community[valid], V(g)$ctn[valid])
        
        ari_results <- rbind(
          ari_results,
          data.frame(year = year, adjusted_rand_index_lcc = ari_lcc)
        )
      }
    }
    
  # View ARI values over time
    print(ari_results)
    
  # plot
    plot(
      ari_results$year,
      ari_results$adjusted_rand_index_lcc,
      type = "b", pch = 16,
      xlab = "Year", ylab = "Adjusted Rand Index",
      main = "Similarity Between Firm Communities and Cantons (in LCC)"
    )
    
### Comparing communities and langage borders in LCC    
    
    for (year in years) {
      g <- get(paste0("g_firms_lcc_", year))
      
      # Only if both attributes are present
      if (!is.null(V(g)$community) && !is.null(V(g)$language)) {
        tab_firm_language <- table(V(g)$community, V(g)$language)
        assign(paste0("tab_firm_language_", year), tab_firm_language)
      }
    }
    
  #Adjusted Rand Index
    if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
    library(mclust)  # For adjustedRandIndex()
    
    ari_results <- data.frame(
      year = integer(),
      adjusted_rand_index = numeric()
    )
    
    for (year in years) {
      g_name <- paste0("g_firms_lcc_", year)
      g <- get(g_name)
      
  #Ensure both attributes are present
      if (!is.null(V(g)$community) && !is.null(V(g)$language)) {
        # Remove nodes with missing canton
        valid <- !is.na(V(g)$community) & !is.na(V(g)$language)
        ari_lcc <- adjustedRandIndex(V(g)$community[valid], V(g)$language[valid])
        
        ari_results <- rbind(
          ari_results,
          data.frame(year = year, adjusted_rand_index_lcc = ari_lcc)
        )
      }
    }
    
  # View ARI values over time
    print(ari_results)
    
  # plot
    plot(
      ari_results$year,
      ari_results$adjusted_rand_index_lcc,
      type = "b", pch = 16,
      xlab = "Year", ylab = "Adjusted Rand Index",
      main = "Similarity Between Firm Communities and Languages (in LCC)"
    )    

### Label-shuffling as benchmark for ARI
   

################################################################################
### Subsamples 

#   Maps - Cantons
    
  # Maps with shade of colours
    map <- tm_shape(mun) +
      tm_fill("dist_to_canton_border_km", palette = "viridis", style = "quantile") +
      tm_borders(lwd = 0.3, col = "grey") +
      tm_shape(internal_borders) +
      tm_lines(lwd = 1, col = "red") +
      tm_shape(lake)+
      tm_fill("lightblue")+
      tm_borders("darkblue", lwd = 0.2)+
      tm_layout(frame = FALSE)
    
    print(map)
    tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/distance_to_canton_borders.pdf", width = 11.69, height = 11.69)
    
  # Median dist map
    map <- tm_shape(mun) +
      tm_fill("dummy_median", 
              palette = c("white", "darksalmon"), 
              title = "Within Median Distance of Border") +
      tm_borders(lwd = 0.3, col = "grey") +
      tm_shape(internal_borders) +
      tm_lines(lwd = 1, col = "black") +
      tm_shape(lake)+
      tm_fill("lightblue")+
      tm_borders("darkblue", lwd = 0.2)+
      tm_layout(frame = FALSE)      
    print(map)
    tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/median_to_canton_borders.pdf", width = 11.69, height = 11.69)
    
  # Median dist + no agglo map
    map <- tm_shape(mun) +
      tm_fill("dummy_median_noagglo", 
              palette = c("white", "darksalmon"), 
              title = "Within Median Distance of Border, not in agglo") +
      tm_borders(lwd = 0.3, col = "grey") +
      tm_shape(internal_borders) +
      tm_lines(lwd = 1, col = "black") +
      tm_shape(lake)+
      tm_fill("lightblue")+
      tm_borders("darkblue", lwd = 0.2)+
      tm_layout(frame = FALSE) 
    print(map)
    tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/median_to_canton_borders_noagglo.pdf", width = 11.69, height = 11.69)
    
  # Agglo
    map <- tm_shape(mun) +
      tm_fill("dummy_agglo", 
              palette = c("white", "darksalmon"), 
              title = "In agglomeration") +
      tm_borders(lwd = 0.3, col = "grey") +
      tm_shape(internal_borders) +
      tm_lines(lwd = 1, col = "black") +
      tm_shape(lake)+
      tm_fill("lightblue")+
      tm_borders("darkblue", lwd = 0.2)+
      tm_layout(frame = FALSE)
    print(map)
    tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/agglo.pdf", width = 11.69, height = 11.69)
    
  # Bins from cantonal borders
    # Collapse into one categorical variable for distance bins
    mun$border_distance_group <- case_when(
      mun$dist_to_canton_border <= 5000 ~ "0–5 km",
      mun$dist_to_canton_border > 5000 & mun$dist_to_canton_border <= 10000 ~ "5–10 km",
      mun$dist_to_canton_border > 10000 & mun$dist_to_canton_border <= 15000 ~ "10–15 km",
      TRUE ~ "More than 15 km"
    )
    
    # Make it a factor to control legend order
      mun$border_distance_group <- factor(
        mun$border_distance_group,
        levels = c("0–5 km", "5–10 km", "10–15 km", "More than 15 km")
      )
      
    # Map with different colors per group
      map <- tm_shape(mun) +
        tm_fill("border_distance_group",
                palette = c("darksalmon", "tomato", "firebrick", "white"),
                title = "Distance to cantonal border") +
        tm_borders(lwd = 0.3, col = "grey") +
        tm_shape(internal_borders) +
        tm_lines(lwd = 1, col = "black") +
        tm_shape(lake) +
        tm_fill("lightblue") +
        tm_borders("darkblue", lwd = 0.2) +
        tm_layout(frame = FALSE)
      
      print(map)
      tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/bins_to_canton_borders.pdf", width = 11.69, height = 11.69)
      
  #Bins from cantonal borders w/o agglo
    # Distance bins excluding agglomeration municipalities
    mun$border_distance_group_noagglo <- case_when(
      mun$dist_to_canton_border <= 5000 & mun$dummy_agglo == 0~ "0–5 km",
      mun$dist_to_canton_border > 5000 & mun$dist_to_canton_border <= 10000 & mun$dummy_agglo == 0 ~ "5–10 km",
      mun$dist_to_canton_border > 10000 & mun$dist_to_canton_border <= 15000 & mun$dummy_agglo == 0 ~ "10–15 km",
      TRUE ~ "More than 15 km"
    )
    
    # Factor to control legend order
    mun$border_distance_group_noagglo <- factor(
      mun$border_distance_group_noagglo,
      levels = c("0–5 km", "5–10 km", "10–15 km", "More than 15 km")
    )
    
    # Map excluding agglomeration municipalities
    map_noagglo <- tm_shape(mun) +
      tm_fill("border_distance_group_noagglo",
              palette = c("darksalmon", "tomato", "firebrick", "white"),
              title = "Distance to cantonal border (excl. agglomerations)") +
      tm_borders(lwd = 0.3, col = "grey") +
      tm_shape(internal_borders) +
      tm_lines(lwd = 1, col = "black") +
      tm_shape(lake) +
      tm_fill("lightblue") +
      tm_borders("darkblue", lwd = 0.2) +
      tm_layout(frame = FALSE)
    
    print(map_noagglo)   
    tmap_save(map_noagglo, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/bins_to_canton_borders.pdf", width = 11.69, height = 11.69)

# Maps - Langage
    
    #  Maps with shade of colours
    map <- tm_shape(mun) +
      tm_fill("dist_to_rosti_border_km", palette = "viridis", style = "quantile") +
      tm_borders(lwd = 0.3, col = "grey") +
      tm_shape(internal_borders) +
      tm_lines(lwd = 2, col = "black") +
      tm_shape(lake)+
      tm_fill("lightblue")+
      tm_borders("darkblue", lwd = 0.2)+
      tm_shape(rostigraben_within_canton) +
      tm_lines(lwd = 2, col = "red") +
      tm_layout(frame = FALSE)
    print(map)
    tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/distance_to_language_borders.pdf", width = 11.69, height = 11.69)
  
    #  10km from rosti
    map <- tm_shape(mun) +
      tm_fill("dummy_10km_rosti", 
              palette = c("white", "darksalmon"), 
              title = "Touching Röstigraben") +
      tm_borders(lwd = 0.3, col = "grey") +
      tm_shape(internal_borders) +
      tm_lines(lwd = 1, col = "black") +
      tm_shape(lake)+
      tm_fill("lightblue")+
      tm_borders("darkblue", lwd = 0.2)+
      tm_shape(rostigraben_within_canton) +
      tm_lines(lwd = 2, col = "red") +
      tm_layout(frame = FALSE) 
    print(map)
    tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/10km_language_borders.pdf", width = 11.69, height = 11.69)
    
    # Bins
    
      # Create combined categories: distance bin + agglo flag
      mun$rosti_distance_group <- case_when(
        mun$dummy_5km_rosti == 1 ~ "0–5 km",
        mun$dummy_5km_rosti == 0 & mun$dummy_10km_rosti == 1 ~ "5–10 km",
        TRUE ~ "More than 10 km"
      )
      
      # Factor to control legend order
      mun$rosti_distance_group <- factor(
        mun$rosti_distance_group,
        levels = c("0–5 km", "5–10 km", "More than 10 km")
      )
      
      # Map with different colors per group
      map_rosti <- tm_shape(mun) +
        tm_fill("rosti_distance_group",
                palette = c("tomato", "darksalmon", "white"),
                title = "Distance to language border") +
        tm_borders(lwd = 0.3, col = "grey") +
        tm_shape(internal_borders) +
        tm_lines(lwd = 0.5, col = "black") +
        tm_shape(lake)+
        tm_fill("lightblue")+
        tm_borders("darkblue", lwd = 0.2)+
        tm_shape(rostigraben_within_canton) +
        tm_lines(lwd = 2, col = "black") +
        tm_layout(frame = FALSE)
      
      print(map_rosti)
      
      # Save
      tmap_save(
        map_rosti,
        "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/bins_to_rosti_border.pdf",
        width = 11.69, height = 11.69
      )
      
################################################################################
### Subsamples 
      ### 0-15km
      ##  Direct edge counts  (share within/across cantons)
      
      # Empty list to hold results
      for (year in years) {
        g <- get(paste0("g_firms_lcc_", year))
        
        # Induce subgraph: only firms in agglomeration municipalities
        g_sub <- induced_subgraph(g, vids = which(V(g)$dist_to_canton_border <= 10000))
        
        # Save back with a new name
        assign(paste0("g_firms_lcc_10km_", year), g_sub)
        
        # Optional: print quick stats
        cat("Year:", year,
            "- nodes:", vcount(g_sub),
            "- edges:", ecount(g_sub), "\n")
      }
      
      cross_summary <- list()
      
      for (year in years) {
        g <- get(paste0("g_firms_lcc_10km_", year))   # call graph
        
        if (ecount(g) == 0 || vcount(g) == 0) {
          if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
          cross_summary[[as.character(year)]] <- data.frame(
            year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
          )
          assign(paste0("g_firms_lcc_10km_", year), g)
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
        assign(paste0("g_firms_lcc_15km_", year), g)
        
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
        select(year, share_cross) %>%
        pivot_longer(-year, names_to = "measure", values_to = "value")
      
      p_both <- ggplot(share_long, aes(x = year, y = value, linetype = measure)) +
        geom_line(na.rm = TRUE) +
        geom_point(na.rm = TRUE) +
        scale_linetype_manual(values = c("solid"),
                              labels = c("Share cross-canton")) +
        labs(x = "Year", y = "Share of edges",
             title = "Share of cross-canton edges over time",
             linetype = NULL) +
        theme_minimal(base_size = 14)
      
      print(p_both)
      ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/bins/fig_within_vs_cross_10km.pdf", p_both, width = 11.69, height = 11.69)
      
      
  ### 0-15km
      ##  Direct edge counts  (share within/across cantons)
      
      # Empty list to hold results
      for (year in years) {
        g <- get(paste0("g_firms_lcc_", year))
        
        # Induce subgraph: only firms in agglomeration municipalities
        g_sub <- induced_subgraph(g, vids = which(V(g)$dist_to_canton_border <= 15000))
        
        # Save back with a new name
        assign(paste0("g_firms_lcc_15km_", year), g_sub)
        
        # Optional: print quick stats
        cat("Year:", year,
            "- nodes:", vcount(g_sub),
            "- edges:", ecount(g_sub), "\n")
      }
      
      cross_summary <- list()
      
      for (year in years) {
        g <- get(paste0("g_firms_lcc_15km_", year))   # call graph
        
        if (ecount(g) == 0 || vcount(g) == 0) {
          if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
          cross_summary[[as.character(year)]] <- data.frame(
            year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
          )
          assign(paste0("g_firms_lcc_15km_", year), g)
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
        assign(paste0("g_firms_lcc_15km_", year), g)
        
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
        select(year, share_cross) %>%
        pivot_longer(-year, names_to = "measure", values_to = "value")
      
      p_both <- ggplot(share_long, aes(x = year, y = value, linetype = measure)) +
        geom_line(na.rm = TRUE) +
        geom_point(na.rm = TRUE) +
        scale_linetype_manual(values = c("solid"),
                              labels = c("Share cross-canton")) +
        labs(x = "Year", y = "Share of edges",
             title = "Share of cross-canton edges over time",
             linetype = NULL) +
        theme_minimal(base_size = 14)
      
      print(p_both)
      ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/bins/fig_within_vs_cross_15km.pdf", p_both, width = 11.69, height = 11.69)

  ### Splint by bins            
   
      node_bin_series <- list()
      
      # helper to interpret your 0/1 dummies robustly (numeric, char, or factor)
      is_one <- function(x) {
        if (is.factor(x) || is.character(x)) {
          as.character(x) == "1"
        } else if (is.logical(x)) {
          x
        } else {
          !is.na(x) & as.numeric(x) == 1
        }
      }
      
      for (year in years) {
        g_sub <- get(paste0("g_firms_lcc_15km_", year))  # ≤ 15 km subgraph
        
        if (vcount(g_sub) == 0 || ecount(g_sub) == 0) {
          node_bin_series[[as.character(year)]] <- data.frame(
            year = year, bin = c("0–5 km","5–10 km","10–15 km"),
            avg_share = NA_real_, n_nodes = 0L
          )
          next
        }
        
        # --- bin membership from dummies on vertices ---
        d0  <- V(g_sub)$dummy_0_5km
        d5  <- V(g_sub)$dummy_5_10km
        d10 <- V(g_sub)$dummy_10_15km
        
        V(g_sub)$bin <- dplyr::case_when(
          as.numeric(as.character(V(g_sub)$dummy_0_5km)) == 1 ~ "0–5 km",
          as.numeric(as.character(V(g_sub)$dummy_5_10km)) == 1 ~ "5–10 km",
          as.numeric(as.character(V(g_sub)$dummy_10_15km)) == 1 ~ "10–15 km",
          TRUE ~ NA_character_
        )
        
        V(g_sub)$bin <- factor(V(g_sub)$bin, levels = c("0–5 km","5–10 km","10–15 km"))
        
        # --- node-level share of cross-canton neighbors ---
        ctn <- V(g_sub)$ctn  # assumed already cleaned upstream
        
        node_share <- sapply(seq_len(vcount(g_sub)), function(i) {
          # skip if own canton missing
          if (is.na(ctn[i]) || ctn[i] == "") return(NA_real_)
          nbrs <- neighbors(g_sub, i, mode = "all")
          degi <- length(nbrs)
          if (degi == 0) return(NA_real_)  # isolate
          nbr_idx  <- as.integer(nbrs)
          nbr_ctn  <- ctn[nbr_idx]
          valid    <- !is.na(nbr_ctn) & nbr_ctn != ""
          if (!any(valid)) return(NA_real_)
          mean(nbr_ctn[valid] != ctn[i])
        })
        
        df_nodes <- data.frame(
          bin   = V(g_sub)$bin,
          share = node_share
        )
        
        # average per bin over nodes (drop NAs from isolates / missing canton)
        bin_stats <- df_nodes %>%
          filter(!is.na(bin)) %>%
          group_by(bin) %>%
          summarise(
            avg_share = mean(share, na.rm = TRUE),
            n_nodes   = sum(!is.na(share)),
            .groups   = "drop"
          ) %>%
          mutate(year = year)
        
        # ensure all three bins are present each year (even if empty)
        bin_stats <- tidyr::complete(
          bin_stats,
          bin = factor(c("0–5 km","5–10 km","10–15 km"), levels = c("0–5 km","5–10 km","10–15 km")),
          year = year,
          fill = list(avg_share = NA_real_, n_nodes = 0L)
        )
        
        node_bin_series[[as.character(year)]] <- bin_stats
      }
      
      # Combine across years
      node_bin_df <- bind_rows(node_bin_series)
      
      # Plot: average node-level share of cross-canton neighbors per bin
      p_nodes <- ggplot(node_bin_df, aes(x = year, y = avg_share, color = bin)) +
        geom_line(linewidth = 1, na.rm = TRUE) +
        geom_point(size = 2, na.rm = TRUE) +
        scale_color_manual(values = c("0–5 km" = "firebrick", "5–10 km" = "orange", "10–15 km" = "salmon")) +
        labs(x = "Year",
             y = "Avg. node share of cross-canton neighbors",
             title = "Share of cross-cantonal edges by distance-to-border bin",
             color = "Bin") +
        theme_minimal(base_size = 14)
      
      print(p_nodes)
      ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/bins/fig_within_vs_cross_3bins.pdf", p_nodes, width = 11.69, height = 11.69)
      
### 0-15km / no agglo
##  Direct edge counts  (share within/across cantons)
      
      # Empty list to hold results
      for (year in years) {
        g <- get(paste0("g_firms_lcc_", year))
        
        # Induce subgraph: only firms not in agglomeration municipalities
        g_sub <- induced_subgraph(
          g,
          vids = which(V(g)$dist_to_canton_border <= 15000 & V(g)$dummy_agglo == 0)
        )        
        # Save back with a new name
        assign(paste0("g_firms_lcc_15km_no_agglo", year), g_sub)
        
        # Optional: print quick stats
        cat("Year:", year,
            "- nodes:", vcount(g_sub),
            "- edges:", ecount(g_sub), "\n")
      }
      
      cross_summary <- list()
      
      for (year in years) {
        g <- get(paste0("g_firms_lcc_15km_no_agglo", year))   # call graph
        
        if (ecount(g) == 0 || vcount(g) == 0) {
          if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
          cross_summary[[as.character(year)]] <- data.frame(
            year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
          )
          assign(paste0("g_firms_lcc_15km_no_agglo", year), g)
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
        assign(paste0("g_firms_lcc_15km_no_agglo", year), g)
        
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
        select(year, share_cross) %>%
        pivot_longer(-year, names_to = "measure", values_to = "value")
      
      p_both <- ggplot(share_long, aes(x = year, y = value, linetype = measure)) +
        geom_line(na.rm = TRUE) +
        geom_point(na.rm = TRUE) +
        scale_linetype_manual(values = c("solid"),
                              labels = c("Share cross-canton")) +
        labs(x = "Year", y = "Share of edges",
             title = "Share of cross-canton edges over time",
             linetype = NULL) +
        theme_minimal(base_size = 14)
      
      print(p_both)
      ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/bins/fig_within_vs_cross_15km.pdf", p_both, width = 11.69, height = 11.69)
      
      ### Splint by bins            
      
      node_bin_series <- list()
      
      # helper to interpret your 0/1 dummies robustly (numeric, char, or factor)
      is_one <- function(x) {
        if (is.factor(x) || is.character(x)) {
          as.character(x) == "1"
        } else if (is.logical(x)) {
          x
        } else {
          !is.na(x) & as.numeric(x) == 1
        }
      }
      
      for (year in years) {
        g_sub <- get(paste0("g_firms_lcc_15km_no_agglo", year))  # ≤ 15 km subgraph
        
        if (vcount(g_sub) == 0 || ecount(g_sub) == 0) {
          node_bin_series[[as.character(year)]] <- data.frame(
            year = year, bin = c("0–5 km","5–10 km","10–15 km"),
            avg_share = NA_real_, n_nodes = 0L
          )
          next
        }
        
        # --- bin membership from dummies on vertices ---
        d0  <- V(g_sub)$dummy_0_5km
        d5  <- V(g_sub)$dummy_5_10km
        d10 <- V(g_sub)$dummy_10_15km
        
        V(g_sub)$bin <- dplyr::case_when(
          as.numeric(as.character(V(g_sub)$dummy_0_5km)) == 1 ~ "0–5 km",
          as.numeric(as.character(V(g_sub)$dummy_5_10km)) == 1 ~ "5–10 km",
          as.numeric(as.character(V(g_sub)$dummy_10_15km)) == 1 ~ "10–15 km",
          TRUE ~ NA_character_
        )
        
        V(g_sub)$bin <- factor(V(g_sub)$bin, levels = c("0–5 km","5–10 km","10–15 km"))
        
        # --- node-level share of cross-canton neighbors ---
        ctn <- V(g_sub)$ctn  # assumed already cleaned upstream
        
        node_share <- sapply(seq_len(vcount(g_sub)), function(i) {
          # skip if own canton missing
          if (is.na(ctn[i]) || ctn[i] == "") return(NA_real_)
          nbrs <- neighbors(g_sub, i, mode = "all")
          degi <- length(nbrs)
          if (degi == 0) return(NA_real_)  # isolate
          nbr_idx  <- as.integer(nbrs)
          nbr_ctn  <- ctn[nbr_idx]
          valid    <- !is.na(nbr_ctn) & nbr_ctn != ""
          if (!any(valid)) return(NA_real_)
          mean(nbr_ctn[valid] != ctn[i])
        })
        
        df_nodes <- data.frame(
          bin   = V(g_sub)$bin,
          share = node_share
        )
        
        # average per bin over nodes (drop NAs from isolates / missing canton)
        bin_stats <- df_nodes %>%
          filter(!is.na(bin)) %>%
          group_by(bin) %>%
          summarise(
            avg_share = mean(share, na.rm = TRUE),
            n_nodes   = sum(!is.na(share)),
            .groups   = "drop"
          ) %>%
          mutate(year = year)
        
        # ensure all three bins are present each year (even if empty)
        bin_stats <- tidyr::complete(
          bin_stats,
          bin = factor(c("0–5 km","5–10 km","10–15 km"), levels = c("0–5 km","5–10 km","10–15 km")),
          year = year,
          fill = list(avg_share = NA_real_, n_nodes = 0L)
        )
        
        node_bin_series[[as.character(year)]] <- bin_stats
      }
      
      # Combine across years
      node_bin_df <- bind_rows(node_bin_series)
      
      # Plot: average node-level share of cross-canton neighbors per bin
      p_nodes <- ggplot(node_bin_df, aes(x = year, y = avg_share, color = bin)) +
        geom_line(linewidth = 1, na.rm = TRUE) +
        geom_point(size = 2, na.rm = TRUE) +
        scale_color_manual(values = c("0–5 km" = "firebrick", "5–10 km" = "orange", "10–15 km" = "salmon")) +
        labs(x = "Year",
             y = "Avg. node share of cross-canton neighbors",
             title = "Share of cross-cantonal edges by distance-to-border bin",
             color = "Bin") +
        theme_minimal(base_size = 14)
      
      print(p_nodes)
      ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/bins/fig_within_vs_cross_3bins_noagglo.pdf", p_nodes, width = 11.69, height = 11.69)
      
################################################################################