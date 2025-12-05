library(xtable)
library(openxlsx)
library(haven)
library(dplyr)
library(igraph)
library(readstata13)
#library(xlsx)
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
#years <- c(2000, 2001, 2002, 2003)
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
  
  # Compute macro-level stats (unweighted)
  net_sum <- list(
    year = year,
    number_of_edges = gsize(g_firm),
    number_of_nodes = gorder(g_firm),
    size_LCC = max(components(g_firm)$csize),
    share_LCC = max(components(g_firm)$csize) / gorder(g_firm),
    size_nonLCC = gorder(g_firm) - max(components(g_firm)$csize),
    share_isolates = sum(degree(g_firm) == 0) / gorder(g_firm),
    density = edge_density(g_firm),
    average_degree = mean(degree(g_firm)),
    max_degree = max(degree(g_firm)),
    is_connected = is_connected(g_firm),
    number_of_components = components(g_firm)$no
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
    density = round(density, 5),
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
      file = "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/firm_LCC_summary.tex",
      booktabs = TRUE)

################################################################################
################################################################################
### Firm-to-firm : node-level measures

# Create storage lists
all_full_ranks <- list()
all_top25 <- list()
all_avg_rank <- list()

# Output directory
outdir <- "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Centrality measures - Loop over all years    
for (year in years) {
  giant_g <- get(paste0("g_firms_lcc_", year))
  
  # Centralities
  V(giant_g)$degree      <- degree(giant_g)
  V(giant_g)$betweenness <- betweenness(giant_g, normalized = TRUE)
  V(giant_g)$closeness   <- closeness(giant_g, normalized = TRUE)
  V(giant_g)$eigenvector <- eigen_centrality(giant_g)$vector
  V(giant_g)$pagerank    <- page_rank(giant_g)$vector
  
  centrality_df <- data.frame(
    firm_id     = V(giant_g)$name,
    firm_name   = V(giant_g)$cname,
    degree      = V(giant_g)$degree,
    betweenness = V(giant_g)$betweenness,
    closeness   = V(giant_g)$closeness,
    eigenvector = V(giant_g)$eigenvector,
    pagerank    = V(giant_g)$pagerank,
    stringsAsFactors = FALSE
  )
  
  # Rank all firms by a measure
  make_full_rank <- function(var) {
    centrality_df %>%
      arrange(desc(.data[[var]])) %>%
      mutate(rank = row_number(),
             centrality = var,
             year = year)
  }
  
  # Stack all measures, with full rankings
  full_rank_df <- bind_rows(
    make_full_rank("degree"),
    make_full_rank("betweenness"),
    make_full_rank("closeness"),
    make_full_rank("eigenvector"),
    make_full_rank("pagerank")
  )
  
  # Top 25 per measure
  top25_df <- full_rank_df %>%
    group_by(year, centrality) %>%
    slice_head(n = 25) %>%
    ungroup()
  
  # Average rank across measures (all firms)
  avg_rank_df <- full_rank_df %>%
    group_by(year, firm_id, firm_name) %>%
    summarise(
      avg_rank = mean(rank, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(year, avg_rank)
  
  # Top 25 by average rank
  top25_avg_rank <- avg_rank_df %>%
    group_by(year) %>%
    slice_head(n = 25) %>%
    ungroup()
  
  # Extract top 25 per measure in wide format
  make_top25_wide <- function(var) {
    full_rank_df %>%
      filter(centrality == var) %>%
      slice_head(n = 25) %>%
      select(rank, firm_name, value = .data[[var]])   # keep value too
  }
  
  top25_degree_wide      <- make_top25_wide("degree")
  top25_betweenness_wide <- make_top25_wide("betweenness")
  top25_closeness_wide   <- make_top25_wide("closeness")
  top25_eigenvector_wide <- make_top25_wide("eigenvector")
  top25_pagerank_wide    <- make_top25_wide("pagerank")
  
  # Store
  all_full_ranks[[as.character(year)]] <- full_rank_df
  all_top25[[as.character(year)]]      <- top25_df
  all_avg_rank[[as.character(year)]]   <- avg_rank_df
  
  # Export LaTeX for top25 per measure
  print(xtable(top25_df,
               caption = paste("Top 25 firms by centrality measures,", year),
               label = paste0("tab:top25_all_", year)),
        include.rownames = FALSE,
        file = file.path(outdir, paste0("top25_all_", year, ".tex")),
        booktabs = TRUE)
  
  # Export LaTeX for top25 by average rank
  print(xtable(top25_avg_rank,
               caption = paste("Top 25 firms by average centrality rank,", year),
               label = paste0("tab:top25_avg_", year)),
        include.rownames = FALSE,
        file = file.path(outdir, paste0("top25_avg_", year, ".tex")),
        booktabs = TRUE)
  
  # Top 25 in Degree    
  print(xtable(top25_degree_wide,
               caption = paste("Top 25 firms by degree centrality,", year),
               label   = paste0("tab:top25_degree_", year),
               digits  = 3),
        include.rownames = FALSE,
        file = file.path(outdir, paste0("top25_degree_", year, ".tex")),
        booktabs = TRUE)
  
  # Top 25 in Betweeness
  print(xtable(top25_betweenness_wide,
               caption = paste("Top 25 firms by betweenness centrality,", year),
               label   = paste0("tab:top25_betweenness_", year),
               digits  = 3),
        include.rownames = FALSE,
        file = file.path(outdir, paste0("top25_betweenness_", year, ".tex")),
        booktabs = TRUE)
  
  # Top 25 in Page Rank
  print(xtable(top25_pagerank_wide,
               caption = paste("Top 25 firms by Page Rank centrality,", year),
               label   = paste0("tab:top25_pagerank_", year),
               digits  = 3),
        include.rownames = FALSE,
        file = file.path(outdir, paste0("top25_pagerank_", year, ".tex")),
        booktabs = TRUE) 
  
  # Top 25 in Eigenvector Centrality
  print(xtable(top25_eigenvector_wide,
               caption = paste("Top 25 firms by eigenvector centrality,", year),
               label   = paste0("tab:top25_eigenvector_", year),
               digits  = 3),
        include.rownames = FALSE,
        file = file.path(outdir, paste0("top25_eigenvector_", year, ".tex")),
        booktabs = TRUE)  
  
  # Top 25 in Closeness
  print(xtable(top25_closeness_wide,
               caption = paste("Top 25 firms by closeness centrality,", year),
               label   = paste0("tab:top25_closeness_", year),
               digits  = 3),
        include.rownames = FALSE,
        file = file.path(outdir, paste0("top25_closeness_", year, ".tex")),
        booktabs = TRUE)         
}

# Combine across years
all_full_ranks_df <- bind_rows(all_full_ranks)
all_top25_df      <- bind_rows(all_top25)
all_avg_rank_df   <- bind_rows(all_avg_rank)

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
    within = within_cnt,
    cross  = cross_cnt,
    share_within = mean_within,
    share_cross  = mean_cross,
    ci95_low_within  = ci95_low_within,
    ci95_high_within = ci95_high_within,
    ci95_low_cross   = ci95_low_cross,
    ci95_high_cross  = ci95_high_cross
  )
}

# Bind the summary into a data.frame for plotting
cross_summary_df <- dplyr::bind_rows(cross_summary)
print(cross_summary_df)

# Reshape for plotting
share_long <- cross_summary_df %>%
  select(year, share_within, share_cross, 
         ci95_low_within, ci95_high_within,
         ci95_low_cross, ci95_high_cross) %>%
  pivot_longer(cols = c(share_within, share_cross),
               names_to = "measure", values_to = "value") %>%
  mutate(
    ci_low  = ifelse(measure == "share_cross", ci95_low_cross, ci95_low_within),
    ci_high = ifelse(measure == "share_cross", ci95_high_cross, ci95_high_within)
  )

# Plot with 95% CI ribbon
p_both <- ggplot(share_long, aes(x = year, y = value, linetype = measure)) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = measure),
              alpha = 0.2, colour = NA) +
  scale_linetype_manual(values = c("solid","dashed"),
                        labels = c("Share cross-canton", "Share within canton")) +
  scale_fill_manual(values = c("grey70","grey40"), guide = "none") +
  labs(x = "Year", y = "Share of edges",
       title = "Within- vs cross-canton edges over time (95% CI)",
       linetype = NULL) +
  theme_minimal(base_size = 14)

print(p_both)
ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_within_vs_cross.pdf", p_both, width = 11.69, height = 11.69)

### Direct weighted edge counts (within/across cantons)

weighted_summary <- list()

for (year in years) {
  g <- get(paste0("g_firms_lcc_", year))  # use your LCC graphs
  
  # ensure cross_canton exists (1 = cross, 0 = within, NA if missing canton) 
  if (is.null(E(g)$cross_canton)) {
    ctn <- V(g)$ctn
    if (ecount(g) > 0 && vcount(g) > 0) {
      ep   <- ends(g, E(g), names = FALSE)
      bad  <- is.na(ctn[ep[,1]]) | is.na(ctn[ep[,2]]) | ctn[ep[,1]] == "" | ctn[ep[,2]] == ""
      diff <- ctn[ep[,1]] != ctn[ep[,2]]
      E(g)$cross_canton <- ifelse(bad, NA_integer_, as.integer(diff))
    } else {
      E(g)$cross_canton <- integer(0)
    }
    assign(paste0("g_firms_lcc_", year), g)
  }
  
  # Weights
  w <- NULL
  if (!is.null(E(g)$weight)) {
    w <- E(g)$weight
  } else if (!is.null(E(g)$n_shared_directors)) {
    w <- E(g)$n_shared_directors
    E(g)$weight <- w  # store for later reuse
    assign(paste0("g_firms_lcc_", year), g)
  } else {
    # no weights defined; treat as 1
    w <- rep(1, ecount(g))
  }
  
  # keep only edges with defined cross_canton and non-missing weights
  if (ecount(g) == 0) {
    weighted_summary[[as.character(year)]] <- data.frame(
      year = year, within_w = 0, cross_w = 0,
      share_within_w = NA_real_, share_cross_w = NA_real_
    )
    next
  }
  cc <- E(g)$cross_canton
  keep <- !is.na(cc) & !is.na(w)
  if (!any(keep)) {
    weighted_summary[[as.character(year)]] <- data.frame(
      year = year, within_w = 0, cross_w = 0,
      share_within_w = NA_real_, share_cross_w = NA_real_
    )
    next
  }
  
  w_keep <- w[keep]
  cc_keep <- cc[keep]
  
  within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
  cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
  tot_w    <- within_w + cross_w
  
  weighted_summary[[as.character(year)]] <- data.frame(
    year = year,
    within_w = within_w,
    cross_w  = cross_w,
    share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
    share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
  )
}

cross_summary_w_df <- bind_rows(weighted_summary)
print(cross_summary_w_df)

# Plot that
share_w_long <- cross_summary_w_df %>%
  select(year, share_within_w, share_cross_w) %>%
  pivot_longer(-year, names_to = "measure", values_to = "value")

p_both_w <- ggplot(share_w_long, aes(x = year, y = value, linetype = measure)) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE) +
  scale_linetype_manual(values = c("solid","dashed"),
                        labels = c("Weighted share cross-canton", "Weighted share within canton")) +
  labs(x = "Year", y = "Weighted share of edges",
       title = "Within vs. cross-canton (weighted)",
       linetype = NULL) +
  theme_minimal(base_size = 14)

print(p_both_w)
ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_within_vs_cross_weighted.pdf", p_both_w, width = 11.69, height = 11.69)

### Direct edge counts  (share within/across languages)

# Storage for summary per year
cross_summary_lang <- list()

for (year in years) {
  g <- get(paste0("g_firms_lcc_", year))   # your LCC graph for that year
  
  if (ecount(g) == 0 || vcount(g) == 0) {
    # no edges or no nodes: set empty attribute + summary
    if (ecount(g) == 0) E(g)$cross_lang <- integer(0)
    cross_summary_lang[[as.character(year)]] <- data.frame(
      year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
    )
    assign(paste0("g_firms_lcc_", year), g)
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
ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_within_vs_cross_lang.pdf", p_both_lang, width = 11.69, height = 11.69)


### Direct weighted edge counts (within/across lang)

weighted_summary_lang <- list()

for (year in years) {
  g <- get(paste0("g_firms_lcc_", year))  # use your LCC graphs
  
  # ensure cross_lang exists (1 = cross, 0 = within, NA if missing lang) 
  if (is.null(E(g)$cross_lang)) {
    lang <- V(g)$language
    if (ecount(g) > 0 && vcount(g) > 0) {
      ep   <- ends(g, E(g), names = FALSE)
      bad  <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]] == "" | lang[ep[,2]] == ""
      diff <- lang[ep[,1]] != lang[ep[,2]]
      E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
    } else {
      E(g)$cross_lang <- integer(0)
    }
    assign(paste0("g_firms_lcc_", year), g)
  }
  
  # Weights
  w <- NULL
  if (!is.null(E(g)$weight)) {
    w <- E(g)$weight
  } else if (!is.null(E(g)$n_shared_directors)) {
    w <- E(g)$n_shared_directors
    E(g)$weight <- w  # store for later reuse
    assign(paste0("g_firms_lcc_", year), g)
  } else {
    # no weights defined; treat as 1
    w <- rep(1, ecount(g))
  }
  
  # keep only edges with defined cross_lang and non-missing weights
  if (ecount(g) == 0) {
    weighted_summary_lang[[as.character(year)]] <- data.frame(
      year = year, within_w = 0, cross_w = 0,
      share_within_w = NA_real_, share_cross_w = NA_real_
    )
    next
  }
  cc <- E(g)$cross_lang
  keep <- !is.na(cc) & !is.na(w)
  if (!any(keep)) {
    weighted_summary_lang[[as.character(year)]] <- data.frame(
      year = year, within_w = 0, cross_w = 0,
      share_within_w = NA_real_, share_cross_w = NA_real_
    )
    next
  }
  
  w_keep <- w[keep]
  cc_keep <- cc[keep]
  
  within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
  cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
  tot_w    <- within_w + cross_w
  
  weighted_summary_lang[[as.character(year)]] <- data.frame(
    year = year,
    within_w = within_w,
    cross_w  = cross_w,
    share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
    share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
  )
}

cross_summary_lang_w_df <- bind_rows(weighted_summary_lang)
print(cross_summary_lang_w_df)

# Plot that
share_w_long_lang <- cross_summary_lang_w_df %>%
  select(year, share_within_w, share_cross_w) %>%
  pivot_longer(-year, names_to = "measure", values_to = "value")

p_both_w <- ggplot(share_w_long_lang, aes(x = year, y = value, linetype = measure)) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE) +
  scale_linetype_manual(values = c("solid","dashed"),
                        labels = c("Weighted share cross-language", "Weighted share within language")) +
  labs(x = "Year", y = "Weighted share of edges",
       title = "Within vs. cross-language (weighted)",
       linetype = NULL) +
  theme_minimal(base_size = 14)

print(p_both_w)
ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_within_vs_cross_lang_weighted.pdf", p_both_w, width = 11.69, height = 11.69)        

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

# Weighted
w_combined <- bind_rows(
  share_w_long %>%
    mutate(family = "canton"),
  share_w_long_lang %>%
    mutate(family = "language")
) %>%
  mutate(
    kind = ifelse(grepl("cross", measure), "cross", "within"),
    color_key = ifelse(family == "canton", "canton", "language"),
    lty_key   = ifelse(kind == "cross", "solid", "dashed")
  )

p_w <- ggplot(w_combined, aes(x = year, y = value,
                              group = interaction(family, kind))) +
  geom_line(aes(color = color_key, linetype = lty_key), na.rm = TRUE) +
  geom_point(aes(color = color_key), na.rm = TRUE, size = 1.8) +
  scale_color_manual(values = c(canton = "red", language = "black"),
                     labels = c("Canton", "Language region"), name = NULL) +
  scale_linetype_identity(guide = "none") +
  labs(x = "Year", y = "Weighted share of edges",
       title = "Within vs. cross edges — Canton (red) vs Language (black), weighted",
       subtitle = "Solid = cross; Dashed = within") +
  theme_minimal(base_size = 13)

print(p_w)
ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/fig_within_cross_canton_language_weighted.pdf", p_w, width = 11.69, height = 11.69)

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
  scale_color_manual(values = c("Canton" = "red", "Language" = "black")) +
  labs(x = "Year", y = "Assortativity",
       title = "Assortativity of firm-to-firm networks",
       subtitle = "Canton vs. language") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())

print(assortativity_graph)

ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/assortativity_graph.pdf", assortativity_graph, width = 11.69, height = 11.69)


################################################################################
################################################################################            

### Node-level homophily      
# Note: Assortativity is a global correlation measure, while node-level homophily is the average of local shares (node-level)

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
################################################################################  
### Ajusted Rand Index : how similar the detected communities are to observed cantons (or languages)
# Cantons

ari_vs_algorithms_canton <- function(g, weighted = FALSE) {
  if (vcount(g) == 0 || ecount(g) == 0) return(rep(NA_real_, 4))
  
  lab <- V(g)$ctn
  keep <- which(!is.na(lab) & lab != "")
  if (length(keep) == 0) return(rep(NA_real_, 4))
  
  g2 <- induced_subgraph(g, keep)
  lab2 <- V(g2)$ctn
  lab_num <- as.integer(factor(lab2))
  
  # Weights
  w <- if (weighted) {
    if (!is.null(E(g2)$weight)) E(g2)$weight
    else if (!is.null(E(g2)$n_shared_directors)) E(g2)$n_shared_directors
    else NULL
  } else NULL
  
  # Communities
  memb_louv <- as.integer(membership(if (is.null(w)) cluster_louvain(g2) else cluster_louvain(g2, weights = w)))
  memb_leid <- as.integer(membership(if (is.null(w)) cluster_leiden(g2) else cluster_leiden(g2, weights = w)))
  memb_walk <- as.integer(membership(if (is.null(w)) cluster_walktrap(g2) else cluster_walktrap(g2, weights = w)))
  memb_fast <- as.integer(membership(if (is.null(w)) cluster_fast_greedy(g2) else cluster_fast_greedy(g2, weights = w)))
  
  # ARIs
  c(
    Louvain    = adjustedRandIndex(memb_louv, lab_num),
    Leiden     = adjustedRandIndex(memb_leid, lab_num),
    Walktrap   = adjustedRandIndex(memb_walk, lab_num),
    FastGreedy = adjustedRandIndex(memb_fast, lab_num)
  )
}

# Languages 
ari_vs_algorithms_language <- function(g, weighted = FALSE) {
  if (vcount(g) == 0 || ecount(g) == 0) return(rep(NA_real_, 4))
  
  lab <- V(g)$language
  keep <- which(!is.na(lab) & lab != "")
  if (length(keep) == 0) return(rep(NA_real_, 4))
  
  g2 <- induced_subgraph(g, keep)
  lab2 <- V(g2)$language
  lab_num <- as.integer(factor(lab2))
  
  # Weights
  w <- if (weighted) {
    if (!is.null(E(g2)$weight)) E(g2)$weight
    else if (!is.null(E(g2)$n_shared_directors)) E(g2)$n_shared_directors
    else NULL
  } else NULL
  
  # Communities
  memb_louv <- as.integer(membership(if (is.null(w)) cluster_louvain(g2) else cluster_louvain(g2, weights = w)))
  memb_leid <- as.integer(membership(if (is.null(w)) cluster_leiden(g2) else cluster_leiden(g2, weights = w)))
  memb_walk <- as.integer(membership(if (is.null(w)) cluster_walktrap(g2) else cluster_walktrap(g2, weights = w)))
  memb_fast <- as.integer(membership(if (is.null(w)) cluster_fast_greedy(g2) else cluster_fast_greedy(g2, weights = w)))
  
  # ARIs
  c(
    Louvain    = adjustedRandIndex(memb_louv, lab_num),
    Leiden     = adjustedRandIndex(memb_leid, lab_num),
    Walktrap   = adjustedRandIndex(memb_walk, lab_num),
    FastGreedy = adjustedRandIndex(memb_fast, lab_num)
  )
}
# Loop over years      

ari_multi_results <- list()

for (year in years) {
  g <- get(paste0("g_firms_lcc_", year))
  
  ari_ctn <- ari_vs_algorithms_canton(g, weighted = TRUE)
  ari_lang <- ari_vs_algorithms_language(g, weighted = TRUE)
  
  ari_multi_results[[as.character(year)]] <- data.frame(
    year = year,
    partition = rep(c("Canton", "Language"), each = 4),
    algorithm = rep(c("Louvain", "Leiden", "Walktrap", "FastGreedy"), times = 2),
    ARI = c(ari_ctn, ari_lang)
  )
}

ari_multi_df <- dplyr::bind_rows(ari_multi_results)
print(ari_multi_df)

# Two plots (one for cantons, one for langage, with the 4 ARI)
ggplot(ari_multi_df, aes(x = year, y = ARI, color = algorithm)) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(~partition, ncol = 1) +
  scale_color_manual(values = c(
    "Louvain"    = "red",
    "Leiden"     = "black",
    "Walktrap"   = "blue",
    "FastGreedy" = "darkgreen"
  )) +
  labs(
    x = "Year",
    y = "Adjusted Rand Index (ARI)",
    title = "Similarity between detected communities and canton/language partitions",
    subtitle = "ARI: 0 = random similarity, 1 = perfect match"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())

################################################################################
################################################################################      
### Louvain community detection of firm-to-firm networks w/o single nodes    

# Drop single nodes
for (year in years) {
  g <- get(paste0("g_firms_", year))       # get graph for that year
  
  # Drop isolates (nodes with degree == 0)
  g_no_isolates <- delete_vertices(g, which(degree(g) == 0))
  
  # Save back
  assign(paste0("g_firms_noisolates_", year), g_no_isolates)
}

# Louvain Communities
for (year in years) {
  g_firm <- get(paste0("g_firms_noisolates_", year))
  if (gorder(g_firm) > 0 && gsize(g_firm) > 0) {
    comm_firm <- cluster_louvain(g_firm, weights = E(g_firm)$weight)
    V(g_firm)$community <- membership(comm_firm)
    assign(paste0("g_firms_", year), g_firm)
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
  g_name <- paste0("g_firms_", year)
  
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

################################################################################    
# Comparing communities and cantonal borders in Network w/o single node

for (year in years) {
  g <- get(paste0("g_firms_", year))
  
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
  g_name <- paste0("g_firms_", year)
  g <- get(g_name)
  
  # Ensure both attributes are present
  if (!is.null(V(g)$community) && !is.null(V(g)$ctn)) {
    # Remove nodes with missing canton
    valid <- !is.na(V(g)$community) & !is.na(V(g)$ctn)
    ari_nosingle <- adjustedRandIndex(V(g)$community[valid], V(g)$ctn[valid])
    
    ari_results <- rbind(
      ari_results,
      data.frame(year = year, adjusted_rand_index_nosingle = ari_nosingle)
    )
  }
}

# View ARI values over time
print(ari_results)

# plot
plot(
  ari_results$year,
  ari_results$adjusted_rand_index_nosingle,
  type = "b", pch = 16,
  xlab = "Year", ylab = "Adjusted Rand Index",
  main = "Similarity Between Firm Communities and Cantons (w/o single nodes)"
)

################################################################################    
# Comparing communities and langage borders in Network w/o single node   

for (year in years) {
  g <- get(paste0("g_firms_", year))
  
  # Only if both attributes are present
  if (!is.null(V(g)$community) && !is.null(V(g)$language)) {
    tab_firm_language <- table(V(g)$community, V(g)$language)
    assign(paste0("tab_firm_language_", year), tab_firm_language)
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
  g_name <- paste0("g_firms_", year)
  g <- get(g_name)
  
  # Ensure both attributes are present
  if (!is.null(V(g)$community) && !is.null(V(g)$language)) {
    # Remove nodes with missing canton
    valid <- !is.na(V(g)$community) & !is.na(V(g)$language)
    ari_nosingle <- adjustedRandIndex(V(g)$community[valid], V(g)$language[valid])
    
    ari_results <- rbind(
      ari_results,
      data.frame(year = year, adjusted_rand_index_nosingle = ari_nosingle)
    )
  }
}

# View ARI values over time
print(ari_results)

# plot
plot(
  ari_results$year,
  ari_results$adjusted_rand_index_nosingle,
  type = "b", pch = 16,
  xlab = "Year", ylab = "Adjusted Rand Index",
  main = "Similarity Between Firm Communities and Languages (w/o single nodes)"
)
######    
### Louvain community detection on LCC of firm-to-firm networks    

# Same w/ LCC
for (year in years) {
  g_firm <- get(paste0("g_firms_lcc_", year))
  if (gorder(g_firm) > 0 && gsize(g_firm) > 0) {
    comm_firm <- cluster_louvain(g_firm, weights = E(g_firm)$weight)
    V(g_firm)$community <- membership(comm_firm)
    assign(paste0("g_firms_", year), g_firm)
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
  g_name <- paste0("g_firms_", year)
  
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

######################################################################  
# Comparing communities and cantonal borders in LCC

for (year in years) {
  g <- get(paste0("g_firms_", year))
  
  # Only if both attributes are present
  if (!is.null(V(g)$community) && !is.null(V(g)$ctn)) {
    tab_firm_canton <- table(V(g)$community, V(g)$ctn)
    assign(paste0("tab_firm_canton_", year), tab_firm_canton)
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
  g_name <- paste0("g_firms_", year)
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

######################################################################    
# Comparing communities and langage borders in LCC    

for (year in years) {
  g <- get(paste0("g_firms_", year))
  
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
  g_name <- paste0("g_firms_", year)
  g <- get(g_name)
  
  # Ensure both attributes are present
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

################################################################################
################################################################################

### Other community detection algorithm and ARI plots



################################################################################
################################################################################
### Subsample of mun close to borders
### Maps
#  Maps with shade of colours
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

# 10km map
map <- tm_shape(mun) +
  tm_fill("dummy_10km", 
          palette = c("white", "darksalmon"), 
          title = "Within 10 km of border") +
  tm_borders(lwd = 0.3, col = "grey") +
  tm_shape(internal_borders) +
  tm_lines(lwd = 1, col = "black") +
  tm_shape(lake)+
  tm_fill("lightblue")+
  tm_borders("darkblue", lwd = 0.2)+
  tm_layout(frame = FALSE)
print(map)
tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/10km_to_canton_borders.pdf", width = 11.69, height = 11.69)

# 5km map
map <- tm_shape(mun) +
  tm_fill("dummy_5km", 
          palette = c("white", "darksalmon"), 
          title = "Within 5 km of border") +
  tm_borders(lwd = 0.3, col = "grey") +
  tm_shape(internal_borders) +
  tm_lines(lwd = 1, col = "black") +
  tm_shape(lake)+
  tm_fill("lightblue")+
  tm_borders("darkblue", lwd = 0.2)+
  tm_layout(frame = FALSE)        
print(map)
tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/5km_to_canton_borders.pdf", width = 11.69, height = 11.69)

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

# Touching borders map
map <- tm_shape(mun) +
  tm_fill("dummy_touch_border", 
          palette = c("white", "darksalmon"), 
          title = "Touching cantonal borders") +
  tm_borders(lwd = 0.3, col = "grey") +
  tm_shape(internal_borders) +
  tm_lines(lwd = 1, col = "black") +
  tm_shape(lake)+
  tm_fill("lightblue")+
  tm_borders("darkblue", lwd = 0.2)+
  tm_layout(frame = FALSE)
print(map)
tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/touch_canton_borders.pdf", width = 11.69, height = 11.69)

#Same, without agglo

map <- tm_shape(mun) +
  tm_fill("dummy_10km_noagglo", 
          palette = c("white", "darksalmon"), 
          title = "Within 10 km of border, not in agglo") +
  tm_borders(lwd = 0.3, col = "grey") +
  tm_shape(internal_borders) +
  tm_lines(lwd = 1, col = "black") +
  tm_shape(lake)+
  tm_fill("lightblue")+
  tm_borders("darkblue", lwd = 0.2)+
  tm_layout(frame = FALSE) 
print(map)
tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/10km_to_canton_borders_noagglo.pdf", width = 11.69, height = 11.69)

map <- tm_shape(mun) +
  tm_fill("dummy_5km_noagglo", 
          palette = c("white", "darksalmon"), 
          title = "Within 5 km of border, not in agglo") +
  tm_borders(lwd = 0.3, col = "grey") +
  tm_shape(internal_borders) +
  tm_lines(lwd = 1, col = "black") +
  tm_shape(lake)+
  tm_fill("lightblue")+
  tm_borders("darkblue", lwd = 0.2)+
  tm_layout(frame = FALSE) 
print(map)
tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/5km_to_canton_borders_noagglo.pdf", width = 11.69, height = 11.69)


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

map <- tm_shape(mun) +
  tm_fill("dummy_touch_border_noagglo", 
          palette = c("white", "darksalmon"), 
          title = "Touching cantonal borders, not in agglo") +
  tm_borders(lwd = 0.3, col = "grey") +
  tm_shape(internal_borders) +
  tm_lines(lwd = 1, col = "black") +
  tm_shape(lake)+
  tm_fill("lightblue")+
  tm_borders("darkblue", lwd = 0.2)+
  tm_layout(frame = FALSE)
print(map)
tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/touch_canton_borders_noagglo.pdf", width = 11.69, height = 11.69)

#Non border / Agglo
map <- tm_shape(mun) +
  tm_fill("dummy_nonborder_agglo", 
          palette = c("white", "darksalmon"), 
          title = "Not Touching cantonal borders, in agglo") +
  tm_borders(lwd = 0.3, col = "grey") +
  tm_shape(internal_borders) +
  tm_lines(lwd = 1, col = "black") +
  tm_shape(lake)+
  tm_fill("lightblue")+
  tm_borders("darkblue", lwd = 0.2)+
  tm_layout(frame = FALSE)
print(map)
tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/nonborders_agglo.pdf", width = 11.69, height = 11.69)

#Agglo
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

### Maps with language borders

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

#  Touching
map <- tm_shape(mun) +
  tm_fill("dummy_touch_rosti_border", 
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
tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/touch_language_borders.pdf", width = 11.69, height = 11.69)

#  5km from rosti
map <- tm_shape(mun) +
  tm_fill("dummy_5km_rosti", 
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
tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/5km_language_borders.pdf", width = 11.69, height = 11.69)

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

#  15km from rosti
map <- tm_shape(mun) +
  tm_fill("dummy_15km_rosti", 
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
tmap_save(map, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/subsamples_map/15km_language_borders.pdf", width = 11.69, height = 11.69)


################################################################################
#### Within/Across cantonal borders in different subsamples
### AGGLO
## Direct edge counts  (share within/across cantons)

  # Empty list to hold results
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    
    # Induce subgraph: only firms in agglomeration municipalities
    g_sub <- induced_subgraph(g, vids = which(V(g)$dummy_agglo == 1))
    
    # Save back with a new name
    assign(paste0("g_firms_lcc_agglo_", year), g_sub)
    
    # Optional: print quick stats
    cat("Year:", year,
        "- nodes:", vcount(g_sub),
        "- edges:", ecount(g_sub), "\n")
  }
  
  cross_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_", year))   # call graph
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
      cross_summary[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_agglo_", year), g)
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
  assign(paste0("g_firms_lcc_agglo_", year), g)
  
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
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo/fig_within_vs_cross.pdf", p_both, width = 11.69, height = 11.69)

  ## Direct weighted edge counts (within/across cantons)

  weighted_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_", year))  # use your LCC graphs
    
    # ensure cross_canton exists (1 = cross, 0 = within, NA if missing canton) 
    if (is.null(E(g)$cross_canton)) {
      ctn <- V(g)$ctn
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(ctn[ep[,1]]) | is.na(ctn[ep[,2]]) | ctn[ep[,1]] == "" | ctn[ep[,2]] == ""
        diff <- ctn[ep[,1]] != ctn[ep[,2]]
        E(g)$cross_canton <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_canton <- integer(0)
      }
      assign(paste0("g_firms_lcc_agglo_", year), g)
    }
  
  # Weights
  w <- NULL
  if (!is.null(E(g)$weight)) {
    w <- E(g)$weight
  } else if (!is.null(E(g)$n_shared_directors)) {
    w <- E(g)$n_shared_directors
    E(g)$weight <- w  # store for later reuse
    assign(paste0("g_firms_lcc_agglo_", year), g)
  } else {
    # no weights defined; treat as 1
    w <- rep(1, ecount(g))
  }
  
  # keep only edges with defined cross_canton and non-missing weights
  if (ecount(g) == 0) {
    weighted_summary[[as.character(year)]] <- data.frame(
      year = year, within_w = 0, cross_w = 0,
      share_within_w = NA_real_, share_cross_w = NA_real_
    )
    next
  }
  cc <- E(g)$cross_canton
  keep <- !is.na(cc) & !is.na(w)
  if (!any(keep)) {
    weighted_summary[[as.character(year)]] <- data.frame(
      year = year, within_w = 0, cross_w = 0,
      share_within_w = NA_real_, share_cross_w = NA_real_
    )
    next
  }
  
  w_keep <- w[keep]
  cc_keep <- cc[keep]
  
  within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
  cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
  tot_w    <- within_w + cross_w
  
  weighted_summary[[as.character(year)]] <- data.frame(
    year = year,
    within_w = within_w,
    cross_w  = cross_w,
    share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
    share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }

  cross_summary_w_df <- bind_rows(weighted_summary)
  print(cross_summary_w_df)
  
  # Plot that
  share_w_long <- cross_summary_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-canton")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-canton (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo/fig_within_vs_cross_weighted.pdf", p_both_w, width = 11.69, height = 11.69)

  ## Direct edge counts  (share within/across languages)

  # Storage for summary per year
  cross_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_", year))   # your LCC graph for that year
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      # no edges or no nodes: set empty attribute + summary
      if (ecount(g) == 0) E(g)$cross_lang <- integer(0)
      cross_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_agglo_", year), g)
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
    assign(paste0("g_firms_lcc_agglo_", year), g)
    
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
    select(year, share_cross_lang) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_lang <- ggplot(share_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Share within language region")) +
    labs(x = "Year", y = "Share of edges",
         title = "Across language regions edges over time",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_lang)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo/fig_within_vs_cross_lang.pdf", p_both_lang, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across lang)

  weighted_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_", year))  # use your LCC graphs
    
    # ensure cross_lang exists (1 = cross, 0 = within, NA if missing lang) 
    if (is.null(E(g)$cross_lang)) {
      lang <- V(g)$language
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]] == "" | lang[ep[,2]] == ""
        diff <- lang[ep[,1]] != lang[ep[,2]]
        E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_lang <- integer(0)
      }
      assign(paste0("g_firms_lcc_agglo_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_agglo_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_lang and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_lang
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary_lang[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_lang_w_df <- bind_rows(weighted_summary_lang)
  print(cross_summary_lang_w_df)
  
  # Plot that
  share_w_long_lang <- cross_summary_lang_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-language")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-language (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo/fig_within_vs_cross_lang_weighted.pdf", p_both_w, width = 11.69, height = 11.69)        

  ## Comparing within/across cantons & within/across language regions 
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
         title = "Cross edges — Canton (red) vs Language (black)") +
    theme_minimal(base_size = 13)
  
  print(p_unw)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo/fig_within_cross_canton_language_unweighted.pdf", p_unw, width = 11.69, height = 11.69)
  
  # Weighted
  w_combined <- bind_rows(
    share_w_long %>%
      mutate(family = "canton"),
    share_w_long_lang %>%
      mutate(family = "language")
  ) %>%
    mutate(
      kind = ifelse(grepl("cross", measure), "cross", "within"),
      color_key = ifelse(family == "canton", "canton", "language"),
      lty_key   = ifelse(kind == "cross", "solid", "dashed")
    )
  
  p_w <- ggplot(w_combined, aes(x = year, y = value,
                                group = interaction(family, kind))) +
    geom_line(aes(color = color_key, linetype = lty_key), na.rm = TRUE) +
    geom_point(aes(color = color_key), na.rm = TRUE, size = 1.8) +
    scale_color_manual(values = c(canton = "red", language = "black"),
                       labels = c("Canton", "Language region"), name = NULL) +
    scale_linetype_identity(guide = "none") +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross edges — Canton (red) vs Language (black), weighted") +
    theme_minimal(base_size = 13)
  
  print(p_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo/fig_within_cross_canton_language_weighted.pdf", p_w, width = 11.69, height = 11.69)

################################################################################
  ### NO AGGLO
  ## Direct edge counts  (share within/across cantons)
  
  # Empty list to hold results
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    
    # Induce subgraph: only firms in agglomeration municipalities
    g_sub <- induced_subgraph(g, vids = which(V(g)$dummy_agglo == 0))
    
    # Save back with a new name
    assign(paste0("g_firms_lcc_noagglo_", year), g_sub)
    
    # Optional: print quick stats
    cat("Year:", year,
        "- nodes:", vcount(g_sub),
        "- edges:", ecount(g_sub), "\n")
  }
  
  cross_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_noagglo_", year))   # call graph
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
      cross_summary[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_noagglo_", year), g)
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
    assign(paste0("g_firms_lcc_noagglo_", year), g)
    
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
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/noagglo/fig_within_vs_cross.pdf", p_both, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across cantons)
  
  weighted_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_noagglo_", year))  # use your LCC graphs
    
    # ensure cross_canton exists (1 = cross, 0 = within, NA if missing canton) 
    if (is.null(E(g)$cross_canton)) {
      ctn <- V(g)$ctn
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(ctn[ep[,1]]) | is.na(ctn[ep[,2]]) | ctn[ep[,1]] == "" | ctn[ep[,2]] == ""
        diff <- ctn[ep[,1]] != ctn[ep[,2]]
        E(g)$cross_canton <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_canton <- integer(0)
      }
      assign(paste0("g_firms_lcc_noagglo_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_noagglo_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_canton and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_canton
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_w_df <- bind_rows(weighted_summary)
  print(cross_summary_w_df)
  
  # Plot that
  share_w_long <- cross_summary_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-canton")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-canton (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/noagglo/fig_within_vs_cross_weighted.pdf", p_both_w, width = 11.69, height = 11.69)
  
  ## Direct edge counts  (share within/across languages)
  
  # Storage for summary per year
  cross_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_noagglo_", year))   # your LCC graph for that year
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      # no edges or no nodes: set empty attribute + summary
      if (ecount(g) == 0) E(g)$cross_lang <- integer(0)
      cross_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_noagglo_", year), g)
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
    assign(paste0("g_firms_lcc_noagglo_", year), g)
    
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
    select(year, share_cross_lang) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_lang <- ggplot(share_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Share within language region")) +
    labs(x = "Year", y = "Share of edges",
         title = "Across language regions edges over time",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_lang)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/noagglo/fig_within_vs_cross_lang.pdf", p_both_lang, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across lang)
  
  weighted_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_noagglo_", year))  # use your LCC graphs
    
    # ensure cross_lang exists (1 = cross, 0 = within, NA if missing lang) 
    if (is.null(E(g)$cross_lang)) {
      lang <- V(g)$language
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]] == "" | lang[ep[,2]] == ""
        diff <- lang[ep[,1]] != lang[ep[,2]]
        E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_lang <- integer(0)
      }
      assign(paste0("g_firms_lcc_noagglo_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_noagglo_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_lang and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_lang
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary_lang[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_lang_w_df <- bind_rows(weighted_summary_lang)
  print(cross_summary_lang_w_df)
  
  # Plot that
  share_w_long_lang <- cross_summary_lang_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-language")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-language (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/noagglo/fig_within_vs_cross_lang_weighted.pdf", p_both_w, width = 11.69, height = 11.69)        
  
  ## Comparing within/across cantons & within/across language regions 
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
         title = "Cross edges — Canton (red) vs Language (black)") +
    theme_minimal(base_size = 13)
  
  print(p_unw)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/noagglo/fig_within_cross_canton_language_unweighted.pdf", p_unw, width = 11.69, height = 11.69)
  
  # Weighted
  w_combined <- bind_rows(
    share_w_long %>%
      mutate(family = "canton"),
    share_w_long_lang %>%
      mutate(family = "language")
  ) %>%
    mutate(
      kind = ifelse(grepl("cross", measure), "cross", "within"),
      color_key = ifelse(family == "canton", "canton", "language"),
      lty_key   = ifelse(kind == "cross", "solid", "dashed")
    )
  
  p_w <- ggplot(w_combined, aes(x = year, y = value,
                                group = interaction(family, kind))) +
    geom_line(aes(color = color_key, linetype = lty_key), na.rm = TRUE) +
    geom_point(aes(color = color_key), na.rm = TRUE, size = 1.8) +
    scale_color_manual(values = c(canton = "red", language = "black"),
                       labels = c("Canton", "Language region"), name = NULL) +
    scale_linetype_identity(guide = "none") +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross edges — Canton (red) vs Language (black), weighted") +
    theme_minimal(base_size = 13)
  
  print(p_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/noagglo/fig_within_cross_canton_language_weighted.pdf", p_w, width = 11.69, height = 11.69)  
  
################################################################################
  ### Non Border Agglo
  ## Direct edge counts  (share within/across cantons)
  
  # Empty list to hold results
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    
    # Induce subgraph: only firms in agglomeration municipalities
    g_sub <- induced_subgraph(g, vids = which(V(g)$dummy_agglo == 0))
    
    # Save back with a new name
    assign(paste0("g_firms_lcc_agglo_nonborder_", year), g_sub)
    
    # Optional: print quick stats
    cat("Year:", year,
        "- nodes:", vcount(g_sub),
        "- edges:", ecount(g_sub), "\n")
  }
  
  cross_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_nonborder_", year))   # call graph
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
      cross_summary[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_agglo_nonborder_", year), g)
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
    assign(paste0("g_firms_lcc_agglo_nonborder_", year), g)
    
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
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo_nonborder/fig_within_vs_cross.pdf", p_both, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across cantons)
  
  weighted_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_nonborder_", year))  # use your LCC graphs
    
    # ensure cross_canton exists (1 = cross, 0 = within, NA if missing canton) 
    if (is.null(E(g)$cross_canton)) {
      ctn <- V(g)$ctn
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(ctn[ep[,1]]) | is.na(ctn[ep[,2]]) | ctn[ep[,1]] == "" | ctn[ep[,2]] == ""
        diff <- ctn[ep[,1]] != ctn[ep[,2]]
        E(g)$cross_canton <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_canton <- integer(0)
      }
      assign(paste0("g_firms_lcc_agglo_nonborder_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_agglo_nonborder_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_canton and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_canton
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_w_df <- bind_rows(weighted_summary)
  print(cross_summary_w_df)
  
  # Plot that
  share_w_long <- cross_summary_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-canton")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-canton (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo_nonborder/fig_within_vs_cross_weighted.pdf", p_both_w, width = 11.69, height = 11.69)
  
  ## Direct edge counts  (share within/across languages)
  
  # Storage for summary per year
  cross_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_nonborder_", year))   # your LCC graph for that year
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      # no edges or no nodes: set empty attribute + summary
      if (ecount(g) == 0) E(g)$cross_lang <- integer(0)
      cross_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_agglo_nonborder_", year), g)
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
    assign(paste0("g_firms_lcc_agglo_nonborder_", year), g)
    
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
    select(year, share_cross_lang) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_lang <- ggplot(share_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Share within language region")) +
    labs(x = "Year", y = "Share of edges",
         title = "Across language regions edges over time",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_lang)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo_nonborder/fig_within_vs_cross_lang.pdf", p_both_lang, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across lang)
  
  weighted_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_nonborder_", year))  # use your LCC graphs
    
    # ensure cross_lang exists (1 = cross, 0 = within, NA if missing lang) 
    if (is.null(E(g)$cross_lang)) {
      lang <- V(g)$language
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]] == "" | lang[ep[,2]] == ""
        diff <- lang[ep[,1]] != lang[ep[,2]]
        E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_lang <- integer(0)
      }
      assign(paste0("g_firms_lcc_agglo_nonborder_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_agglo_nonborder_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_lang and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_lang
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary_lang[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_lang_w_df <- bind_rows(weighted_summary_lang)
  print(cross_summary_lang_w_df)
  
  # Plot that
  share_w_long_lang <- cross_summary_lang_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-language")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-language (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo_nonborder/fig_within_vs_cross_lang_weighted.pdf", p_both_w, width = 11.69, height = 11.69)        
  
  ## Comparing within/across cantons & within/across language regions 
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
         title = "Cross edges — Canton (red) vs Language (black)") +
    theme_minimal(base_size = 13)
  
  print(p_unw)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo_nonborder/fig_within_cross_canton_language_unweighted.pdf", p_unw, width = 11.69, height = 11.69)
  
  # Weighted
  w_combined <- bind_rows(
    share_w_long %>%
      mutate(family = "canton"),
    share_w_long_lang %>%
      mutate(family = "language")
  ) %>%
    mutate(
      kind = ifelse(grepl("cross", measure), "cross", "within"),
      color_key = ifelse(family == "canton", "canton", "language"),
      lty_key   = ifelse(kind == "cross", "solid", "dashed")
    )
  
  p_w <- ggplot(w_combined, aes(x = year, y = value,
                                group = interaction(family, kind))) +
    geom_line(aes(color = color_key, linetype = lty_key), na.rm = TRUE) +
    geom_point(aes(color = color_key), na.rm = TRUE, size = 1.8) +
    scale_color_manual(values = c(canton = "red", language = "black"),
                       labels = c("Canton", "Language region"), name = NULL) +
    scale_linetype_identity(guide = "none") +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross edges — Canton (red) vs Language (black), weighted") +
    theme_minimal(base_size = 13)
  
  print(p_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo_nonborder/fig_within_cross_canton_language_weighted.pdf", p_w, width = 11.69, height = 11.69)  
  
################################################################################
  ### Border Agglo
  ##  Direct edge counts  (share within/across cantons)
  
  # Empty list to hold results
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    
    # Induce subgraph: only firms in agglomeration municipalities
    g_sub <- induced_subgraph(g, vids = which(V(g)$dummy_border_agglo == 1))
    
    # Save back with a new name
    assign(paste0("g_firms_lcc_agglo_border_", year), g_sub)
    
    # Optional: print quick stats
    cat("Year:", year,
        "- nodes:", vcount(g_sub),
        "- edges:", ecount(g_sub), "\n")
  }
  
  cross_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_border_", year))   # call graph
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
      cross_summary[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_agglo_border_", year), g)
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
    assign(paste0("g_firms_lcc_agglo_border_", year), g)
    
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
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo_border/fig_within_vs_cross.pdf", p_both, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across cantons)
  
  weighted_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_border_", year))  # use your LCC graphs
    
    # ensure cross_canton exists (1 = cross, 0 = within, NA if missing canton) 
    if (is.null(E(g)$cross_canton)) {
      ctn <- V(g)$ctn
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(ctn[ep[,1]]) | is.na(ctn[ep[,2]]) | ctn[ep[,1]] == "" | ctn[ep[,2]] == ""
        diff <- ctn[ep[,1]] != ctn[ep[,2]]
        E(g)$cross_canton <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_canton <- integer(0)
      }
      assign(paste0("g_firms_lcc_agglo_border_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_agglo_border_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_canton and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_canton
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_w_df <- bind_rows(weighted_summary)
  print(cross_summary_w_df)
  
  # Plot that
  share_w_long <- cross_summary_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-canton")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-canton (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo_border/fig_within_vs_cross_weighted.pdf", p_both_w, width = 11.69, height = 11.69)
  
  ## Direct edge counts  (share within/across languages)
  
  # Storage for summary per year
  cross_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_border_", year))   # your LCC graph for that year
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      # no edges or no nodes: set empty attribute + summary
      if (ecount(g) == 0) E(g)$cross_lang <- integer(0)
      cross_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_agglo_border_", year), g)
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
    assign(paste0("g_firms_lcc_agglo_border_", year), g)
    
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
    select(year, share_cross_lang) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_lang <- ggplot(share_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Share within language region")) +
    labs(x = "Year", y = "Share of edges",
         title = "Across language regions edges over time",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_lang)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo_border/fig_within_vs_cross_lang.pdf", p_both_lang, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across lang)
  
  weighted_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_border_", year))  # use your LCC graphs
    
    # ensure cross_lang exists (1 = cross, 0 = within, NA if missing lang) 
    if (is.null(E(g)$cross_lang)) {
      lang <- V(g)$language
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]] == "" | lang[ep[,2]] == ""
        diff <- lang[ep[,1]] != lang[ep[,2]]
        E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_lang <- integer(0)
      }
      assign(paste0("g_firms_lcc_agglo_border_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_agglo_border_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_lang and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_lang
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary_lang[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_lang_w_df <- bind_rows(weighted_summary_lang)
  print(cross_summary_lang_w_df)
  
  # Plot that
  share_w_long_lang <- cross_summary_lang_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-language")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-language (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo_border/fig_within_vs_cross_lang_weighted.pdf", p_both_w, width = 11.69, height = 11.69)        
  
  ## Comparing within/across cantons & within/across language regions 
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
         title = "Cross edges — Canton (red) vs Language (black)") +
    theme_minimal(base_size = 13)
  
  print(p_unw)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo_border/fig_within_cross_canton_language_unweighted.pdf", p_unw, width = 11.69, height = 11.69)
  
  # Weighted
  w_combined <- bind_rows(
    share_w_long %>%
      mutate(family = "canton"),
    share_w_long_lang %>%
      mutate(family = "language")
  ) %>%
    mutate(
      kind = ifelse(grepl("cross", measure), "cross", "within"),
      color_key = ifelse(family == "canton", "canton", "language"),
      lty_key   = ifelse(kind == "cross", "solid", "dashed")
    )
  
  p_w <- ggplot(w_combined, aes(x = year, y = value,
                                group = interaction(family, kind))) +
    geom_line(aes(color = color_key, linetype = lty_key), na.rm = TRUE) +
    geom_point(aes(color = color_key), na.rm = TRUE, size = 1.8) +
    scale_color_manual(values = c(canton = "red", language = "black"),
                       labels = c("Canton", "Language region"), name = NULL) +
    scale_linetype_identity(guide = "none") +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross edges — Canton (red) vs Language (black), weighted") +
    theme_minimal(base_size = 13)
  
  print(p_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo_border/fig_within_cross_canton_language_weighted.pdf", p_w, width = 11.69, height = 11.69)  
  
################################################################################
  ### Median distance from borders
  ##  Direct edge counts  (share within/across cantons)
  
  # Empty list to hold results
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    
    # Induce subgraph: only firms in agglomeration municipalities
    g_sub <- induced_subgraph(g, vids = which(V(g)$dummy_median == 1))
    
    # Save back with a new name
    assign(paste0("g_firms_lcc_median_", year), g_sub)
    
    # Optional: print quick stats
    cat("Year:", year,
        "- nodes:", vcount(g_sub),
        "- edges:", ecount(g_sub), "\n")
  }
  
  cross_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_median_", year))   # call graph
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
      cross_summary[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_median_", year), g)
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
    assign(paste0("g_firms_lcc_median_", year), g)
    
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
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median/fig_within_vs_cross.pdf", p_both, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across cantons)
  
  weighted_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_median_", year))  # use your LCC graphs
    
    # ensure cross_canton exists (1 = cross, 0 = within, NA if missing canton) 
    if (is.null(E(g)$cross_canton)) {
      ctn <- V(g)$ctn
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(ctn[ep[,1]]) | is.na(ctn[ep[,2]]) | ctn[ep[,1]] == "" | ctn[ep[,2]] == ""
        diff <- ctn[ep[,1]] != ctn[ep[,2]]
        E(g)$cross_canton <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_canton <- integer(0)
      }
      assign(paste0("g_firms_lcc_median_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_median_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_canton and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_canton
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_w_df <- bind_rows(weighted_summary)
  print(cross_summary_w_df)
  
  # Plot that
  share_w_long <- cross_summary_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-canton")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-canton (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median/fig_within_vs_cross_weighted.pdf", p_both_w, width = 11.69, height = 11.69)
  
  ## Direct edge counts  (share within/across languages)
  
  # Storage for summary per year
  cross_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_median_", year))   # your LCC graph for that year
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      # no edges or no nodes: set empty attribute + summary
      if (ecount(g) == 0) E(g)$cross_lang <- integer(0)
      cross_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_median_", year), g)
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
    assign(paste0("g_firms_lcc_median_", year), g)
    
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
    select(year, share_cross_lang) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_lang <- ggplot(share_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Share within language region")) +
    labs(x = "Year", y = "Share of edges",
         title = "Across language regions edges over time",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_lang)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median/fig_within_vs_cross_lang.pdf", p_both_lang, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across lang)
  
  weighted_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_median_", year))  # use your LCC graphs
    
    # ensure cross_lang exists (1 = cross, 0 = within, NA if missing lang) 
    if (is.null(E(g)$cross_lang)) {
      lang <- V(g)$language
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]] == "" | lang[ep[,2]] == ""
        diff <- lang[ep[,1]] != lang[ep[,2]]
        E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_lang <- integer(0)
      }
      assign(paste0("g_firms_lcc_median_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_median_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_lang and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_lang
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary_lang[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_lang_w_df <- bind_rows(weighted_summary_lang)
  print(cross_summary_lang_w_df)
  
  # Plot that
  share_w_long_lang <- cross_summary_lang_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-language")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-language (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median/fig_within_vs_cross_lang_weighted.pdf", p_both_w, width = 11.69, height = 11.69)        
  
  ## Comparing within/across cantons & within/across language regions 
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
         title = "Cross edges — Canton (red) vs Language (black)") +
    theme_minimal(base_size = 13)
  
  print(p_unw)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median/fig_within_cross_canton_language_unweighted.pdf", p_unw, width = 11.69, height = 11.69)
  
  # Weighted
  w_combined <- bind_rows(
    share_w_long %>%
      mutate(family = "canton"),
    share_w_long_lang %>%
      mutate(family = "language")
  ) %>%
    mutate(
      kind = ifelse(grepl("cross", measure), "cross", "within"),
      color_key = ifelse(family == "canton", "canton", "language"),
      lty_key   = ifelse(kind == "cross", "solid", "dashed")
    )
  
  p_w <- ggplot(w_combined, aes(x = year, y = value,
                                group = interaction(family, kind))) +
    geom_line(aes(color = color_key, linetype = lty_key), na.rm = TRUE) +
    geom_point(aes(color = color_key), na.rm = TRUE, size = 1.8) +
    scale_color_manual(values = c(canton = "red", language = "black"),
                       labels = c("Canton", "Language region"), name = NULL) +
    scale_linetype_identity(guide = "none") +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross edges — Canton (red) vs Language (black), weighted") +
    theme_minimal(base_size = 13)
  
  print(p_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median/fig_within_cross_canton_language_weighted.pdf", p_w, width = 11.69, height = 11.69)  
  
################################################################################
  ### Median distance from borders & not in agglo
  ##  Direct edge counts  (share within/across cantons)
  
  # Empty list to hold results
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    
    # Induce subgraph: only firms in agglomeration municipalities
    g_sub <- induced_subgraph(g, vids = which(V(g)$dummy_median_noagglo == 1))
    
    # Save back with a new name
    assign(paste0("g_firms_lcc_median_noagglo_", year), g_sub)
    
    # Optional: print quick stats
    cat("Year:", year,
        "- nodes:", vcount(g_sub),
        "- edges:", ecount(g_sub), "\n")
  }
  
  cross_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_median_noagglo_", year))   # call graph
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
      cross_summary[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_median_noagglo_", year), g)
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
    assign(paste0("g_firms_lcc_median_noagglo_", year), g)
    
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
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median_noagglo/fig_within_vs_cross.pdf", p_both, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across cantons)
  
  weighted_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_median_noagglo_", year))  # use your LCC graphs
    
    # ensure cross_canton exists (1 = cross, 0 = within, NA if missing canton) 
    if (is.null(E(g)$cross_canton)) {
      ctn <- V(g)$ctn
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(ctn[ep[,1]]) | is.na(ctn[ep[,2]]) | ctn[ep[,1]] == "" | ctn[ep[,2]] == ""
        diff <- ctn[ep[,1]] != ctn[ep[,2]]
        E(g)$cross_canton <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_canton <- integer(0)
      }
      assign(paste0("g_firms_lcc_median_noagglo_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_median_noagglo_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_canton and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_canton
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_w_df <- bind_rows(weighted_summary)
  print(cross_summary_w_df)
  
  # Plot that
  share_w_long <- cross_summary_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-canton")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-canton (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median_noagglo/fig_within_vs_cross_weighted.pdf", p_both_w, width = 11.69, height = 11.69)
  
  ## Direct edge counts  (share within/across languages)
  
  # Storage for summary per year
  cross_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_median_noagglo_", year))   # your LCC graph for that year
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      # no edges or no nodes: set empty attribute + summary
      if (ecount(g) == 0) E(g)$cross_lang <- integer(0)
      cross_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_median_noagglo_", year), g)
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
    assign(paste0("g_firms_lcc_median_noagglo_", year), g)
    
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
    select(year, share_cross_lang) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_lang <- ggplot(share_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Share within language region")) +
    labs(x = "Year", y = "Share of edges",
         title = "Across language regions edges over time",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_lang)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median_noagglo/fig_within_vs_cross_lang.pdf", p_both_lang, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across lang)
  
  weighted_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_median_noagglo_", year))  # use your LCC graphs
    
    # ensure cross_lang exists (1 = cross, 0 = within, NA if missing lang) 
    if (is.null(E(g)$cross_lang)) {
      lang <- V(g)$language
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]] == "" | lang[ep[,2]] == ""
        diff <- lang[ep[,1]] != lang[ep[,2]]
        E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_lang <- integer(0)
      }
      assign(paste0("g_firms_lcc_median_noagglo_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_median_noagglo_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_lang and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_lang
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary_lang[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_lang_w_df <- bind_rows(weighted_summary_lang)
  print(cross_summary_lang_w_df)
  
  # Plot that
  share_w_long_lang <- cross_summary_lang_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-language")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-language (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median_noagglo/fig_within_vs_cross_lang_weighted.pdf", p_both_w, width = 11.69, height = 11.69)        
  
  ## Comparing within/across cantons & within/across language regions 
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
         title = "Cross edges — Canton (red) vs Language (black)") +
    theme_minimal(base_size = 13)
  
  print(p_unw)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median_noagglo/fig_within_cross_canton_language_unweighted.pdf", p_unw, width = 11.69, height = 11.69)
  
  # Weighted
  w_combined <- bind_rows(
    share_w_long %>%
      mutate(family = "canton"),
    share_w_long_lang %>%
      mutate(family = "language")
  ) %>%
    mutate(
      kind = ifelse(grepl("cross", measure), "cross", "within"),
      color_key = ifelse(family == "canton", "canton", "language"),
      lty_key   = ifelse(kind == "cross", "solid", "dashed")
    )
  
  p_w <- ggplot(w_combined, aes(x = year, y = value,
                                group = interaction(family, kind))) +
    geom_line(aes(color = color_key, linetype = lty_key), na.rm = TRUE) +
    geom_point(aes(color = color_key), na.rm = TRUE, size = 1.8) +
    scale_color_manual(values = c(canton = "red", language = "black"),
                       labels = c("Canton", "Language region"), name = NULL) +
    scale_linetype_identity(guide = "none") +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross edges — Canton (red) vs Language (black), weighted") +
    theme_minimal(base_size = 13)
  
  print(p_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median_noagglo/fig_within_cross_canton_language_weighted.pdf", p_w, width = 11.69, height = 11.69)  

################################################################################
  ### Median distance from borders & in agglo
  ##  Direct edge counts  (share within/across cantons)
  
  # Empty list to hold results
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    
    # Induce subgraph: only firms in agglomeration municipalities
    g_sub <- induced_subgraph(g, vids = which(V(g)$dummy_median_agglo == 1))
    
    # Save back with a new name
    assign(paste0("g_firms_lcc_median_agglo_", year), g_sub)
    
    # Optional: print quick stats
    cat("Year:", year,
        "- nodes:", vcount(g_sub),
        "- edges:", ecount(g_sub), "\n")
  }
  
  cross_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_median_agglo_", year))   # call graph
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
      cross_summary[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_median_agglo_", year), g)
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
    assign(paste0("g_firms_lcc_median_agglo_", year), g)
    
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
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median_agglo/fig_within_vs_cross.pdf", p_both, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across cantons)
  
  weighted_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_median_agglo_", year))  # use your LCC graphs
    
    # ensure cross_canton exists (1 = cross, 0 = within, NA if missing canton) 
    if (is.null(E(g)$cross_canton)) {
      ctn <- V(g)$ctn
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(ctn[ep[,1]]) | is.na(ctn[ep[,2]]) | ctn[ep[,1]] == "" | ctn[ep[,2]] == ""
        diff <- ctn[ep[,1]] != ctn[ep[,2]]
        E(g)$cross_canton <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_canton <- integer(0)
      }
      assign(paste0("g_firms_lcc_median_agglo_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_median_agglo_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_canton and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_canton
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_w_df <- bind_rows(weighted_summary)
  print(cross_summary_w_df)
  
  # Plot that
  share_w_long <- cross_summary_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-canton")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-canton (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median_agglo/fig_within_vs_cross_weighted.pdf", p_both_w, width = 11.69, height = 11.69)
  
  ## Direct edge counts  (share within/across languages)
  
  # Storage for summary per year
  cross_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_median_agglo_", year))   
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      # no edges or no nodes: set empty attribute + summary
      if (ecount(g) == 0) E(g)$cross_lang <- integer(0)
      cross_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_median_agglo_", year), g)
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
    assign(paste0("g_firms_lcc_median_agglo_", year), g)
    
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
    select(year, share_cross_lang) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_lang <- ggplot(share_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Share within language region")) +
    labs(x = "Year", y = "Share of edges",
         title = "Across language regions edges over time",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_lang)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median_agglo/fig_within_vs_cross_lang.pdf", p_both_lang, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across lang)
  
  weighted_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_median_agglo_", year))  # use your LCC graphs
    
    # ensure cross_lang exists (1 = cross, 0 = within, NA if missing lang) 
    if (is.null(E(g)$cross_lang)) {
      lang <- V(g)$language
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]] == "" | lang[ep[,2]] == ""
        diff <- lang[ep[,1]] != lang[ep[,2]]
        E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_lang <- integer(0)
      }
      assign(paste0("g_firms_lcc_median_agglo_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_median_agglo_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_lang and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_lang
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary_lang[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_lang_w_df <- bind_rows(weighted_summary_lang)
  print(cross_summary_lang_w_df)
  
  # Plot that
  share_w_long_lang <- cross_summary_lang_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-language")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-language (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median_agglo/fig_within_vs_cross_lang_weighted.pdf", p_both_w, width = 11.69, height = 11.69)        
  
  ## Comparing within/across cantons & within/across language regions 
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
         title = "Cross edges — Canton (red) vs Language (black)") +
    theme_minimal(base_size = 13)
  
  print(p_unw)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median_agglo/fig_within_cross_canton_language_unweighted.pdf", p_unw, width = 11.69, height = 11.69)
  
  # Weighted
  w_combined <- bind_rows(
    share_w_long %>%
      mutate(family = "canton"),
    share_w_long_lang %>%
      mutate(family = "language")
  ) %>%
    mutate(
      kind = ifelse(grepl("cross", measure), "cross", "within"),
      color_key = ifelse(family == "canton", "canton", "language"),
      lty_key   = ifelse(kind == "cross", "solid", "dashed")
    )
  
  p_w <- ggplot(w_combined, aes(x = year, y = value,
                                group = interaction(family, kind))) +
    geom_line(aes(color = color_key, linetype = lty_key), na.rm = TRUE) +
    geom_point(aes(color = color_key), na.rm = TRUE, size = 1.8) +
    scale_color_manual(values = c(canton = "red", language = "black"),
                       labels = c("Canton", "Language region"), name = NULL) +
    scale_linetype_identity(guide = "none") +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross edges — Canton (red) vs Language (black), weighted") +
    theme_minimal(base_size = 13)
  
  print(p_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median_agglo/fig_within_cross_canton_language_weighted.pdf", p_w, width = 11.69, height = 11.69)  

################################################################################
  ### Bordering
  ##  Direct edge counts  (share within/across cantons)
  
  # Empty list to hold results
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    
    # Induce subgraph: only firms in agglomeration municipalities
    g_sub <- induced_subgraph(g, vids = which(V(g)$dummy_touch_border == 1))
    
    # Save back with a new name
    assign(paste0("g_firms_lcc_border_", year), g_sub)
    
    # Optional: print quick stats
    cat("Year:", year,
        "- nodes:", vcount(g_sub),
        "- edges:", ecount(g_sub), "\n")
  }
  
  cross_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_border_", year))   # call graph
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
      cross_summary[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_border_", year), g)
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
    assign(paste0("g_firms_lcc_border_", year), g)
    
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
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/border/fig_within_vs_cross.pdf", p_both, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across cantons)
  
  weighted_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_border_", year))  
    
    # ensure cross_canton exists (1 = cross, 0 = within, NA if missing canton) 
    if (is.null(E(g)$cross_canton)) {
      ctn <- V(g)$ctn
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(ctn[ep[,1]]) | is.na(ctn[ep[,2]]) | ctn[ep[,1]] == "" | ctn[ep[,2]] == ""
        diff <- ctn[ep[,1]] != ctn[ep[,2]]
        E(g)$cross_canton <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_canton <- integer(0)
      }
      assign(paste0("g_firms_lcc_border_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_border_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_canton and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_canton
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_w_df <- bind_rows(weighted_summary)
  print(cross_summary_w_df)
  
  # Plot that
  share_w_long <- cross_summary_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-canton")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-canton (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/border/fig_within_vs_cross_weighted.pdf", p_both_w, width = 11.69, height = 11.69)
  
  ## Direct edge counts  (share within/across languages)
  
  # Storage for summary per year
  cross_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_border_", year))   # your LCC graph for that year
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      # no edges or no nodes: set empty attribute + summary
      if (ecount(g) == 0) E(g)$cross_lang <- integer(0)
      cross_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_border_", year), g)
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
    assign(paste0("g_firms_lcc_border_", year), g)
    
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
    select(year, share_cross_lang) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_lang <- ggplot(share_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Share within language region")) +
    labs(x = "Year", y = "Share of edges",
         title = "Across language regions edges over time",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_lang)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/border/fig_within_vs_cross_lang.pdf", p_both_lang, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across lang)
  
  weighted_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_border_", year))  # use your LCC graphs
    
    # ensure cross_lang exists (1 = cross, 0 = within, NA if missing lang) 
    if (is.null(E(g)$cross_lang)) {
      lang <- V(g)$language
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]] == "" | lang[ep[,2]] == ""
        diff <- lang[ep[,1]] != lang[ep[,2]]
        E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_lang <- integer(0)
      }
      assign(paste0("g_firms_lcc_border_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_border_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_lang and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_lang
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary_lang[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_lang_w_df <- bind_rows(weighted_summary_lang)
  print(cross_summary_lang_w_df)
  
  # Plot that
  share_w_long_lang <- cross_summary_lang_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-language")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-language (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/border/fig_within_vs_cross_lang_weighted.pdf", p_both_w, width = 11.69, height = 11.69)        
  
  ## Comparing within/across cantons & within/across language regions 
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
         title = "Cross edges — Canton (red) vs Language (black)") +
    theme_minimal(base_size = 13)
  
  print(p_unw)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/border/fig_within_cross_canton_language_unweighted.pdf", p_unw, width = 11.69, height = 11.69)
  
  # Weighted
  w_combined <- bind_rows(
    share_w_long %>%
      mutate(family = "canton"),
    share_w_long_lang %>%
      mutate(family = "language")
  ) %>%
    mutate(
      kind = ifelse(grepl("cross", measure), "cross", "within"),
      color_key = ifelse(family == "canton", "canton", "language"),
      lty_key   = ifelse(kind == "cross", "solid", "dashed")
    )
  
  p_w <- ggplot(w_combined, aes(x = year, y = value,
                                group = interaction(family, kind))) +
    geom_line(aes(color = color_key, linetype = lty_key), na.rm = TRUE) +
    geom_point(aes(color = color_key), na.rm = TRUE, size = 1.8) +
    scale_color_manual(values = c(canton = "red", language = "black"),
                       labels = c("Canton", "Language region"), name = NULL) +
    scale_linetype_identity(guide = "none") +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross edges — Canton (red) vs Language (black), weighted") +
    theme_minimal(base_size = 13)
  
  print(p_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/border/fig_within_cross_canton_language_weighted.pdf", p_w, width = 11.69, height = 11.69)  

  ################################################################################
  ### Bordering & non agglo
  ##  Direct edge counts  (share within/across cantons)
  
  # Empty list to hold results
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    
    # Induce subgraph: only firms in agglomeration municipalities
    g_sub <- induced_subgraph(g, vids = which(V(g)$dummy_touch_border_noagglo == 1))
    
    # Save back with a new name
    assign(paste0("g_firms_lcc_border_noagglo_", year), g_sub)
    
    # Optional: print quick stats
    cat("Year:", year,
        "- nodes:", vcount(g_sub),
        "- edges:", ecount(g_sub), "\n")
  }
  
  cross_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_border_noagglo_", year))   # call graph
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      if (ecount(g) == 0) E(g)$cross_canton <- integer(0)
      cross_summary[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_border_noagglo_", year), g)
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
    assign(paste0("g_firms_lcc_border_noagglo_", year), g)
    
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
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/border_noagglo/fig_within_vs_cross.pdf", p_both, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across cantons)
  
  weighted_summary <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_border_noagglo_", year))  # use your LCC graphs
    
    # ensure cross_canton exists (1 = cross, 0 = within, NA if missing canton) 
    if (is.null(E(g)$cross_canton)) {
      ctn <- V(g)$ctn
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(ctn[ep[,1]]) | is.na(ctn[ep[,2]]) | ctn[ep[,1]] == "" | ctn[ep[,2]] == ""
        diff <- ctn[ep[,1]] != ctn[ep[,2]]
        E(g)$cross_canton <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_canton <- integer(0)
      }
      assign(paste0("g_firms_lcc_border_noagglo_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_border_noagglo_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_canton and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_canton
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_w_df <- bind_rows(weighted_summary)
  print(cross_summary_w_df)
  
  # Plot that
  share_w_long <- cross_summary_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-canton")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-canton (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/border_noagglo/fig_within_vs_cross_weighted.pdf", p_both_w, width = 11.69, height = 11.69)
  
  ## Direct edge counts  (share within/across languages)
  
  # Storage for summary per year
  cross_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_border_noagglo_", year))   # your LCC graph for that year
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      # no edges or no nodes: set empty attribute + summary
      if (ecount(g) == 0) E(g)$cross_lang <- integer(0)
      cross_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_border_noagglo_", year), g)
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
    assign(paste0("g_firms_lcc_border_noagglo_", year), g)
    
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
    select(year, share_cross_lang) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_lang <- ggplot(share_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Share within language region")) +
    labs(x = "Year", y = "Share of edges",
         title = "Across language regions edges over time",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_lang)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/border_noagglo/fig_within_vs_cross_lang.pdf", p_both_lang, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across lang)
  
  weighted_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_border_noagglo_", year))  # use your LCC graphs
    
    # ensure cross_lang exists (1 = cross, 0 = within, NA if missing lang) 
    if (is.null(E(g)$cross_lang)) {
      lang <- V(g)$language
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]] == "" | lang[ep[,2]] == ""
        diff <- lang[ep[,1]] != lang[ep[,2]]
        E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_lang <- integer(0)
      }
      assign(paste0("g_firms_lcc_border_noagglo_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_border_noagglo_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_lang and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_lang
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary_lang[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_lang_w_df <- bind_rows(weighted_summary_lang)
  print(cross_summary_lang_w_df)
  
  # Plot that
  share_w_long_lang <- cross_summary_lang_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-language")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-language (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/border_noagglo/fig_within_vs_cross_lang_weighted.pdf", p_both_w, width = 11.69, height = 11.69)        
  
  ## Comparing within/across cantons & within/across language regions 
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
         title = "Cross edges — Canton (red) vs Language (black)") +
    theme_minimal(base_size = 13)
  
  print(p_unw)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/border_noagglo/fig_within_cross_canton_language_unweighted.pdf", p_unw, width = 11.69, height = 11.69)
  
  # Weighted
  w_combined <- bind_rows(
    share_w_long %>%
      mutate(family = "canton"),
    share_w_long_lang %>%
      mutate(family = "language")
  ) %>%
    mutate(
      kind = ifelse(grepl("cross", measure), "cross", "within"),
      color_key = ifelse(family == "canton", "canton", "language"),
      lty_key   = ifelse(kind == "cross", "solid", "dashed")
    )
  
  p_w <- ggplot(w_combined, aes(x = year, y = value,
                                group = interaction(family, kind))) +
    geom_line(aes(color = color_key, linetype = lty_key), na.rm = TRUE) +
    geom_point(aes(color = color_key), na.rm = TRUE, size = 1.8) +
    scale_color_manual(values = c(canton = "red", language = "black"),
                       labels = c("Canton", "Language region"), name = NULL) +
    scale_linetype_identity(guide = "none") +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross edges — Canton (red) vs Language (black), weighted") +
    theme_minimal(base_size = 13)
  
  print(p_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/border_noagglo/fig_within_cross_canton_language_weighted.pdf", p_w, width = 11.69, height = 11.69)  
  
################################################################################
### Language borders and within/across borders link
  
  ### Create subgraph with close to border nodes only 
  
  # Empty list to hold results
  for (year in years) {
    g <- get(paste0("g_firms_lcc_", year))
    
    # Induce subgraph: only firms in agglomeration municipalities
    g_sub <- induced_subgraph(g, vids = which(V(g)$dummy_10km_rosti == 1))
    
    # Save back with a new name
    assign(paste0("g_firms_lcc_10km_rosti_", year), g_sub)
    
    # Optional: print quick stats
    cat("Year:", year,
        "- nodes:", vcount(g_sub),
        "- edges:", ecount(g_sub), "\n")
  }
  
  ## Direct edge counts  (share within/across languages)
  
  # Storage for summary per year
  cross_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_10km_rosti_", year))   # your LCC graph for that year
    
    if (ecount(g) == 0 || vcount(g) == 0) {
      # no edges or no nodes: set empty attribute + summary
      if (ecount(g) == 0) E(g)$cross_lang <- integer(0)
      cross_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within = 0L, cross = 0L, share_within = NA_real_, share_cross = NA_real_
      )
      assign(paste0("g_firms_lcc_10km_rosti_", year), g)
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
    assign(paste0("g_firms_lcc_10km_rosti_", year), g)
    
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
    select(year, share_cross_lang) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_lang <- ggplot(share_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Share across language region")) +
    labs(x = "Year", y = "Share of edges",
         title = "Across language regions edges over time",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_lang)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/rosti/fig_within_vs_cross_lang.pdf", p_both_lang, width = 11.69, height = 11.69)
  
  ## Direct weighted edge counts (within/across lang)
  
  weighted_summary_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_10km_rosti_", year))  # use your LCC graphs
    
    # ensure cross_lang exists (1 = cross, 0 = within, NA if missing lang) 
    if (is.null(E(g)$cross_lang)) {
      lang <- V(g)$language
      if (ecount(g) > 0 && vcount(g) > 0) {
        ep   <- ends(g, E(g), names = FALSE)
        bad  <- is.na(lang[ep[,1]]) | is.na(lang[ep[,2]]) | lang[ep[,1]] == "" | lang[ep[,2]] == ""
        diff <- lang[ep[,1]] != lang[ep[,2]]
        E(g)$cross_lang <- ifelse(bad, NA_integer_, as.integer(diff))
      } else {
        E(g)$cross_lang <- integer(0)
      }
      assign(paste0("g_firms_lcc_10km_rosti_", year), g)
    }
    
    # Weights
    w <- NULL
    if (!is.null(E(g)$weight)) {
      w <- E(g)$weight
    } else if (!is.null(E(g)$n_shared_directors)) {
      w <- E(g)$n_shared_directors
      E(g)$weight <- w  # store for later reuse
      assign(paste0("g_firms_lcc_10km_rosti_", year), g)
    } else {
      # no weights defined; treat as 1
      w <- rep(1, ecount(g))
    }
    
    # keep only edges with defined cross_lang and non-missing weights
    if (ecount(g) == 0) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    cc <- E(g)$cross_lang
    keep <- !is.na(cc) & !is.na(w)
    if (!any(keep)) {
      weighted_summary_lang[[as.character(year)]] <- data.frame(
        year = year, within_w = 0, cross_w = 0,
        share_within_w = NA_real_, share_cross_w = NA_real_
      )
      next
    }
    
    w_keep <- w[keep]
    cc_keep <- cc[keep]
    
    within_w <- sum(w_keep[cc_keep == 0], na.rm = TRUE)
    cross_w  <- sum(w_keep[cc_keep == 1], na.rm = TRUE)
    tot_w    <- within_w + cross_w
    
    weighted_summary_lang[[as.character(year)]] <- data.frame(
      year = year,
      within_w = within_w,
      cross_w  = cross_w,
      share_within_w = if (tot_w > 0) within_w / tot_w else NA_real_,
      share_cross_w  = if (tot_w > 0) cross_w  / tot_w else NA_real_
    )
  }
  
  cross_summary_lang_w_df <- bind_rows(weighted_summary_lang)
  print(cross_summary_lang_w_df)
  
  # Plot that
  share_w_long_lang <- cross_summary_lang_w_df %>%
    select(year, share_cross_w) %>%
    pivot_longer(-year, names_to = "measure", values_to = "value")
  
  p_both_w <- ggplot(share_w_long_lang, aes(x = year, y = value, linetype = measure)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_linetype_manual(values = c("solid"),
                          labels = c("Weighted share cross-language")) +
    labs(x = "Year", y = "Weighted share of edges",
         title = "Cross-language (weighted)",
         linetype = NULL) +
    theme_minimal(base_size = 14)
  
  print(p_both_w)
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/rosti/fig_within_vs_cross_lang_weighted.pdf", p_both_w, width = 11.69, height = 11.69)        

################################################################################            
### E-I (External-Internal) Index in subsamples   
  # Median
  
  ei_results <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_median_", year))
    
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
    g <- get(paste0("g_firms_lcc_median_", year))
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
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/median/fig_ei_unw_vs_w.pdf", p_ei, width = 11.69, height = 11.69)
  
  # Notes
  #EI=−1: all ties are internal (only within cantons).
  #EI=+1: all ties are external (only across cantons).
  #EI=0: balance of internal and external ties (in relative terms).  
  
  ###############
  ### 10km rostigraben
  ### E-I (External-Internal) Index for Language   
  
  ei_results_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_10km_rosti_", year))
    
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
    g <- get(paste0("g_firms_lcc_10km_rosti_", year))
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
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/rosti/fig_ei_unw_vs_w.pdf", p_ei, width = 11.69, height = 11.69)
  
  ###
  # Agglo
  ei_results <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_", year))
    
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
    g <- get(paste0("g_firms_lcc_agglo_", year))
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
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo/fig_ei_unw_vs_w.pdf", p_ei, width = 11.69, height = 11.69)
  
  # Notes
  #EI=−1: all ties are internal (only within cantons).
  #EI=+1: all ties are external (only across cantons).
  #EI=0: balance of internal and external ties (in relative terms)
  
  ### E-I (External-Internal) Index for Language   
  
  ei_results_lang <- list()
  
  for (year in years) {
    g <- get(paste0("g_firms_lcc_agglo_", year))
    
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
    g <- get(paste0("g_firms_lcc_agglo_", year))
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
  ggsave("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/agglo/fig_ei_unw_vs_w.pdf", p_ei, width = 11.69, height = 11.69)
  