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
    
    for (year in years) {
      g_bi <- get(paste0("g_bi_", year))
      
      # Check directors (type == TRUE)
      missing_dir <- sum(is.na(V(g_bi)$ctn[V(g_bi)$type == TRUE]))
      total_dir   <- sum(V(g_bi)$type == TRUE)
      
      # Check firms (type == FALSE)
      missing_firm <- sum(is.na(V(g_bi)$ctn[V(g_bi)$type == FALSE]))
      total_firm   <- sum(V(g_bi)$type == FALSE)
      
      cat("Year:", year, "\n")
      cat("  Directors:", total_dir, "nodes,", missing_dir, "missing ctn\n")
      cat("  Firms:    ", total_firm, "nodes,", missing_firm, "missing ctn\n\n")
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
  
################################################################################
################################################################################
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
  