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

#Export
write.csv(all_full_ranks_df,
          "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/31_Place_Based_Policies/data/derived/all_full_ranks_df.csv",
          row.names = FALSE)

write_dta(all_full_ranks_df, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/31_Place_Based_Policies/data/derived/all_full_ranks_df.dta")

################################################################################
################################################################################
### Dir-to-dir : node-level measures

# Create storage lists
all_full_ranks <- list()
all_top25 <- list()
all_avg_rank <- list()

# Output directory
outdir <- "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/noSI/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Centrality measures - Loop over all years    
for (year in years) {
  giant_g_dir <- get(paste0("g_directors_lcc_", year))
  
  # Centralities
  V(giant_g_dir)$degree      <- degree(giant_g_dir)
  V(giant_g_dir)$betweenness <- betweenness(giant_g_dir, normalized = TRUE)
  V(giant_g_dir)$closeness   <- closeness(giant_g_dir, normalized = TRUE)
  V(giant_g_dir)$eigenvector <- eigen_centrality(giant_g_dir)$vector
  V(giant_g_dir)$pagerank    <- page_rank(giant_g_dir)$vector
  
  centrality_df <- data.frame(
    firm_id     = V(giant_g_dir)$name,
    firm_name   = V(giant_g_dir)$cname,
    degree      = V(giant_g_dir)$degree,
    betweenness = V(giant_g_dir)$betweenness,
    closeness   = V(giant_g_dir)$closeness,
    eigenvector = V(giant_g_dir)$eigenvector,
    pagerank    = V(giant_g_dir)$pagerank,
    stringsAsFactors = FALSE
  )
  
  # Rank all firms by a measure
  make_full_rank_dir <- function(var) {
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

#Export
write.csv(all_full_ranks_df,
          "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/31_Place_Based_Policies/data/derived/all_full_ranks_df.csv",
          row.names = FALSE)

write_dta(all_full_ranks_df, "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/31_Place_Based_Policies/data/derived/all_full_ranks_df.dta")

