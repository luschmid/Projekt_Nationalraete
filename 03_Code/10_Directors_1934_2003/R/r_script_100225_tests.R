install.packages(c("kableExtra", "knitr", "openxlsx", "devtools", "igraph", "bipartite", "asnipe", "assortnet", "ggplot2", "ggmap", "rnetcarto", "ecodist", "igraphdata", "statnet", "RColorBrewer", "tidyverse", "xlsx"))
install.packages("readstata13")
install.packages(c("tidygraph", "ggraph"))
install.packages("rnaturalearth")

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

# load map data from Fretz Parchet Robert-Nicoud
load("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/20_Highway_Firm/data/raw/maps/maps.RData")

#Path
setwd("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/14_Network_OneMode")

################################################################################
# Creation of bipartite network, firm-to-firm, director-to-director, mun-to-mun : all years

  base_path <- "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/"

#Store all edgelists for all years
  years <- c(1934, 1943, 1960, 1962, 1963, 1964, 1965, 1966, 1969,
             1972, 1975, 1979, 1980, 1981, 1982, 1983, 1984, 1985,
             1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994,
             1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003)
  
  for (year in years) {
    # 1. Load edgelist
    path <- paste0(base_path, "edgeliste_", year, "_bipartite.dta")
    edgelist <- read_stata(path)
    assign(paste0("edgeliste_", year, "_bi"), edgelist)
    
    # 2. Create bipartite graph
    g_bi <- graph_from_data_frame(edgelist, directed = FALSE)
    V(g_bi)$type <- bipartite_mapping(g_bi)$type
    assign(paste0("g_bi_", year), g_bi)
    
    # 3. Create firm–firm and director–director projections
    projections <- bipartite_projection(g_bi)
    
    g_firms <- projections$proj1
    E(g_firms)$weight <- count_multiple(g_firms)
    assign(paste0("g_firms_", year), g_firms)
    
    g_directors <- projections$proj2
    E(g_directors)$weight <- count_multiple(g_directors)
    assign(paste0("g_directors_", year), g_directors)
  }
  
# Add nodes' attributes

# Firms newtork
  base_attr_path <- "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/Node_attribute_comp/"
  
  for (year in years) {
    # 1. Read the firm attribute file for the year
    attr_path <- paste0(base_attr_path, "comp_", year, "_attribute_complete.dta")
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

    # 5. Optionally reassign updated graph back into environment
    assign(graph_name, g_firms)
  }

#Directors network
  base_attr_path_dir <- "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/Node_attribute_pers/"
  
  for (year in years) {
    # 1. Read the director attribute file
    attr_path <- paste0(base_attr_path_dir, "pers_", year, "_attribute_complete.dta")
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

#Bipartite
  for (year in years) {
    graph_name <- paste0("g_bi_", year)
    g_bi <- get(graph_name)
    
    # -----------------------------------
    # Director attributes (type == TRUE)
    # -----------------------------------
    dir_attrs <- read_stata(paste0(base_attr_path_dir, "pers_", year, "_attribute_complete.dta"))
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
    firm_attrs <- read_stata(paste0(base_attr_path, "comp_", year, "_attribute_complete.dta"))
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
  
#Municipality Coordinates
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
  ### Macro-level description of network over time
  ################################################################################    
  
  ### Firm to firm
  ## Full Firm-to-Firm Network
  # Create empty list to store results
  firm_network_summaries <- list()
  
  # Loop over years
  for (year in years) {
    # Get the firm–firm graph object
    g_firm <- get(paste0("g_firms_", year))
    
    # Compute macro-level stats
    net_sum <- list(
      year = year,
      number_of_edges = gsize(g_firm),
      number_of_nodes = gorder(g_firm),
      density = edge_density(g_firm),
      is_connected = is_connected(g_firm),
      number_of_components = components(g_firm)$no,
      #diameter = if (is_connected(g_firm)) diameter(g_firm) else NA,
      #average_path_length = if (is_connected(g_firm)) mean_distance(g_firm) else NA,
      transitivity_global = transitivity(g_firm, type = "global"),
      assortativity_degree = assortativity_degree(g_firm)
    )
    
    # Store
    firm_network_summaries[[as.character(year)]] <- net_sum
  }
  
  # Convert to data.frame for summary or plotting
  library(dplyr)
  firm_network_summary_df <- bind_rows(firm_network_summaries)
  
  # View
  print(firm_network_summary_df) 
  
  # Cleaning of table
  firm_network_summary_df_clean <- firm_network_summary_df %>%
    mutate(
      year = as.integer(year),
      density = round(density, 3),
      transitivity_global = round(transitivity_global, 3),
      assortativity_degree = round(assortativity_degree, 3),
      number_of_components = as.integer(number_of_components),
      number_of_edges = as.integer(number_of_edges),
      number_of_nodes = as.integer(number_of_nodes)
    )
  
  # Convert to LaTeX
  latex_table <- xtable(firm_network_summary_df_clean,
                        caption = "Macro-level measures of the firm network",
                        label = "tab:firm_net_summary")
  
  # Print to .tex file
  print(latex_table,
        include.rownames = FALSE,  # cleaner
        file = "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/firm_net_summary.tex",
        booktabs = TRUE)           # nicer style
  
  ### In Largest Connected Component (LCC)
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
  
  # List to store results
    lcc_summary_list <- list()
  
    for (year in years) {
      # Get the LCC graph object from the environment
      g_lcc <- get(paste0("g_firms_lcc_", year))
      
      # Compute macro-level measures
      summary_stats <- list(
        year = year,
        number_of_edges = gsize(g_lcc),
        number_of_nodes = gorder(g_lcc),
        density = edge_density(g_lcc),
        is_connected = is_connected(g_lcc),
        number_of_components = components(g_lcc)$no,
        transitivity_global = transitivity(g_lcc, type = "global"),
        assortativity_degree = assortativity_degree(g_lcc)
      )
      
      # Store
      lcc_summary_list[[as.character(year)]] <- summary_stats
      }
  
    # Convert to data frame
    lcc_summary_df_firm <- bind_rows(lcc_summary_list)
  
    # View or export
    print(lcc_summary_df_firm)  
  
    # Cleaning of table
    lcc_summary_df_firm_clean <- lcc_summary_df_firm %>%
      mutate(
        year = as.integer(year),
        density = round(density, 3),
        transitivity_global = round(transitivity_global, 3),
        assortativity_degree = round(assortativity_degree, 3),
        number_of_edges = as.integer(number_of_edges),
        number_of_nodes = as.integer(number_of_nodes)
      )
    
    # Convert to LaTeX
    latex_table <- xtable(lcc_summary_df_firm_clean,
                          caption = "Macro-level measures of the firm network (LCC)",
                          label = "tab:lcc_summary")
    
    # Print to .tex file
    print(latex_table,
          include.rownames = FALSE,  # cleaner
          file = "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/lcc_summary_table.tex",
          booktabs = TRUE)           # nicer style

    ################################################################################    
    ### Micro-level description of network over time : centrality
    ################################################################################   
    # 7.4 Centrality at node-level 1943, 1960, 1972, 1980, 1990, 2000, 2003
    
    ##1943
    # 1. Identify components
    comp_1943 <- components(g_firms_1943)
    
    # 2. Extract the largest connected component (LCC)
    giant_g_1943 <- induced_subgraph(g_firms_1943, which(comp_1943$membership == which.max(comp_1943$csize)))
    
    # 3. Compute centrality measures on the LCC
    V(giant_g_1943)$degree      <- degree(giant_g_1943)
    V(giant_g_1943)$betweenness <- betweenness(giant_g_1943, normalized = TRUE)
    V(giant_g_1943)$closeness   <- closeness(giant_g_1943, normalized = TRUE)
    V(giant_g_1943)$eigenvector <- eigen_centrality(giant_g_1943)$vector
    V(giant_g_1943)$pagerank    <- page_rank(giant_g_1943)$vector
    
    # 4. Into a dataframe
    centrality_df_1943 <- data.frame(
      firm_id = V(giant_g_1943)$name,
      firm_name = V(giant_g_1943)$cname,
      degree = V(giant_g_1943)$degree,
      betweenness = V(giant_g_1943)$betweenness,
      closeness = V(giant_g_1943)$closeness,
      eigenvector = V(giant_g_1943)$eigenvector,
      pagerank = V(giant_g_1943)$pagerank,
      stringsAsFactors = FALSE
    )
    
    # 5. Get top 25 nodes for each measure
    top25_degree_1943 <- centrality_df_1943 %>% arrange(desc(degree)) %>% slice_head(n = 25)
    top25_betweenness_1943 <- centrality_df_1943 %>% arrange(desc(betweenness)) %>% slice_head(n = 25)
    top25_closeness_1943 <- centrality_df_1943 %>% arrange(desc(closeness)) %>% slice_head(n = 25)
    top25_eigenvector_1943 <- centrality_df_1943 %>% arrange(desc(eigenvector)) %>% slice_head(n = 25)
    top25_pagerank_1943 <- centrality_df_1943 %>% arrange(desc(pagerank)) %>% slice_head(n = 25)
    top25_all_1943 <- bind_rows(top25_degree_1943, top25_betweenness_1943, top25_closeness_1943, top25_eigenvector_1943, top25_pagerank_1943)
    
    # Add rank and centrality type to each
    top_degree_df_1943 <- top25_degree_1943 %>% mutate(rank = row_number(), centrality = "degree")
    top_bet_df_1943 <- top25_betweenness_1943 %>% mutate(rank = row_number(), centrality = "betweenness")
    top_close_df_1943 <- top25_closeness_1943 %>% mutate(rank = row_number(), centrality = "closeness")
    top_eigen_df_1943 <- top25_eigenvector_1943 %>% mutate(rank = row_number(), centrality = "eigenvector")
    top_pr_df_1943 <- top25_pagerank_1943 %>% mutate(rank = row_number(), centrality = "pagerank")
    
    # Combine all into one table
    top25_all_1943 <- bind_rows(top_degree_df_1943, top_bet_df_1943, top_close_df_1943, top_eigen_df_1943, top_pr_df_1943)
    
    #Wide table
    top25_wide_1943 <- list(
      degree = top25_degree_1943$firm_name,
      betweenness = top25_betweenness_1943$firm_name,
      closeness = top25_closeness_1943$firm_name,
      eigenvector = top25_eigenvector_1943$firm_name,
      pagerank = top25_pagerank_1943$firm_name
    ) %>% as.data.frame()
    
    # View the result
    print(top25_all_1943)
    print(top25_wide_1943)
    
    library(xtable)
    
    print(xtable(top25_all_1943,
                 caption = "Top 25 firms by centrality measures, 1943",
                 label = "tab:top25_all_1943"),
          include.rownames = FALSE,
          file = "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/top25_all_1943.tex",
          booktabs = TRUE)
    
    print(xtable(top25_wide_1943,
                 caption = "Top 25 firms (wide format), 1943",
                 label = "tab:top25_wide_1943"),
          include.rownames = FALSE,
          file = "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/top25_wide_1943.tex",
          booktabs = TRUE)
    
    ##1960
    # 1. Identify components
    comp_1960 <- components(g_firms_1960)
    
    # 2. Extract the largest connected component (LCC)
    giant_g_1960 <- induced_subgraph(g_firms_1960, which(comp_1960$membership == which.max(comp_1960$csize)))
    
    # 3. Compute centrality measures on the LCC
    V(giant_g_1960)$degree      <- degree(giant_g_1960)
    V(giant_g_1960)$betweenness <- betweenness(giant_g_1960, normalized = TRUE)
    V(giant_g_1960)$closeness   <- closeness(giant_g_1960, normalized = TRUE)
    V(giant_g_1960)$eigenvector <- eigen_centrality(giant_g_1960)$vector
    V(giant_g_1960)$pagerank    <- page_rank(giant_g_1960)$vector
    
    # 4. Into a dataframe
    centrality_df_1960 <- data.frame(
      firm_id = V(giant_g_1960)$name,
      firm_name = V(giant_g_1960)$cname,
      degree = V(giant_g_1960)$degree,
      betweenness = V(giant_g_1960)$betweenness,
      closeness = V(giant_g_1960)$closeness,
      eigenvector = V(giant_g_1960)$eigenvector,
      pagerank = V(giant_g_1960)$pagerank,
      stringsAsFactors = FALSE
    )
    
    # 5. Get top 25 nodes for each measure
    top25_degree_1960 <- centrality_df_1960 %>% arrange(desc(degree)) %>% slice_head(n = 25)
    top25_betweenness_1960 <- centrality_df_1960 %>% arrange(desc(betweenness)) %>% slice_head(n = 25)
    top25_closeness_1960 <- centrality_df_1960 %>% arrange(desc(closeness)) %>% slice_head(n = 25)
    top25_eigenvector_1960 <- centrality_df_1960 %>% arrange(desc(eigenvector)) %>% slice_head(n = 25)
    top25_pagerank_1960 <- centrality_df_1960 %>% arrange(desc(pagerank)) %>% slice_head(n = 25)
    top25_all_1960 <- bind_rows(top25_degree_1960, top25_betweenness_1960, top25_closeness_1960, top25_eigenvector_1960, top25_pagerank_1960)
    
    # Add rank and centrality type to each
    top_degree_df_1960 <- top25_degree_1960 %>% mutate(rank = row_number(), centrality = "degree")
    top_bet_df_1960 <- top25_betweenness_1960 %>% mutate(rank = row_number(), centrality = "betweenness")
    top_close_df_1960 <- top25_closeness_1960 %>% mutate(rank = row_number(), centrality = "closeness")
    top_eigen_df_1960 <- top25_eigenvector_1960 %>% mutate(rank = row_number(), centrality = "eigenvector")
    top_pr_df_1960 <- top25_pagerank_1960 %>% mutate(rank = row_number(), centrality = "pagerank")
    
    # Combine all into one table
    top25_all_1960 <- bind_rows(top_degree_df_1960, top_bet_df_1960, top_close_df_1960, top_eigen_df_1960, top_pr_df_1960)
    
    #Wide table
    top25_wide_1960 <- list(
      degree = top25_degree_1960$firm_name,
      betweenness = top25_betweenness_1960$firm_name,
      closeness = top25_closeness_1960$firm_name,
      eigenvector = top25_eigenvector_1960$firm_name,
      pagerank = top25_pagerank_1960$firm_name
    ) %>% as.data.frame()
    
    # View the result
    print(top25_all_1960)
    print(top25_wide_1960)
    
    library(xtable)
    
    print(xtable(top25_all_1960,
                 caption = "Top 25 firms by centrality measures, 1960",
                 label = "tab:top25_all_1960"),
          include.rownames = FALSE,
          file = "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/top25_all_1960.tex",
          booktabs = TRUE)
    
    print(xtable(top25_wide_1960,
                 caption = "Top 25 firms (wide format), 1960",
                 label = "tab:top25_wide_1960"),
          include.rownames = FALSE,
          file = "C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/32_Network_descriptive_paper/results/top25_wide_1960.tex",
          booktabs = TRUE)
    

    
################################################################################    
# Louvain community detection on full firm-to-firm networks
################################################################################    

    for (year in years) {
      g_firm <- get(paste0("g_firms_", year))
      if (gorder(g_firm) > 0 && gsize(g_firm) > 0) {
        comm_firm <- cluster_louvain(g_firm)
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
    
    ############################################################################    
    #Comparing communities and cantonal borders in full network
    
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
        ari <- adjustedRandIndex(V(g)$community[valid], V(g)$ctn[valid])
        
        ari_results <- rbind(
          ari_results,
          data.frame(year = year, adjusted_rand_index = ari)
        )
      }
    }
    
    # View ARI values over time
    print(ari_results)
    
    # plot
    plot(
      ari_results$year,
      ari_results$adjusted_rand_index,
      type = "b", pch = 16,
      xlab = "Year", ylab = "Adjusted Rand Index",
      main = "Similarity Between Firm Communities and Cantons"
    )
    
    ################################################################################    
    # Comparing communities and langage borders in full network
    
    
  
    ################################################################################    
    # Louvain community detection on LCC of firm-to-firm networks
    ################################################################################ 
    
    for (year in years) {
      g_firm <- get(paste0("g_firms_", year))
      if (gorder(g_firm) > 0 && gsize(g_firm) > 0) {
        comm_firm <- cluster_louvain(g_firm)
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
      main = "Similarity Between Firm Communities and Cantons"
    )
    
    ################################################################################    
    # Comparing communities and langage borders in LCC    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


    

    
    
################################################################################    
### Mun-to-Mun aggregations
  
 ## Firm to firm
  # Step 1: Get the edgelist with firm IDs and municipalities
  for (year in years) {
    g_firm <- get(paste0("g_firms_", year))
    
    edge_df <- as_data_frame(g_firm, what = "edges") %>%
      mutate(
        from_mun = V(g_firm)$gdenr[match(from, V(g_firm)$name)],
        to_mun   = V(g_firm)$gdenr[match(to,   V(g_firm)$name)]
      )
    
    assign(paste0("edge_df_", year), edge_df)
  }
  
  # Step 2: Aggregate to municipality-municipality links
  for (year in years) {
    # Get the edge_df with municipality info
    edge_df <- get(paste0("edge_df_", year))
    
    # A. Without self-links
    edge_mun_no_self <- edge_df %>%
      filter(from_mun != to_mun) %>%
      count(from_mun, to_mun, name = "weight")
    
    # B. With self-links (within-municipality links)
    edge_mun_with_self <- edge_df %>%
      count(from_mun, to_mun, name = "weight")
    
    # Assign to global environment
    assign(paste0("edge_mun_no_self_", year), edge_mun_no_self)
    assign(paste0("edge_mun_with_self_", year), edge_mun_with_self)
  }
  
  # Define helper function to rescale firm counts
  rescale_size <- function(x, min_size = 1, max_size = 30) {
    rng <- range(x, na.rm = TRUE)
    if (rng[1] == rng[2]) return(rep(min_size, length(x)))  # avoid division by zero
    scaled <- (x - rng[1]) / (rng[2] - rng[1])
    return(min_size + scaled * (max_size - min_size))
  }
  
  # Step 3 : motherfucking graphs
  for (year in years) {
    # A. Create igraph objects from edge lists
    edge_mun_no_self <- get(paste0("edge_mun_no_self_", year))
    edge_mun_with_self <- get(paste0("edge_mun_with_self_", year))
    
    g_mun_no_self <- graph_from_data_frame(edge_mun_no_self, directed = FALSE)
    g_mun_with_self <- graph_from_data_frame(edge_mun_with_self, directed = FALSE)
    
    assign(paste0("g_mun_no_self_", year), g_mun_no_self)
    assign(paste0("g_mun_with_self_", year), g_mun_with_self)
    
    # B. Count how many firms are located in each municipality
    g_firm <- get(paste0("g_firms_", year))
    firm_counts <- as.data.frame(table(V(g_firm)$gdenr))
    colnames(firm_counts) <- c("gdenr", "firm_count")
    firm_counts$gdenr <- as.character(firm_counts$gdenr)  # match vertex names if needed
    
    # C. Assign firm_count to g_mun_with_self
    g_with <- get(paste0("g_mun_with_self_", year))
    V(g_with)$name <- as.character(V(g_with)$name)
    idx_with <- match(V(g_with)$name, firm_counts$gdenr)
    V(g_with)$firm_count <- firm_counts$firm_count[idx_with]
    V(g_with)$firm_count[is.na(V(g_with)$firm_count)] <- 0
    V(g_with)$size_scaled <- rescale_size(V(g_with)$firm_count)
    assign(paste0("g_mun_with_self_", year), g_with)
    
    # D. Assign firm_count to g_mun_no_self
    g_no <- get(paste0("g_mun_no_self_", year))
    V(g_no)$name <- as.character(V(g_no)$name)
    idx_no <- match(V(g_no)$name, firm_counts$gdenr)
    V(g_no)$firm_count <- firm_counts$firm_count[idx_no]
    V(g_no)$firm_count[is.na(V(g_no)$firm_count)] <- 0
    V(g_no)$size_scaled <- rescale_size(V(g_no)$firm_count)
    assign(paste0("g_mun_no_self_", year), g_no)
  }
  
  for (year in years) {
    # Get the municipality-to-municipality graph
    g <- get(paste0("g_mun_no_self_", year))
    
      # A. Rescale firm count
        V(g)$vertex_size <- rescale_size(V(g)$firm_count)
    
      # B. Add geographical coordinates
        muni_coords_year <- muni_coords_lv95 %>%
        filter(year == year) %>%
        select(gdenr, lon_wgs84, lat_wgs84) %>%
        mutate(gdenr = as.character(gdenr))
  }
  
    # E. Graph with CH layout
 
    for (year in years) {
      g <- get(paste0("g_mun_no_self_", year))
      V(g)$vertex_size <- rescale_size(V(g)$firm_count)
      assign(paste0("g_mun_no_self_", year), g)
    }
    
    for (year in years) {
        coords_year <- muni_coords_lv95 %>%
          filter(year == year) %>%
          select(gdenr, lon_wgs84, lat_wgs84)
        
        assign(paste0("muni_coords_", year), coords_year)
      }
  
    for (year in years) {
      # Retrieve objects
      g <- get(paste0("g_mun_no_self_", year))
      coords <- get(paste0("muni_coords_", year))
      
      # Ensure gdenr codes are character for matching
      V(g)$gdenr <- as.character(V(g)$name)
      coords$gdenr <- as.character(coords$gdenr)
      
      # Match coordinates
      match_idx <- match(V(g)$gdenr, coords$gdenr)
      V(g)$lon <- coords$lon_wgs84[match_idx]
      V(g)$lat <- coords$lat_wgs84[match_idx]
      
      # Layout matrix
      layout_coords <- cbind(V(g)$lon, V(g)$lat)
      
      # Filter nodes with valid coordinates
      valid_nodes <- which(complete.cases(layout_coords))
      g_valid <- induced_subgraph(g, vids = valid_nodes)
      layout_valid <- layout_coords[valid_nodes, ]
      
      # Save results
      assign(paste0("g_mun_no_self_", year), g)
      assign(paste0("g_mun_valid_", year), g_valid)
      assign(paste0("layout_valid_", year), layout_valid)
    }

###
 ## Testing different plots using mun-to-mun
  
    # Plot for 2003 mun-to-mun firm network
    plot(
      g_mun_valid_2003,
      layout = layout_valid_2003,
      vertex.size = V(g_mun_no_self_2003)$vertex_size,
      vertex.label = NA,
      edge.width = E(g_mun_valid_2003)$weight / max(E(g_mun_valid_2003)$weight) * 3,
      edge.color = "gray60",
      main = "Municipality Network (2003, WGS84)"
    )
   
    # Hiding edges with weight < 20
    plot(
      g_mun_valid_2003,
      layout = layout_valid_2003,
      vertex.size = V(g_mun_no_self_2003)$vertex_size,
      vertex.label = NA,
      edge.width = ifelse(E(g_mun_valid_2003)$weight > 20,
                          E(g_mun_valid_2003)$weight / max(E(g_mun_valid_2003)$weight) * 3,
                          0),
      edge.color = ifelse(E(g_mun_valid_2003)$weight > 20, "gray60", NA),
      main = "Municipality Network of Connected Firms (2003, WGS84)"
    )
    
    #Showing top 1% edges wrt weight
    # Compute weight threshold for top 1%
    weights <- E(g_mun_valid_2003)$weight
    threshold <- quantile(weights, 0.99, na.rm = TRUE)  # top 1%
    
    # Build edge.width and edge.color vectors
    edge_width <- ifelse(weights >= threshold, weights / max(weights) * 3, 0)
    edge_color <- ifelse(weights >= threshold, "gray60", NA)
    
    # Plot with only top 1% visible
    plot(
      g_mun_valid_2003,
      layout = layout_valid_2003,
      vertex.size = V(g_mun_no_self_2003)$vertex_size,
      vertex.label = NA,
      edge.width = edge_width,
      edge.color = edge_color,
      main = "Top 1% Inter-Municipality Connections (2003)"
    )
    
### Directors to Directors
 ## Full D-to-D Network
  # Create empty list to store results
    director_network_summaries <- list()
    
    # Loop over years
    for (year in years) {
      # Get the firm–firm graph object
      g_director <- get(paste0("g_directors_", year))
      
      # Compute macro-level stats
      net_sum <- list(
        year = year,
        number_of_edges = gsize(g_director),
        number_of_nodes = gorder(g_director),
        density = edge_density(g_director),
        is_connected = is_connected(g_director),
        number_of_components = components(g_director)$no,
        diameter = if (is_connected(g_director)) diameter(g_director) else NA,
        average_path_length = if (is_connected(g_director)) mean_distance(g_director) else NA,
        transitivity_global = transitivity(g_director, type = "global"),
        assortativity_degree = assortativity_degree(g_director)
      )
      
      # Store
      director_network_summaries[[as.character(year)]] <- net_sum
    }
    
    # Convert to data.frame for summary or plotting
    library(dplyr)
    director_network_summary_df <- bind_rows(director_network_summaries)
    
    # View
    print(director_network_summary_df)
    
 ## In Largest Connected Component (LCC)
    for (year in years) {
      # Get the firm network for the year
      g_director <- get(paste0("g_directors_", year))
      
      # Identify the largest component
      comp <- components(g_director)
      giant_id_dir <- which.max(comp$csize)
      lcc_nodes_dir <- which(comp$membership == giant_id_dir)
      
      # Induce subgraph (largest connected component)
      g_lcc_dir <- induced_subgraph(g_director, vids = lcc_nodes_dir)
      
      # Assign to environment with name like g_firms_lcc_xxxx
      assign(paste0("g_directors_lcc_", year), g_lcc_dir)
    }
    
    # List to store results
    lcc_summary_list_dir <- list()
    
    for (year in years) {
      # Get the LCC graph object from the environment
      g_lcc_dir <- get(paste0("g_directors_lcc_", year))
      
      # Compute macro-level measures
      summary_stats <- list(
        year = year,
        number_of_edges = gsize(g_lcc_dir),
        number_of_nodes = gorder(g_lcc_dir),
        density = edge_density(g_lcc_dir),
        transitivity_global = transitivity(g_lcc_dir, type = "global"),
        assortativity_degree = assortativity_degree(g_lcc_dir)
      )
      
      # Store
      lcc_summary_list_dir[[as.character(year)]] <- summary_stats
    }
    
    # Convert to data frame
    lcc_summary_df_dir <- bind_rows(lcc_summary_list_dir)
    
    # View or export
    print(lcc_summary_df)
  
### Bipartite
 ## Full Bipartite Network
  # Create empty list to store results
    bi_network_summaries <- list()
    
    # Loop over years
    for (year in years) {
      # Get the firm–firm graph object
      g_bi <- get(paste0("g_bi_", year))
      
      # Compute macro-level stats => Search if this make sense for bipartite graphs !
      net_sum <- list(
        year = year,
        number_of_edges = gsize(g_bi),
        number_of_nodes = gorder(g_bi),
        density = edge_density(g_bi),
        is_connected = is_connected(g_bi),
        number_of_components = components(g_bi)$no,
        diameter = if (is_connected(g_bi)) diameter(g_bi) else NA,
        average_path_length = if (is_connected(g_bi)) mean_distance(g_bi) else NA,
        transitivity_global = transitivity(g_bi, type = "global"),
        assortativity_degree = assortativity_degree(g_bi)
      )
      
      # Store
      bi_network_summaries[[as.character(year)]] <- net_sum
    }
    
    # Convert to data.frame for summary or plotting
    library(dplyr)
    bi_network_summary_df <- bind_rows(bi_network_summaries)
    
    # View
    print(bi_network_summary_df)
    
    ## In Largest Connected Component (LCC)
    for (year in years) {
      # Get the firm network for the year
      g_bi <- get(paste0("g_bi_", year))
      
      # Identify the largest component
      comp <- components(g_bi)
      giant_id_bi <- which.max(comp$csize)
      lcc_nodes_bi <- which(comp$membership == giant_id_bi)
      
      # Induce subgraph (largest connected component)
      g_lcc_bi <- induced_subgraph(g_bi, vids = lcc_nodes_bi)
      
      # Assign to environment with name like g_firms_lcc_xxxx
      assign(paste0("g_bi_lcc_", year), g_lcc_bi)
    }
    
    # List to store results
    lcc_summary_list_bi <- list()
    
    for (year in years) {
      # Get the LCC graph object from the environment
      g_lcc_bi <- get(paste0("g_bi_lcc_", year))
      
      # Compute macro-level measures
      summary_stats <- list(
        year = year,
        number_of_edges = gsize(g_lcc_bi),
        number_of_nodes = gorder(g_lcc_bi),
        density = edge_density(g_lcc_bi),
        transitivity_global = transitivity(g_lcc_bi, type = "global"),
        assortativity_degree = assortativity_degree(g_lcc_bi)
      )
      
      # Store
      lcc_summary_list_bi[[as.character(year)]] <- summary_stats
    }
    
    # Convert to data frame
    lcc_summary_df_bi <- bind_rows(lcc_summary_list_bi)
    
    # View or export
    print(lcc_summary_df_bi)

################################################################################
    # -----------------------------
    # Louvain community detection for each year
    # -----------------------------
    for (year in years) {
      g <- get(paste0("g_mun_no_self_", year))
      if (gorder(g) > 0 && gsize(g) > 0) {
        comm <- cluster_louvain(g)
        V(g)$community <- membership(comm)
        assign(paste0("g_mun_no_self_", year), g)
        assign(paste0("comm_mun_", year), comm)
      }
    }
    
    # -----------------------------    
    # Number of communities over time
    # -----------------------------
    
    # Initialize empty data frame
    community_evolution <- data.frame(
      year = integer(),
      num_communities = integer()
    )
    
    # Loop over years and extract number of communities
    for (year in years) {
      comm_name <- paste0("comm_mun_", year)
      if (exists(comm_name)) {
        comm_obj <- get(comm_name)
        num_comm <- length(unique(membership(comm_obj)))
        community_evolution <- rbind(community_evolution, data.frame(year = year, num_communities = num_comm))
      }
    }
    
    # Display results
    print(community_evolution)
    
    # Optional: plot
    plot(
      community_evolution$year,
      community_evolution$num_communities,
      type = "b", pch = 16,
      xlab = "Year", ylab = "Number of communities",
      main = "Evolution of Municipality Network Communities (Louvain)"
    )
    
    # --------------------------------------------      
    # Visualization of mun. communities over CH layout
    # --------------------------------------------

    # Generate a distinct color per community
    n_comm <- length(unique(V(g_mun_no_self_1960)$community))
    pal <- rainbow(n_comm, alpha = 0.7)  # semi-transparent colors
    comm_colors <- pal[V(g_mun_no_self_1960)$community]
    
    n_comm <- length(unique(V(g_mun_no_self_1992)$community))
    pal <- rainbow(n_comm, alpha = 0.7)  # semi-transparent colors
    comm_colors <- pal[V(g_mun_no_self_1992)$community]
    
    n_comm <- length(unique(V(g_mun_no_self_2003)$community))
    pal <- rainbow(n_comm, alpha = 0.7)  # semi-transparent colors
    comm_colors <- pal[V(g_mun_no_self_2003)$community]
    
    # Plot using geographic layout
    
    plot(
      g_mun_no_self_1960,
      layout = layout_valid_1960,
      vertex.color = comm_colors,
      vertex.size = V(g_mun_no_self_1960)$vertex_size,
      vertex.label = NA,
      edge.width = 0.5,
      edge.color = "gray85",
      main = "Louvain Communities in Municipality Network (1960)"
    )
    
    plot(
      g_mun_no_self_1992,
      layout = layout_valid_1992,
      vertex.color = comm_colors,
      vertex.size = V(g_mun_no_self_1992)$vertex_size,
      vertex.label = NA,
      edge.width = 0.5,
      edge.color = "gray85",
      main = "Louvain Communities in Municipality Network (1992)"
    )
    
    plot(
      g_mun_no_self_2003,
      layout = layout_valid_2003,
      vertex.color = comm_colors,
      vertex.size = V(g_mun_no_self_2003)$vertex_size,
      vertex.label = NA,
      edge.width = 0.5,
      edge.color = "gray85",
      main = "Louvain Communities in Municipality Network (2003)"
    )

    
    
    # Louvain community detection on director-to-director networks
    # -----------------------------
    for (year in years) {
      g_dir <- get(paste0("g_directors_", year))
      if (gorder(g_dir) > 0 && gsize(g_dir) > 0) {
        comm_dir <- cluster_louvain(g_dir)
        V(g_dir)$community <- membership(comm_dir)
        assign(paste0("g_directors_", year), g_dir)
        assign(paste0("comm_directors_", year), comm_dir)
      }
    }
    
    # -----------------------------
    # Evolution of director communities (normalized)
    # -----------------------------
    director_community_evolution <- data.frame(
      year = integer(),
      num_directors = integer(),
      num_communities = integer(),
      communities_per_1000_directors = numeric()
    )
    
    for (year in years) {
      comm_name <- paste0("comm_directors_", year)
      g_name <- paste0("g_directors_", year)
      
      if (exists(comm_name) && exists(g_name)) {
        comm <- get(comm_name)
        g <- get(g_name)
        
        num_comm <- length(unique(membership(comm)))
        num_nodes <- gorder(g)
        comm_per_1000 <- num_comm / num_nodes * 1000
        
        director_community_evolution <- rbind(
          director_community_evolution,
          data.frame(
            year = year,
            num_directors = num_nodes,
            num_communities = num_comm,
            communities_per_1000_directors = comm_per_1000
          )
        )
      }
    }
    
    # -----------------------------
    # ARI between director communities and cantons
    # -----------------------------
    if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
    library(mclust)
    
    ari_results_directors <- data.frame(
      year = integer(),
      adjusted_rand_index = numeric()
    )
    
    for (year in years) {
      g_name <- paste0("g_directors_", year)
      g <- get(g_name)
      
      if (!is.null(V(g)$community) && !is.null(V(g)$ctn)) {
        valid <- !is.na(V(g)$community) & !is.na(V(g)$ctn)
        ari <- adjustedRandIndex(V(g)$community[valid], V(g)$ctn[valid])
        
        ari_results_directors <- rbind(
          ari_results_directors,
          data.frame(year = year, adjusted_rand_index = ari)
        )
      }
    }
    
    # Plotting results
    plot(
      director_community_evolution$year,
      director_community_evolution$communities_per_1000_directors,
      type = "b", pch = 16,
      xlab = "Year", ylab = "Communities per 1000 directors",
      main = "Normalized Director Community Count (Louvain)"
    )
    
    plot(
      ari_results_directors$year,
      ari_results_directors$adjusted_rand_index,
      type = "b", pch = 16,
      xlab = "Year", ylab = "Adjusted Rand Index",
      main = "Similarity Between Director Communities and Cantons"
    )    # -----------------------------

    
    # -----------------------------
    # ARI between director communities and langage
    # -----------------------------
    
    # Create data frame to store ARI results
    ari_results_lang_directors <- data.frame(
      year = integer(),
      adjusted_rand_index = numeric()
    )
    
    for (year in years) {
      g_name <- paste0("g_directors_", year)
      g <- get(g_name)
      
      if (!is.null(V(g)$community) && !is.null(V(g)$language_w)) {
        valid <- !is.na(V(g)$community) & !is.na(V(g)$language_w)
        ari <- adjustedRandIndex(
          V(g)$community[valid],
          V(g)$language_w[valid]
        )
        
        ari_results_lang_directors <- rbind(
          ari_results_lang_directors,
          data.frame(year = year, adjusted_rand_index = ari)
        )
      }
    }
    
    # Plot 
    plot(
      ari_results_lang_directors$year,
      ari_results_lang_directors$adjusted_rand_index,
      type = "b", pch = 16,
      xlab = "Year", ylab = "Adjusted Rand Index",
      main = "Similarity Between Director Communities and Language Regions"
    )
    
    # Replace with your director-based municipal graph

    top_dir_edges <- as_data_frame(g_directors_2000, what = "edges") %>%
      mutate(
        from_lang = V(g_directors_2000)$language_w[match(from, V(g_directors_2000)$name)],
        to_lang   = V(g_directors_2000)$language_w[match(to,   V(g_directors_2000)$name)],
        lang_relation = ifelse(from_lang == to_lang, "Same Language", "Different Language")
      ) %>%
      arrange(desc(weight)) %>%
      slice(1:20) %>%
      select(from, to, weight, lang_relation)
    
    print(top_dir_edges)
    
    
################################################################################    
# 5. Firm Network Metrics

  # Metrics
  degree(g_firms_1943)
  degree_distribution(g_firms_1943, cumulative = TRUE)
  plot(degree_distribution(g_firms_1943, cumulative = TRUE), col="orange")
  hist(degree(g_firms_1943))
  abline(v=mean(degree(g_firms_1943)), col="red")
  mean(degree(g_firms_1943))
################################################################################


# 7. Community Detection test

  k1_1943 = cluster_walktrap(g_mun_no_self_1943)
  k2_1943 = cluster_infomap(g_mun_no_self_1943) #M. Rosvall and C. T. Bergstrom : community structure that minimizes the expected description length of a random walker trajectory.
  #k3_1943 = cluster_edge_betweenness(g_mun_no_self_1943) # Girvan-Newman algorithm : betweenness of the edges connecting two communities is typically high
  k5_1943 = cluster_fast_greedy(g_mun_no_self_1943) #see A Clauset, MEJ Newman, C Moore: Finding community structure in very large networks, http://www.arxiv.org/abs/cond-mat/0408187 for the details.
  
