install.packages(c("readstata13", "openxlsx", "devtools", "igraph", "bipartite", "asnipe", "assortnet", "ggplot2", "ggmap", "rnetcarto", "ecodist", "igraphdata", "statnet", "RColorBrewer", "tidyverse", "xlsx"))
#?graph_from_adjacency_matrix

library(openxlsx)
library(haven)
library(dplyr)
library(igraph) #load the package
install.packages("readstata13")
library(readstata13)
library(xlsx)
library(sf)

install.packages(c("tidygraph", "ggraph"))
install.packages("rnaturalearth")
library(ggraph)
library(ggplot2)
library(tidygraph)  # Needed for ggraph compatibility

BiocManager::install("RCy3")
library(RCy3)

setwd("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/14_Network_OneMode")

################################################################################
# 1943 : New code (29.04.2025), to update, I am working with 1943
## Creation of bipartite network, firm-to-firm, director-to-director, mun-to-mun : all years

#Store all years
years <- c(1934, 1943, 1960, 1962, 1963, 1964, 1965, 1966, 1969,
           1972, 1975, 1979, 1980, 1981, 1982, 1983, 1984, 1985,
           1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994,
           1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003)

for (year in years) {
  path <- paste0("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_", year, "_bipartite.dta")
  edgelist <- read_stata(path)
  assign(paste0("edgeliste_", year, "_bi"), edgelist)
}

## Macro-level description of network over time

edgeliste_1943_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1943_bipartite.dta")

# 1. Create bipartite graph from edge list
g_bi_1943 <- graph_from_data_frame(edgeliste_1943_bi, directed = FALSE)

# 2. Set 'type' attribute: TRUE for directors, FALSE for firms
# (igraph needs to know which nodes are which type)
V(g_bi_1943)$type <- bipartite_mapping(g_bi_1943)$type
projections_1943 <- bipartite_projection(g_bi_1943)
g_firms_1943 <- projections_1943$proj1  # This will be the firm-firm projection
g_directors_1943 <- projections_1943$proj2  # This will be the firm-firm projection

# 3. Add characteristics and generate mun-to-mun network

# 3.1. Firms
firms_charac_1943 <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/Node_attribute_comp/comp_1943_attribute_complete.dta")

# 3.1.1. Create the index matching graph nodes to rows in the attribute dataframe
match_idx_firms_1943 <- match(V(g_firms_1943)$name, firms_charac_1943$N_UCID)

# 3.1.2. Add attributes to the graph
#Municipality nb
V(g_firms_1943)$gdenr <- firms_charac_1943$gdenr[match_idx_firms_1943]
V(g_firms_1943)$gdenr_2012 <- firms_charac_1943$gdenr_2012[match_idx_firms_1943]
V(g_firms_1943)$gdenr_2018 <- firms_charac_1943$gdenr_2018[match_idx_firms_1943]

#Others ? Capital, Owners, Share female, etc.

# 3.2. Directors
directors_charac_1943 <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/Node_attribute_pers/pers_1943_attribute_complete.dta")
# 3.2.1. Create the index matching graph nodes to rows in the attribute dataframe
match_idx_directors_1943 <- match(V(g_directors_1943)$name,directors_charac_1943$N_id_sug)

# 3.2.2. Add attributes to the graph
#Municipality nb
V(g_directors_1943)$gdenr_2018 <- directors_charac_1943$gdenr_2018[match_idx_directors_1943]

#Others ? Capital, Owners, Share female, etc.

# 3.3. Municipality Coordinates
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


#muni_coords_1943 <- muni_coords %>%
#  filter(year == 1943) %>%
#  select(gdenr, E_CNTR, N_CNTR)

# 4. Bipartite Network Graphs and Metrics
# 4.1 Graph
plot(g_bi_1943, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1943 : Firm Networks and interlocks")

# 4.2 Metrics
degree(g_bi_1943)
degree_distribution(g_bi_1943, cumulative = TRUE)
plot(degree_distribution(g_bi_1943, cumulative = TRUE), col="orange")
hist(degree(g_bi_1943))
abline(v=mean(degree(g_bi_1943)), col="red")
mean(degree(g_bi_1943))

# 5. Firm Network Graphs and Metrics
# 5.1 Graph
plot(g_firms_1943, layout=layout_with_drl, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1943 : Firm Networks")

# 5.2 Metrics
degree(g_firms_1943)
degree_distribution(g_firms_1943, cumulative = TRUE)
plot(degree_distribution(g_firms_1943, cumulative = TRUE), col="orange")
hist(degree(g_firms_1943))
abline(v=mean(degree(g_firms_1943)), col="red")
mean(degree(g_firms_1943))

# 6. Directors Network Graphs and Metrics
# 6.1 Graph
plot(g_directors_1943, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1943 : Firm Networks")

# 6.2 Degree distribution
degree(g_directors_1943)
degree_distribution(g_directors_1943, cumulative = TRUE)
plot(degree_distribution(g_directors_1943, cumulative = TRUE), col="orange")
hist(degree(g_directors_1943))
abline(v=mean(degree(g_directors_1943)), col="red")
mean(degree(g_directors_1943))

# 7. Tests

# 7.1 Macro-level network summary of firm network
network_summary_1943 <- list(
  number_of_edges = gsize(g_firms_1943),
  number_of_nodes = gorder(g_firms_1943),
  density = edge_density(g_firms_1943),
  is_connected = is_connected(g_firms_1943),
  number_of_components = components(g_firms_1943)$no,
  diameter = if (is_connected(g_firms_1943)) diameter(g_firms_1943) else NA,
  average_path_length = if (is_connected(g_firms_1943)) mean_distance(g_firms_1943) else NA,
  transitivity_global = transitivity(g_firms_1943, type = "global"), #This is simply the ratio of the count of triangles and connected triples in the graph.
  assortativity_degree = assortativity_degree(g_firms_1943),
  reciprocity = if (is_directed(g_firms_1943)) reciprocity(g_firms_1943) else NA
)

# Print the summary
print(network_summary_1943)

###Other measures
#Average nearest neighbor degree: mean of degree from direct neighbors (how connected are your neighbors)
knn(
  g_firms_1943,
  neighbor.degree.mode = c("all"),
  weights = NULL
) #unweighted

#Transitivity of graph
transitivity_local_1943 = transitivity(g_firms_1943, type = "local") #The local transitivity of a vertex is the ratio of the count of triangles connected to the vertex and the triples centered on the vertex.


# 7.2 Random Graphs

#Random Graphs with nodes and edges charac (Erdős-Rényi model)
ER_graph_1943 <-sample_gnm(network_summary_1943$number_of_nodes, network_summary_1943$number_of_edges, directed = FALSE, loops = FALSE)
#graph1 <- sample_grg(network_summary_1943$number_of_nodes, radius?!)
plot(ER_graph_1943, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1943 : Random Graph")

network_summary_random_1943 <- list(
  number_of_edges = gsize(ER_graph_1943),
  number_of_nodes = gorder(ER_graph_1943),
  density = edge_density(ER_graph_1943),
  is_connected = is_connected(ER_graph_1943),
  number_of_components = components(ER_graph_1943)$no,
  diameter = if (is_connected(ER_graph_1943)) diameter(ER_graph_1943) else NA,
  average_path_length = if (is_connected(ER_graph_1943)) mean_distance(ER_graph_1943) else NA,
  transitivity_global = transitivity(ER_graph_1943, type = "global"),
  assortativity_degree = assortativity_degree(ER_graph_1943),
)
print(network_summary_random_1943)

degree_distribution(ER_graph_1943, cumulative = TRUE)
plot(degree_distribution(ER_graph_1943, cumulative = TRUE), col="orange")
hist(degree(ER_graph_1943))
abline(v=mean(degree(ER_graph_1943)), col="red")
mean(degree(ER_graph_1943))

# 7.3 Community Detection

k1_1943 = cluster_walktrap(g_firms_1943)
k2_1943 = cluster_infomap(g_firms_1943) #M. Rosvall and C. T. Bergstrom : community structure that minimizes the expected description length of a random walker trajectory.
#k3_1943 = cluster_edge_betweenness(g_firms_1943) # Girvan-Newman algorithm : betweenness of the edges connecting two communities is typically high
k4_1943 = cluster_louvain(g_firms_1943) #see VD Blondel, J-L Guillaume, R Lambiotte and E Lefebvre: Fast unfolding of community hierarchies in large networks, https://arxiv.org/abs/0803.0476 for the details.
k5_1943 = cluster_fast_greedy(g_firms_1943) #see A Clauset, MEJ Newman, C Moore: Finding community structure in very large networks, http://www.arxiv.org/abs/cond-mat/0408187 for the details.

#plot(k4_1943,g_firms_1943)
k4_1943$membership

# 7.4 Centrality at node-level

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
  degree = top25_degree_1943$firm_id,
  betweenness = top25_betweenness_1943$firm_id,
  closeness = top25_closeness_1943$firm_id,
  eigenvector = top25_eigenvector_1943$firm_id,
  pagerank = top25_pagerank_1943$firm_id
) %>% as.data.frame()

# View the result
print(top25_all_1943)
print(top25_wide_1943)

################################################################################
# Mun-to-Mun network

# Step 1: Get the edgelist with firm IDs and municipalities
edge_df_1943 <- get.data.frame(g_firms_1943, what = "edges") %>%
  mutate(
    from_mun = V(g_firms_1943)$gdenr[match(from, V(g_firms_1943)$name)],
    to_mun   = V(g_firms_1943)$gdenr[match(to, V(g_firms_1943)$name)]
  )

# Step 2: Aggregate to municipality-municipality links

# A. Without self-links
edge_mun_no_self_1943 <- edge_df_1943 %>%
  filter(from_mun != to_mun) %>%
  count(from_mun, to_mun, name = "weight")

# B. With self-links (within-municipality links)
edge_mun_with_self_1943 <- edge_df_1943 %>%
  count(from_mun, to_mun, name = "weight")

# Step 3 : Graphs
# A. Version A: No self-links
g_mun_no_self_1943 <- graph_from_data_frame(edge_mun_no_self_1943, directed = FALSE)
g_mun_with_self_1943 <- graph_from_data_frame(edge_mun_with_self_1943, directed = FALSE)

# B. Basic Visualization
plot(g_mun_no_self_1943, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1943 : Mun-to-Mun, no self-loop")
plot(g_mun_with_self_1943, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1943 : Mun-to-Mun, with self-loop")

# C.Count how many firms are in each municipality
firm_counts_1943 <- as.data.frame(table(V(g_firms_1943)$gdenr))
colnames(firm_counts_1943) <- c("gdenr", "firm_count")

# D. Self loop graph
# Match municipality names in the graph to the count table
match_idx_mun_1943 <- match(V(g_mun_with_self_1943)$name, firm_counts_1943$gdenr)

# Assign firm counts to nodes
V(g_mun_with_self_1943)$firm_count <- firm_counts_1943$firm_count[match_idx_mun_1943]

# Deal with NAs
V(g_mun_with_self_1943)$firm_count[is.na(V(g_mun_with_self_1943)$firm_count)] <- 0


plot(
  g_mun_with_self_1943,
  vertex.size = sqrt(V(g_mun_with_self_1943)$firm_count + 1),  # square root to avoid huge size disparities
  vertex.label = NA,
  main = "Municipality Network"
)

# E. No self loop graph
# Match municipality names in the graph to the count table
match_idx_mun_1943 <- match(V(g_mun_no_self_1943)$name, firm_counts_1943$gdenr)

# Assign firm counts to nodes
V(g_mun_no_self_1943)$firm_count <- firm_counts_1943$firm_count[match_idx_mun_1943]

#Visualization
plot(
  g_mun_no_self_1943,
  vertex.size = sqrt(V(g_mun_no_self_1943)$firm_count + 1),  # square root to avoid huge size disparities
  vertex.label = NA,
  main = "Municipality Network"
)

# Rescale firm count to range [1, 30]
rescale_size <- function(x, min_size = 1, max_size = 30) {
  rng <- range(x, na.rm = TRUE)
  scaled <- (x - rng[1]) / (rng[2] - rng[1])
  return(min_size + scaled * (max_size - min_size))
}

V(g_mun_no_self_1943)$vertex_size <- rescale_size(V(g_mun_no_self_1943)$firm_count)

plot(
  g_mun_no_self_1943,
  vertex.size = V(g_mun_no_self_1943)$vertex_size,
  vertex.label = NA,
  main = "Municipality Network (normalized size)"
)

# Add geographical info
muni_coords_1943 <- muni_coords_lv95 %>%
  filter(year == 1943) %>%
  select(gdenr, lon_wgs84, lat_wgs84)


# Ensure municipality codes are strings for matching
V(g_mun_no_self_1943)$gdenr <- as.character(V(g_mun_no_self_1943)$name)
muni_coords_1943$gdenr <- as.character(muni_coords_1943$gdenr)

# Match coordinates to nodes in the graph
match_idx_mun_1943 <- match(V(g_mun_no_self_1943)$gdenr, muni_coords_1943$gdenr)

# Assign coordinates (with NA for unmatched ones)
V(g_mun_no_self_1943)$lon <- muni_coords_1943$lon_wgs84[match_idx_mun_1943]
V(g_mun_no_self_1943)$lat <- muni_coords_1943$lat_wgs84[match_idx_mun_1943]

#Plot that shit
  layout_coords <- cbind(V(g_mun_no_self_1943)$lon, V(g_mun_no_self_1943)$lat)
  
  # Filter out nodes without valid coordinates
  valid_nodes <- which(complete.cases(layout_coords))
  g_mun_valid <- induced_subgraph(g_mun_no_self_1943, vids = valid_nodes)
  layout_valid <- layout_coords[valid_nodes, ]
  
  # Plot
  plot(
    g_mun_valid,
    layout = layout_valid,
    vertex.size = V(g_mun_no_self_1943)$vertex_size,
    vertex.label = NA,
    edge.width = E(g_mun_valid)$weight / max(E(g_mun_valid)$weight) * 3,
    edge.color = "gray60",
    main = "Municipality Network (1943, WGS84)"
  )

  
################################################################################
################################################################################
################################################################################

# 1943 : Old code, keep it for export of xlsx and csv

write.xlsx(pers_unimode_weight_1943, file = "pers_unimode_w_1943.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode_1943, file = "pers_unimode_1943.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight_1943, file = "pers_unimode_w_1943.csv")
#write.csv(pers_unimode_1943, file = "pers_unimode_1943.csv")

pers_eg_1943=graph_from_data_frame(pers_unimode_1943, directed=FALSE)
pers_eg_1943
pers_el <- as_edgelist(pers_eg_1943)



################################################################################
# Loop over other years : 

edgeliste_1960_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1960_bipartite.dta")

# 1. Create bipartite graph from edge list
g_bi_1960 <- graph_from_data_frame(edgeliste_1960_bi, directed = FALSE)

# 2. Set 'type' attribute: TRUE for directors, FALSE for firms
# (igraph needs to know which nodes are which type)
V(g_bi_1960)$type <- bipartite_mapping(g_bi_1960)$type
projections_1960 <- bipartite_projection(g_bi_1960)
g_firms_1960 <- projections_1960$proj1  # This will be the firm-firm projection
g_directors_1960 <- projections_1960$proj2  # This will be the firm-firm projection

# 3. Add characteristics and generate mun-to-mun network

# 3.1. Firms
firms_charac_1960 <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/Node_attribute_comp/comp_1960_attribute_complete.dta")

# 3.1.1. Create the index matching graph nodes to rows in the attribute dataframe
match_idx_firms_1960 <- match(V(g_firms_1960)$name, firms_charac_1960$N_UCID)

# 3.1.2. Add attributes to the graph
#Municipality nb
V(g_firms_1960)$gdenr <- firms_charac_1960$gdenr[match_idx_firms_1960]
V(g_firms_1960)$gdenr_2012 <- firms_charac_1960$gdenr_2012[match_idx_firms_1960]
V(g_firms_1960)$gdenr_2018 <- firms_charac_1960$gdenr_2018[match_idx_firms_1960]

#Others ? Capital, Owners, Share female, etc.

# 3.2. Directors
directors_charac_1960 <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/Node_attribute_pers/pers_1960_attribute_complete.dta")
# 3.2.1. Create the index matching graph nodes to rows in the attribute dataframe
match_idx_directors_1960 <- match(V(g_directors_1960)$name,directors_charac_1960$N_id_sug)

# 3.2.2. Add attributes to the graph
#Municipality nb
V(g_directors_1960)$gdenr_2018 <- directors_charac_1960$gdenr_2018[match_idx_directors_1960]

#Others ? Capital, Owners, Share female, etc.

# 3.3. Municipality Coordinates
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


#muni_coords_1960 <- muni_coords %>%
#  filter(year == 1960) %>%
#  select(gdenr, E_CNTR, N_CNTR)

# 4. Bipartite Network Graphs and Metrics
# 4.1 Graph
plot(g_bi_1960, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1960 : Firm Networks and interlocks")

# 4.2 Metrics
degree(g_bi_1960)
degree_distribution(g_bi_1960, cumulative = TRUE)
plot(degree_distribution(g_bi_1960, cumulative = TRUE), col="orange")
hist(degree(g_bi_1960))
abline(v=mean(degree(g_bi_1960)), col="red")
mean(degree(g_bi_1960))

# 5. Firm Network Graphs and Metrics
# 5.1 Graph
plot(g_firms_1960, layout=layout_with_drl, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1960 : Firm Networks")

# 5.2 Metrics
degree(g_firms_1960)
degree_distribution(g_firms_1960, cumulative = TRUE)
plot(degree_distribution(g_firms_1960, cumulative = TRUE), col="orange")
hist(degree(g_firms_1960))
abline(v=mean(degree(g_firms_1960)), col="red")
mean(degree(g_firms_1960))

# 6. Directors Network Graphs and Metrics
# 6.1 Graph
plot(g_directors_1960, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1960 : Firm Networks")

# 6.2 Degree distribution
degree(g_directors_1960)
degree_distribution(g_directors_1960, cumulative = TRUE)
plot(degree_distribution(g_directors_1960, cumulative = TRUE), col="orange")
hist(degree(g_directors_1960))
abline(v=mean(degree(g_directors_1960)), col="red")
mean(degree(g_directors_1960))

# 7. Tests

# 7.1 Macro-level network summary of firm network
network_summary_1960 <- list(
  number_of_edges = gsize(g_firms_1960),
  number_of_nodes = gorder(g_firms_1960),
  density = edge_density(g_firms_1960),
  is_connected = is_connected(g_firms_1960),
  number_of_components = components(g_firms_1960)$no,
  diameter = if (is_connected(g_firms_1960)) diameter(g_firms_1960) else NA,
  average_path_length = if (is_connected(g_firms_1960)) mean_distance(g_firms_1960) else NA,
  transitivity_global = transitivity(g_firms_1960, type = "global"), #This is simply the ratio of the count of triangles and connected triples in the graph.
  assortativity_degree = assortativity_degree(g_firms_1960),
  reciprocity = if (is_directed(g_firms_1960)) reciprocity(g_firms_1960) else NA
)

# Print the summary
print(network_summary_1960)

###Other measures
#Average nearest neighbor degree: mean of degree from direct neighbors (how connected are your neighbors)
knn(
  g_firms_1960,
  neighbor.degree.mode = c("all"),
  weights = NULL
) #unweighted

#Transitivity of graph
transitivity_local_1960 = transitivity(g_firms_1960, type = "local") #The local transitivity of a vertex is the ratio of the count of triangles connected to the vertex and the triples centered on the vertex.


# 7.2 Random Graphs

#Random Graphs with nodes and edges charac (Erdős-Rényi model)
ER_graph_1960 <-sample_gnm(network_summary_1960$number_of_nodes, network_summary_1960$number_of_edges, directed = FALSE, loops = FALSE)
#graph1 <- sample_grg(network_summary_1960$number_of_nodes, radius?!)
plot(ER_graph_1960, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1960 : Random Graph")

network_summary_random_1960 <- list(
  number_of_edges = gsize(ER_graph_1960),
  number_of_nodes = gorder(ER_graph_1960),
  density = edge_density(ER_graph_1960),
  is_connected = is_connected(ER_graph_1960),
  number_of_components = components(ER_graph_1960)$no,
  diameter = if (is_connected(ER_graph_1960)) diameter(ER_graph_1960) else NA,
  average_path_length = if (is_connected(ER_graph_1960)) mean_distance(ER_graph_1960) else NA,
  transitivity_global = transitivity(ER_graph_1960, type = "global"),
  assortativity_degree = assortativity_degree(ER_graph_1960),
)
print(network_summary_random_1960)

degree_distribution(ER_graph_1960, cumulative = TRUE)
plot(degree_distribution(ER_graph_1960, cumulative = TRUE), col="orange")
hist(degree(ER_graph_1960))
abline(v=mean(degree(ER_graph_1960)), col="red")
mean(degree(ER_graph_1960))

# 7.3 Community Detection

k1_1960 = cluster_walktrap(g_firms_1960)
k2_1960 = cluster_infomap(g_firms_1960) #M. Rosvall and C. T. Bergstrom : community structure that minimizes the expected description length of a random walker trajectory.
#k3_1960 = cluster_edge_betweenness(g_firms_1960) # Girvan-Newman algorithm : betweenness of the edges connecting two communities is typically high
k4_1960 = cluster_louvain(g_firms_1960) #see VD Blondel, J-L Guillaume, R Lambiotte and E Lefebvre: Fast unfolding of community hierarchies in large networks, https://arxiv.org/abs/0803.0476 for the details.
k5_1960 = cluster_fast_greedy(g_firms_1960) #see A Clauset, MEJ Newman, C Moore: Finding community structure in very large networks, http://www.arxiv.org/abs/cond-mat/0408187 for the details.

#plot(k4_1960,g_firms_1960)
k4_1960$membership

# 7.4 Centrality at node-level

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
  degree = top25_degree_1960$firm_id,
  betweenness = top25_betweenness_1960$firm_id,
  closeness = top25_closeness_1960$firm_id,
  eigenvector = top25_eigenvector_1960$firm_id,
  pagerank = top25_pagerank_1960$firm_id
) %>% as.data.frame()

# View the result
print(top25_all_1960)
print(top25_wide_1960)

################################################################################
# Mun-to-Mun network

# Step 1: Get the edgelist with firm IDs and municipalities
edge_df_1960 <- get.data.frame(g_firms_1960, what = "edges") %>%
  mutate(
    from_mun = V(g_firms_1960)$gdenr[match(from, V(g_firms_1960)$name)],
    to_mun   = V(g_firms_1960)$gdenr[match(to, V(g_firms_1960)$name)]
  )

# Step 2: Aggregate to municipality-municipality links

# A. Without self-links
edge_mun_no_self_1960 <- edge_df_1960 %>%
  filter(from_mun != to_mun) %>%
  count(from_mun, to_mun, name = "weight")

# B. With self-links (within-municipality links)
edge_mun_with_self_1960 <- edge_df_1960 %>%
  count(from_mun, to_mun, name = "weight")

# Step 3 : Graphs
# A. Version A: No self-links
g_mun_no_self_1960 <- graph_from_data_frame(edge_mun_no_self_1960, directed = FALSE)
g_mun_with_self_1960 <- graph_from_data_frame(edge_mun_with_self_1960, directed = FALSE)

# B. Basic Visualization
plot(g_mun_no_self_1960, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1960 : Mun-to-Mun, no self-loop")
plot(g_mun_with_self_1960, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1960 : Mun-to-Mun, with self-loop")

# C.Count how many firms are in each municipality
firm_counts_1960 <- as.data.frame(table(V(g_firms_1960)$gdenr))
colnames(firm_counts_1960) <- c("gdenr", "firm_count")

# D. Self loop graph
# Match municipality names in the graph to the count table
match_idx_mun_1960 <- match(V(g_mun_with_self_1960)$name, firm_counts_1960$gdenr)

# Assign firm counts to nodes
V(g_mun_with_self_1960)$firm_count <- firm_counts_1960$firm_count[match_idx_mun_1960]

# Deal with NAs
V(g_mun_with_self_1960)$firm_count[is.na(V(g_mun_with_self_1960)$firm_count)] <- 0


plot(
  g_mun_with_self_1960,
  vertex.size = sqrt(V(g_mun_with_self_1960)$firm_count + 1),  # square root to avoid huge size disparities
  vertex.label = NA,
  main = "Municipality Network"
)

# E. No self loop graph
# Match municipality names in the graph to the count table
match_idx_mun_1960 <- match(V(g_mun_no_self_1960)$name, firm_counts_1960$gdenr)

# Assign firm counts to nodes
V(g_mun_no_self_1960)$firm_count <- firm_counts_1960$firm_count[match_idx_mun_1960]

#Visualization
plot(
  g_mun_no_self_1960,
  vertex.size = sqrt(V(g_mun_no_self_1960)$firm_count + 1),  # square root to avoid huge size disparities
  vertex.label = NA,
  main = "Municipality Network"
)

# Rescale firm count to range [1, 30]
rescale_size <- function(x, min_size = 1, max_size = 30) {
  rng <- range(x, na.rm = TRUE)
  scaled <- (x - rng[1]) / (rng[2] - rng[1])
  return(min_size + scaled * (max_size - min_size))
}

V(g_mun_no_self_1960)$vertex_size <- rescale_size(V(g_mun_no_self_1960)$firm_count)

plot(
  g_mun_no_self_1960,
  vertex.size = V(g_mun_no_self_1960)$vertex_size,
  vertex.label = NA,
  main = "Municipality Network (normalized size)"
)

# Add geographical info
muni_coords_1960 <- muni_coords_lv95 %>%
  filter(year == 1960) %>%
  select(gdenr, lon_wgs84, lat_wgs84)


# Ensure municipality codes are strings for matching
V(g_mun_no_self_1960)$gdenr <- as.character(V(g_mun_no_self_1960)$name)
muni_coords_1960$gdenr <- as.character(muni_coords_1960$gdenr)

# Match coordinates to nodes in the graph
match_idx_mun_1960 <- match(V(g_mun_no_self_1960)$gdenr, muni_coords_1960$gdenr)

# Assign coordinates (with NA for unmatched ones)
V(g_mun_no_self_1960)$lon <- muni_coords_1960$lon_wgs84[match_idx_mun_1960]
V(g_mun_no_self_1960)$lat <- muni_coords_1960$lat_wgs84[match_idx_mun_1960]

#Plot that shit
layout_coords <- cbind(V(g_mun_no_self_1960)$lon, V(g_mun_no_self_1960)$lat)

# Filter out nodes without valid coordinates
valid_nodes <- which(complete.cases(layout_coords))
g_mun_valid <- induced_subgraph(g_mun_no_self_1960, vids = valid_nodes)
layout_valid <- layout_coords[valid_nodes, ]

# Plot
plot(
  g_mun_valid,
  layout = layout_valid,
  vertex.size = V(g_mun_no_self_1960)$vertex_size,
  vertex.label = NA,
  edge.width = E(g_mun_valid)$weight / max(E(g_mun_valid)$weight) * 3,
  edge.color = "gray60",
  main = "Municipality Network (1960, WGS84)"
)


################################################################################
################################################################################