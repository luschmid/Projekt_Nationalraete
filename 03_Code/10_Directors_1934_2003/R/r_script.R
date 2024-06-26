install.packages(c("openxlsx", "devtools", "igraph", "bipartite", "asnipe", "assortnet", "ggplot2", "ggmap", "rnetcarto", "ecodist", "igraphdata", "statnet", "RColorBrewer", "tidyverse", "xlsx"))
#?graph_from_adjacency_matrix
setwd("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/14_Network_OneMode")

library(openxlsx)
library(haven)
library(igraph) #load the package
install.packages("readstata13")
library(readstata13)
library(xlsx)

# 1934

edgeliste_1934_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1934_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1934_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1934_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1934.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1934.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1934.csv")
#write.csv(pers_unimode, file = "pers_unimode_1934.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1934_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1934.xlsx", sheetName = "firms",
            colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1934.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1934.csv")
#write.csv(firm_unimode, file = "firm_unimode_1934.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1943

edgeliste_1943_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1943_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1943_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1943_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1943.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1943.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1943.csv")
#write.csv(pers_unimode, file = "pers_unimode_1943.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1943_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1943.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1943.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1943.csv")
#write.csv(firm_unimode, file = "firm_unimode_1943.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1960

edgeliste_1960_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1960_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1960_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1960_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1960.xlsx", sheetName = "pers",
            colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1960.xlsx", sheetName = "pers",
            colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1960.csv")
#write.csv(pers_unimode, file = "pers_unimode_1960.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1960_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1960.xlsx", sheetName = "firms",
            colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1960.xlsx", sheetName = "firms",
            colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1960.csv")
#write.csv(firm_unimode, file = "firm_unimode_1960.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1962

edgeliste_1962_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1962_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1962_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1962_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1962.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1962.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1962.csv")
#write.csv(pers_unimode, file = "pers_unimode_1962.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1962_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1962.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1962.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1962.csv")
#write.csv(firm_unimode, file = "firm_unimode_1962.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1963 

edgeliste_1963_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1963_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1963_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1963_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1963.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1963.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1963.csv")
#write.csv(pers_unimode, file = "pers_unimode_1963.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1963_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1963.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1963.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1963.csv")
#write.csv(firm_unimode, file = "firm_unimode_1963.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1964 

edgeliste_1964_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1964_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1964_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1964_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1964.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1964.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1964.csv")
#write.csv(pers_unimode, file = "pers_unimode_1964.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1964_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1964.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1964.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1964.csv")
#write.csv(firm_unimode, file = "firm_unimode_1964.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1965 

edgeliste_1965_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1965_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1965_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1965_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1965.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1965.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1965.csv")
#write.csv(pers_unimode, file = "pers_unimode_1965.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1965_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1965.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1965.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1965.csv")
#write.csv(firm_unimode, file = "firm_unimode_1965.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1966

edgeliste_1966_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1966_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1966_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1966_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1966.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1966.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1966.csv")
#write.csv(pers_unimode, file = "pers_unimode_1966.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1966_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1966.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1966.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1966.csv")
#write.csv(firm_unimode, file = "firm_unimode_1966.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1969

edgeliste_1969_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1969_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1969_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1969_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1969.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1969.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1969.csv")
#write.csv(pers_unimode, file = "pers_unimode_1969.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1969_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1969.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1969.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1969.csv")
#write.csv(firm_unimode, file = "firm_unimode_1969.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1972

edgeliste_1972_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1972_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1972_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1972_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1972.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1972.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1972.csv")
#write.csv(pers_unimode, file = "pers_unimode_1972.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1972_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1972.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1972.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1972.csv")
#write.csv(firm_unimode, file = "firm_unimode_1972.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1975

edgeliste_1975_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1975_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1975_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1975_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1975.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1975.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1975.csv")
#write.csv(pers_unimode, file = "pers_unimode_1975.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1975_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1975.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1975.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1975.csv")
#write.csv(firm_unimode, file = "firm_unimode_1975.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1979

edgeliste_1979_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1979_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1979_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1979_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1979.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1979.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1979.csv")
#write.csv(pers_unimode, file = "pers_unimode_1979.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1979_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1979.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1979.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1979.csv")
#write.csv(firm_unimode, file = "firm_unimode_1979.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1980

edgeliste_1980_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1980_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1980_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1980_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1980.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1980.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1980.csv")
#write.csv(pers_unimode, file = "pers_unimode_1980.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1980_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1980.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1980.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1980.csv")
#write.csv(firm_unimode, file = "firm_unimode_1980.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1981

edgeliste_1981_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1981_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1981_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1981_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1981.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1981.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1981.csv")
#write.csv(pers_unimode, file = "pers_unimode_1981.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1981_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1981.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1981.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1981.csv")
#write.csv(firm_unimode, file = "firm_unimode_1981.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1982

edgeliste_1982_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1982_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1982_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1982_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1982.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1982.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1982.csv")
#write.csv(pers_unimode, file = "pers_unimode_1982.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1982_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1982.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1982.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1982.csv")
#write.csv(firm_unimode, file = "firm_unimode_1982.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1983

edgeliste_1983_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1983_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1983_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1983_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1983.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1983.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1983.csv")
#write.csv(pers_unimode, file = "pers_unimode_1983.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1983_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1983.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1983.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1983.csv")
#write.csv(firm_unimode, file = "firm_unimode_1983.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1984

edgeliste_1984_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1984_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1984_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1984_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1984.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1984.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1984.csv")
#write.csv(pers_unimode, file = "pers_unimode_1984.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1984_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1984.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1984.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1984.csv")
#write.csv(firm_unimode, file = "firm_unimode_1984.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1985

edgeliste_1985_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1985_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1985_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1985_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1985.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1985.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1985.csv")
#write.csv(pers_unimode, file = "pers_unimode_1985.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1985_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1985.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1985.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1985.csv")
#write.csv(firm_unimode, file = "firm_unimode_1985.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1986

edgeliste_1986_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1986_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1986_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1986_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1986.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1986.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1986.csv")
#write.csv(pers_unimode, file = "pers_unimode_1986.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1986_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1986.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1986.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1986.csv")
#write.csv(firm_unimode, file = "firm_unimode_1986.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1987

edgeliste_1987_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1987_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1987_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1987_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1987.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1987.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1987.csv")
#write.csv(pers_unimode, file = "pers_unimode_1987.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1987_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1987.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1987.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1987.csv")
#write.csv(firm_unimode, file = "firm_unimode_1987.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1988 

edgeliste_1988_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1988_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1988_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1988_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1988.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1988.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1988.csv")
#write.csv(pers_unimode, file = "pers_unimode_1988.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1988_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1988.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1988.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1988.csv")
#write.csv(firm_unimode, file = "firm_unimode_1988.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1989

edgeliste_1989_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1989_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1989_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1989_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1989.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1989.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1989.csv")
#write.csv(pers_unimode, file = "pers_unimode_1989.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1989_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1989.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1989.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1989.csv")
#write.csv(firm_unimode, file = "firm_unimode_1989.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1990

edgeliste_1990_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1990_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1990_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1990_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1990.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1990.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1990.csv")
#write.csv(pers_unimode, file = "pers_unimode_1990.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1990_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1990.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1990.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1990.csv")
#write.csv(firm_unimode, file = "firm_unimode_1990.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1991

edgeliste_1991_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1991_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1991_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1991_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1991.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1991.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1991.csv")
#write.csv(pers_unimode, file = "pers_unimode_1991.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1991_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1991.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1991.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1991.csv")
#write.csv(firm_unimode, file = "firm_unimode_1991.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1992

edgeliste_1992_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1992_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1992_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1992_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1992.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1992.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1992.csv")
#write.csv(pers_unimode, file = "pers_unimode_1992.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1992_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1992.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1992.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1992.csv")
#write.csv(firm_unimode, file = "firm_unimode_1992.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1993

edgeliste_1993_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1993_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1993_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1993_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1993.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1993.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1993.csv")
#write.csv(pers_unimode, file = "pers_unimode_1993.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1993_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1993.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1993.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1993.csv")
#write.csv(firm_unimode, file = "firm_unimode_1993.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1994

edgeliste_1994_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1994_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1994_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1994_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1994.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1994.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1994.csv")
#write.csv(pers_unimode, file = "pers_unimode_1994.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1994_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1994.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1994.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1994.csv")
#write.csv(firm_unimode, file = "firm_unimode_1994.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1995

edgeliste_1995_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1995_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1995_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1995_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1995.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1995.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1995.csv")
#write.csv(pers_unimode, file = "pers_unimode_1995.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1995_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1995.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1995.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1995.csv")
#write.csv(firm_unimode, file = "firm_unimode_1995.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1996

edgeliste_1996_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1996_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1996_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1996_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1996.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1996.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1996.csv")
#write.csv(pers_unimode, file = "pers_unimode_1996.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1996_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1996.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1996.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1996.csv")
#write.csv(firm_unimode, file = "firm_unimode_1996.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1997

edgeliste_1997_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1997_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1997_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1997_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1997.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1997.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1997.csv")
#write.csv(pers_unimode, file = "pers_unimode_1997.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1997_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1997.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1997.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1997.csv")
#write.csv(firm_unimode, file = "firm_unimode_1997.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1998

edgeliste_1998_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1998_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1998_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1998_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1998.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1998.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1998.csv")
#write.csv(pers_unimode, file = "pers_unimode_1998.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1998_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1998.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1998.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1998.csv")
#write.csv(firm_unimode, file = "firm_unimode_1998.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 1999

edgeliste_1999_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_1999_bipartite.dta")
eg=graph_from_data_frame(edgeliste_1999_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_1999_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_1999.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_1999.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_1999.csv")
#write.csv(pers_unimode, file = "pers_unimode_1999.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_1999_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_1999.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_1999.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_1999.csv")
#write.csv(firm_unimode, file = "firm_unimode_1999.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 2000

edgeliste_2000_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_2000_bipartite.dta")
eg=graph_from_data_frame(edgeliste_2000_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_2000_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_2000.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_2000.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_2000.csv")
#write.csv(pers_unimode, file = "pers_unimode_2000.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_2000_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_2000.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_2000.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_2000.csv")
#write.csv(firm_unimode, file = "firm_unimode_2000.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 2001

edgeliste_2001_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_2001_bipartite.dta")
eg=graph_from_data_frame(edgeliste_2001_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_2001_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_2001.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_2001.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_2001.csv")
#write.csv(pers_unimode, file = "pers_unimode_2001.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_2001_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_2001.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_2001.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_2001.csv")
#write.csv(firm_unimode, file = "firm_unimode_2001.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 2002

edgeliste_2002_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_2002_bipartite.dta")
eg=graph_from_data_frame(edgeliste_2002_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_2002_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_2002.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_2002.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_2002.csv")
#write.csv(pers_unimode, file = "pers_unimode_2002.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_2002_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_2002.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_2002.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_2002.csv")
#write.csv(firm_unimode, file = "firm_unimode_2002.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

# 2003

edgeliste_2003_bi <- read_stata("C:/Users/schmutzy/Desktop/Projekt Nationalräte/02_Processed_data/10_Directors_1934_2003/10_Network_Bipartite/edgeliste_2003_bipartite.dta")
eg=graph_from_data_frame(edgeliste_2003_bi, directed=FALSE)
eg
el <- as_edgelist(eg)

#Transforming bipartite network into unimodal
devtools::install_github("mdlincoln/projectoR")
library(projectoR)

pers_unimode_weight <- project_table(edgeliste_2003_bi, joining_col = "N_UCID")
pers_unimode = subset(pers_unimode_weight, select = -c(weight))

write.xlsx(pers_unimode_weight, file = "pers_unimode_w_2003.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(pers_unimode, file = "pers_unimode_2003.xlsx", sheetName = "pers",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(pers_unimode_weight, file = "pers_unimode_w_2003.csv")
#write.csv(pers_unimode, file = "pers_unimode_2003.csv")

pers_eg=graph_from_data_frame(pers_unimode, directed=FALSE)
pers_eg
pers_el <- as_edgelist(pers_eg)

firm_unimode_weight <- project_table(edgeliste_2003_bi, joining_col = "N_id_sug")
firm_unimode = subset(firm_unimode_weight, select = -c(weight))

write.xlsx(firm_unimode_weight, file = "firm_unimode_w_2003.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
write.xlsx(firm_unimode, file = "firm_unimode_2003.xlsx", sheetName = "firms",
           colNames = TRUE, rowNames = TRUE, append = FALSE)
#write.csv(firm_unimode_weight, file = "firm_unimode_w_2003.csv")
#write.csv(firm_unimode, file = "firm_unimode_2003.csv")

firm_eg=graph_from_data_frame(firm_unimode, directed=FALSE)
firm_eg
firm_el <- as_edgelist(firm_eg)

################################################################################

#igraph commands

plot(eg, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1934 : Firm Networks and interlocks")
plot(firm_eg, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1934 : Firm Networks")
plot(pers_eg, edge.arrow.size=0.05, vertex.size = 0.1, vertex.label="", xlab = "1934 : Interlocks")

#Bipartite network info
diameter(eg)
is_simple(eg)
is_connected(eg)
degree_distribution(eg, cumulative = TRUE)
plot(degree_distribution(eg, cumulative = TRUE), col="orange")
hist(degree(eg))
abline(v=mean(degree(eg)), col="red")
mean(degree(eg))

#Firm network info
degree(firm_eg)
gsize(firm_eg) # nb of edges
gorder(firm_eg) # how many nodes ?
diameter(firm_eg)
is_simple(firm_eg)
is_connected(firm_eg)
degree_distribution(firm_eg, cumulative = TRUE)
plot(degree_distribution(firm_eg, cumulative = TRUE), col="orange")
hist(degree(firm_eg))
abline(v=mean(degree(firm_eg)), col="red")
mean(degree(firm_eg))

#Pers network info
degree(pers_eg)
gsize(pers_eg) # nb of edges
gorder(pers_eg) # how many nodes ?
diameter(pers_eg)
is_simple(pers_eg)
is_connected(pers_eg)
degree_distribution(pers_eg, cumulative = TRUE)
plot(degree_distribution(pers_eg, cumulative = TRUE), col="orange")
hist(degree(pers_eg))
abline(v=mean(degree(pers_eg)), col="red")
mean(degree(pers_eg))

################################################################################
#Code copied from UPF crash course : sort and update!

#--------------------- Degree and degree distribution ------------------------------

# if we have a matrice: rowSums(A)
# if we have a matrice: colSums(A)

degree(eg)
gsize(eg) # nb of edges
gorder(eg) # how many nodes ?
degree(firm_eg, normalized = TRUE) # the result is divided by n-1, where n is the number of vertices in the graph.
sum(degree(eg, normalized = TRUE))

degree_distribution(eg)
sum(degree_distribution(eg))

deg <- degree(firm_eg, mode="all")
deg.dist <- degree_distribution(firm_eg, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

farthest_vertices(firm_eg)
C6180_C16860 <- get.shortest.paths(firm_eg, V(firm_eg)[name=="d"],
                                 V(firm_eg)[name=="c"],
                                 mode="out", output="both")