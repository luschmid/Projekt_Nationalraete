ls()

# Functions ----
 
Combine_all_files <- function(contents){
  out <- tibble()
  for (i in 1:length(names(contents))){
    
    out <- bind_rows(out,ExtractResults(raw_text=contents[[i]])%>%
                       mutate(type=names(contents)[i]))
  }
  return(out)
}

ExtractResults <- function(raw_text){
  
  # Step 1: Clean and split the line
  
  estimate_lines <- grep("Estimate", raw_text, value = TRUE)
  
  out <- tibble()
  
  for (i in 1: length(estimate_lines)){
    cleaned_line <- gsub("^Estimate\\s+&\\s*|\\\\\\\\", "", estimate_lines[i])
    values <- strsplit(cleaned_line, "&")[[1]]
    values <- trimws(values)
    
    # Step 2: Extract numeric part (remove \sym{})
    numeric_values <- as.numeric(gsub("\\\\sym\\{.*?\\}", "", values))
    
    # Step 3: Extract significance level
    significance <- sapply(values, function(x) {
      if (grepl("\\\\sym\\{\\*\\*\\*\\}", x)) {
        "sig1"
      } else if (grepl("\\\\sym\\{\\*\\*\\}", x)) {
        "sig5"
      } else if (grepl("\\\\sym\\{\\*\\}", x)) {
        "sig10"
      } else {
        NA
      }
    })
    
    # Step 4: Combine into a data frame
    df_est <- data.frame(
      Estimate = numeric_values,
      Significance = significance,
      stringsAsFactors = FALSE,
      Panel = i
    )
    
    out <- bind_rows(out,df_est)
  }
  return(out %>%
           mutate(horizon = case_when(
             row_number() %in% sort(as.vector(outer(c(0, 18, 36,54,72,90), 1:6, `+`))) ~ "-1",
             row_number() %in% sort(as.vector(outer(c(0, 18, 36,54,72,90), 7:12, `+`)))       ~ "0",
             row_number() %in% sort(as.vector(outer(c(0, 18, 36,54,72,90), 13:18, `+`)))     ~ "1")))

}

PlotResults <- function(df_all_files){
  
  panels <- length(levels(as.factor(df_all_files$Panel)))  
  
  for (i in 1:panels){
    for (j in -1:1){
      
      df_all <- df_all_files %>%
        filter(Panel==i & horizon == j )
      
      df_all$Significance[df_all$Significance == "sig10"] <- "Significant at 10% level"
      df_all$Significance[df_all$Significance == "sig5"]  <- "Significant at 5% level"
      df_all$Significance[df_all$Significance == "sig1"]  <- "Significant at 1% level"
      df_all$Significance[is.na(df_all$Significance)==TRUE]  <- "Not significant"
      
      # Convert back to factor and order the levels for plotting
      df_all$Significance <- factor(df_all$Significance, levels = c(
        "Significant at 10% level",
        "Significant at 5% level",
        "Significant at 1% level",
        "Not significant"
      ))
      
      custom_colors <- c(
        "Significant at 10% level" = brewer.pal(3, "Set1")[1],
        "Significant at 5% level"  = brewer.pal(3, "Set1")[2],
        "Significant at 1% level"  = brewer.pal(3, "Set1")[3],
        "Not significant"          = "black"
      )
      
      ggplot(df_all,aes(x=factor(type_nice),y=Estimate,color=factor(Significance))) + 
        geom_point() + 
        theme_bw(base_size=24) + 
        scale_color_manual(values = custom_colors) +
        ggtitle(paste0("Panel ",LETTERS[i])) + 
        labs(color = NULL, x=NULL)
      
      ggsave(paste0(ol_path,"figures/Rob_Het_","Panel",LETTERS[i],"_",as.character(j),".pdf"),
             height=10,width=15)
      
    }
  }
}

# (A) Set directories and load packages ----


setwd("C:/Schmidlu/Dropbox/Apps/Overleaf/Political_Rents/tables")
ol_path <- "C:/Schmidlu/Dropbox/Apps/Overleaf/Political_Rents/"

library(tidyverse)
library(RColorBrewer)

# (B) Read in all files and combine them to one dataframe ----

tex_files <- list.files(pattern = "^app.*\\.tex$")
tex_contents <- lapply(tex_files, readLines)
names(tex_contents) <- tools::file_path_sans_ext(tex_files)


df_final <- Combine_all_files(contents=tex_contents) %>%
  mutate(type_nice = case_when ( 
    type=="app_c1_poly_1_kern_tri" ~ "main",
    type=="app_c2_poly_1_kern_tri" ~ "c2", 
    type=="app_c1_poly_1_kern_uni" ~ "uni", 
    type=="app_c1_poly_2_kern_tri" ~ "pol2", 
    type=="app_c1_poly_1_kern_tri_trim" ~ "trim", 
    type=="app_c1_poly_1_kern_tri_left" ~ "left", 
    type=="app_c1_poly_1_kern_tri_center" ~ "center", 
    type=="app_c1_poly_1_kern_tri_right" ~ "right",
    type=="app_c1_poly_1_kern_tri_women" ~ "women",
    type=="app_c1_poly_1_kern_tri_men" ~ "men",
    type=="app_c1_poly_1_kern_tri_early" ~ "early",
    type=="app_c1_poly_1_kern_tri_late" ~ "late",
  ))

df_final$type_nice <- factor(df_final$type_nice, 
                                 levels = c("main", "c2", "uni","pol2","trim",
                                            "left","center","right","women",
                                            "men","early","late"))



# (C) Plot results and save in overleaf folder ----

PlotResults(df_final)






# old code
  
# rob_c2 <- readLines("app_c2_poly_1_kern_tri.tex")
# df_rob_c2 <- ExtractResults(raw_text=rob_c2) %>%
#  mutate(type="2_rob_c2")



