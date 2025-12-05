
CreateDFAllBuergerorte <- function(data, data_buergerorte) {
  out_data <- data.frame()
  for (i in 1:dim(data)[1]) {
    df_buergeorte <- data_buergerorte %>% filter(id_0 == data$id_0[i])
    out_data <- bind_rows(out_data, data[i, ] %>%
                            left_join(df_buergeorte, by = c("id_0")))
  }
  return(out_data)
}


LookupDoubleEntries <- function(data1, data2) {
  doubles <- 0
  for (i in 1:dim(data1)[1]) {
    if (max(data2$id_0 == data1[i, ]$id_0) == 1) {
      print(paste0("double entry: ", data1[i, ]$id_0))
      doubles <- doubles + 1
    } else {}
  }
  print(paste0("No. double entries: ", doubles))
}


WriteAllVarsToWideFormat <- function(data, idvar, targetvars) {
  out <- data[1, ]
  for (k in 1:length(targetvars)) {
    for (i in 2:dim(data)[1]) {
      if (data[[idvar]][i] != data[[idvar]][i - 1]) {
        out <- bind_rows(out, data[i, ])
      } else {
        out[[targetvars[k]]][dim(out)[1]] <- paste0(out[[targetvars[k]]][dim(out)[1]], ", ", data[[targetvars[k]]][i])
      }
    }
  }
  return(out)
}



GetResultBisnode <- function(data_target, data_input,
                             max_distance_input = 0.2,
                             return_option = T,
                             append_option = T,
                             progress = F,
                             path = "/02_Processed_data/12_Record_linkage/01_Bisnode/", 
                             filename = "_") {
  
  # Note: max_distance_input[1] is for surname and max_distance_input[2] is for firstname
  
  
  if (length(max_distance_input) == 1) {
    max_distance_input <- rep(max_distance_input, 2)
  }
  out <- data.frame()
  no_results <- 0
  for (i in 1:dim(data_input)[1]) {
    if (progress == T) {
      print(paste0("Round: ", i))
    }
    bigger_result <- data_target[agrepl(data_input$name[i], data_target$name_bis, max.distance = max_distance_input[1]) &
                                   agrepl(data_input$firstname[i], data_target$firstname_bis, max.distance = max_distance_input[2]), ]
    if (dim(bigger_result)[1] > 0) {
      no_results <- no_results + 1
      bigger_result <- bigger_result %>%
        mutate(
          Nachname_Pol = data_input$name[i],
          Vorname_Pol = data_input$firstname[i],
          E_CNTR_w_Pol = data_input$E_CNTR_w[i],
          N_CNTR_w_Pol = data_input$N_CNTR_w[i],
          E_CNTR_b_Pol = data_input$E_CNTR_b[i],
          N_CNTR_b_Pol = data_input$N_CNTR_b[i],
          ID_Pol = data_input$id_0[i],
          W_Distanz = sqrt((as.numeric(as.character(data_input$E_CNTR_w[i])) - e_cntr_w_bis)^2 + (as.numeric(data_input$N_CNTR_w[i]) - n_cntr_w_bis)^2) / 1000,
          B_Distanz = sqrt((as.numeric(as.character(data_input$E_CNTR_b[i])) - e_cntr_b_bis)^2 + (as.numeric(data_input$N_CNTR_b[i]) - n_cntr_b_bis)^2) / 1000,
          Beruf_Pol = data_input$job[i]
        ) %>%
        rename(
          Vorname_VR = firstname_bis,
          Nachname_VR = name_bis,
          ID_VR = id_bis
        ) %>%
        select(
          ID_Pol, ID_VR,
          Vorname_Pol, Nachname_Pol,
          Vorname_VR, Nachname_VR,
          W_Distanz, B_Distanz,
          Beruf_Pol
        ) %>%
        mutate(
          max_distance_input_name = max_distance_input[1],
          max_distance_input_firstname = max_distance_input[2]
        ) %>%
        arrange(ID_Pol, Vorname_VR, Nachname_VR)
      
      # print(bigger_result)
      
      if (return_option == F) {
        append_option <- ifelse(no_results > 1, TRUE, append_option)
        data.table::fwrite(
          x = bigger_result,
          file = paste0(path, filename, ".csv"),
          sep = ";", append = append_option
        )
      }
      
      
      if (return_option == T) {
        out <- rbind(out, bigger_result)
      }
    }
  }
  
  if (return_option == T) {
    return(out)
  }
}





RL_StdSearch <- function(data_target, data_input, max_distances = list(c(0, 0.1), c(0.005, 0.001)),
                         filename, path, append_choice = F, progress = F) {
  
  # This function implements a simple name search (for firstname and lastname)
  # and saves the output in a specified csv file (filename)
  # append_choice should be FALSE when we start a new iteration but TRUE when we
  # add new observations to an existing batch (due to error in loop)
  
  for (i in 1:length(max_distances[[1]])) {
    for (j in 1:length(max_distances[[2]])) {
      if (append_choice == F) {
        append_choice <- ifelse(i == 1 & j == 1, F, T)
      }
      # Note: replace file only for first search input and not if you restart search
      # (then append_choice=T)
      GetResultBisnode(
        data_target = data_target,
        data_input = data_input,
        filename = filename,
        append_option = append_choice,
        path = path,
        progress = progress,
        max_distance_input = c(max_distances[[1]][i], max_distances[[2]][j]),
        return_option = F
      )
    }
  }
}




KeepOnlyMatches <- function(data) {
  return(
    data %>%
      group_by(Cluster_ID) %>%
      mutate(nobs = n()) %>%
      ungroup() %>%
      filter(nobs > 1 & is.na(Link_Score) == F)
  )
}



LookUpConnection <- function(data, data2 = "", entry, type = "pair", idvars = c("id_0", "id_1"),
                             idvars_obs = "") {
  # This function looks entries from the dataset "entry" in the target datset
  # "data". The type "pair" looks up specific entry combinations, the type
  # "observations" looks up all observations related to both ids (id_0 and id_1)
  # in the entry dataset.
  
  if (type == "pair") {
    out <- data %>%
      left_join(entry %>% distinct(id_0, id_1, .keep_all = T) %>%
                  mutate(indi = 1), by = c(idvars)) %>%
      filter(indi == 1) %>%
      distinct(id_0, id_1, .keep_all = T)
  }
  if (type == "observations") {
    entry_id_1 <- unique(entry[[idvars[2]]])
    entry_id_0 <- unique(entry[[idvars[1]]])
    out <- data %>%
      left_join(entry %>% distinct(id_0, id_1, .keep_all = T),
                by = c(idvars)
      ) %>%
      filter(id_0 %in% entry_id_0 | id_1 %in% entry_id_1) %>%
      distinct(id_0, id_1, .keep_all = T)
  }
  
  if (type == "not_matched") {
    entry_id_1 <- unique(entry[[idvars[2]]])
    entry_id_0 <- unique(entry[[idvars[1]]])
    
    out1 <- entry %>%
      left_join(data, by = c("id_0")) %>%
      distinct(id_0, .keep_all = T) %>%
      select(id_0, ends_with("_0"))
    
    out2 <- entry %>%
      left_join(data, by = c("id_1")) %>%
      distinct(id_1, .keep_all = T) %>%
      select(id_1, ends_with("_1"))
    
    out <- entry %>%
      left_join(out1, by = c("id_0")) %>%
      left_join(out2, by = c("id_1")) %>%
      arrange(id_0, id_1) %>%
      select(
        id_0, id_1, name_1, name_0, firstname_1, firstname_0, sex_1, sex_0,
        e_cntr_w_1, e_cntr_w_0, n_cntr_w_1, n_cntr_w_0,
        birthyear_1, birthyear_0, e_cntr_b_1, e_cntr_b_0, n_cntr_b_1,
        tset, vset, gt_data, rl_data, true_positives, true_negatives,
        false_positives, false_negatives
      ) %>%
      mutate(Link_Score = NA, Cluster_ID = NA)
  }
  return(out %>% select(id_0, id_1, everything()))
}