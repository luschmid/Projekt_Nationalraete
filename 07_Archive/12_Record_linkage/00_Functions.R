ReadXlsxFiles <- function(x,path_rawdata){
  df <- read_xlsx( paste(path_rawdata, x, sep = "/")) %>%
    mutate(sourcefile=x)
  return(df)
}

ReadCsvFiles <- function(x,path_rawdata){
  df <- read_delim( paste(path_rawdata, x, sep = "/"),delim=";") %>%
    mutate(sourcefile=x)
  return(df)
}

ReplacePatternLoop <- function(data,variable,pattern,pattern_replace){
  for (i in 1:dim(data)[1]){
    data[[variable]][i] <- gsub(pattern, pattern_replace,data[[variable]][i])  
  }
  return(data)
}

CreateDFAllBuergerorte <- function(data,data_buergerorte){
  out_data <- data.frame()
  for (i in 1:dim(data)[1]){
    df_buergeorte <- data_buergerorte %>% filter(id_0==data$id_0[i])
    out_data <- bind_rows(out_data,data[i,]%>%  
                            left_join(df_buergeorte, by=c("id_0")))
  }
  return(out_data)
}


LookupDoubleEntries <- function(data1,data2){
  doubles=0
  for (i in 1:dim(data1)[1]){
    if (max(data2$id_0==data1[i,]$id_0)==1){
      print(paste0("double entry: ",data1[i,]$id_0 ))
      doubles=doubles+1
    } else{}
  }
  print(paste0("No. double entries: ",doubles ))
}

GetResultBisnode <- function(data_target, data_input,
                             max_distance_input=0.2,
                             return_option=T,
                             append_option=T,
                             progress=F,
                             path="/02_Processed_data/12_Record_linkage/01_Bisnode/",filename="_"){
  
  # Note: max_distance_input[1] is for surname and max_distance_input[2] is for firstname
  
  
  if (length(max_distance_input)==1){max_distance_input=rep(max_distance_input,2)}
  out <- data.frame()
  no_results <- 0
  for (i in 1:dim(data_input)[1]){
    if(progress==T){print(paste0("Round: ",i))}
    bigger_result <- data_target[agrepl(data_input$name[i], data_target$name_bis, max.distance=max_distance_input[1]) & 
                                   agrepl(data_input$firstname[i], data_target$firstname_bis, max.distance=max_distance_input[2]),]
    if (dim(bigger_result)[1]>0){
      no_results <- no_results + 1
      bigger_result <- bigger_result %>% mutate(Nachname_Pol=data_input$name[i],
                                                Vorname_Pol=data_input$firstname[i],
                                                E_CNTR_w_Pol=data_input$E_CNTR_w[i],
                                                N_CNTR_w_Pol=data_input$N_CNTR_w[i],
                                                E_CNTR_b_Pol=data_input$E_CNTR_b[i],
                                                N_CNTR_b_Pol=data_input$N_CNTR_b[i],
                                                ID_Pol=data_input$id_0[i],
                                                W_Distanz=sqrt((as.numeric(as.character(data_input$E_CNTR_w[i]))-e_cntr_w_bis)^2+(as.numeric(data_input$N_CNTR_w[i])-n_cntr_w_bis)^2)/1000,
                                                B_Distanz=sqrt((as.numeric(as.character(data_input$E_CNTR_b[i]))-e_cntr_b_bis)^2+(as.numeric(data_input$N_CNTR_b[i])-n_cntr_b_bis)^2)/1000,
                                                Beruf_Pol=data_input$job[i]) %>% 
        rename(Vorname_VR=firstname_bis ,
               Nachname_VR=name_bis,
               ID_VR=id_bis)%>%
        select(ID_Pol,ID_VR,
               Vorname_Pol,Nachname_Pol,
               Vorname_VR,Nachname_VR,
               W_Distanz,B_Distanz,
               Beruf_Pol)%>%
        mutate(max_distance_input_name=max_distance_input[1],
               max_distance_input_firstname=max_distance_input[2])%>%
        arrange(ID_Pol,Vorname_VR,Nachname_VR)
      
      #print(bigger_result)
      
      if(return_option==F){
        append_option <- ifelse(no_results>1,TRUE,append_option) 
        data.table::fwrite(x=bigger_result,
                           file=paste0(path,filename,".csv"), 
                           sep=";",append=append_option)}
      
      
      if(return_option==T){out <- rbind(out,bigger_result)}
    }
  } 
  
  if(return_option==T){return(out)}
}

WriteAllVarsToWideFormat <- function(data,idvar,targetvars){
  out <- data[1,]
  for (k in 1:length(targetvars)){
    for (i in 2:dim(data)[1]){
      if (data[[idvar]][i]!= data[[idvar]][i-1]){
        out <- bind_rows(out,data[i,])
      } else {
        out[[targetvars[k]]][dim(out)[1]]  <- paste0(out[[targetvars[k]]][dim(out)[1]],", ",data[[targetvars[k]]][i])
      }
    }
  }
  return(out)
}



ReadinResults <- function(file_name,path_rel){
  # This function reads in the record linkage result file. 
  return(fread(paste0(path_rel,file_name), 
               sep = ",", encoding = "UTF-8") %>%
           select(`source file`,`Cluster ID`,`Link Score`,id,name,firstname,birthyear,sex,e_cntr_w,n_cntr_w,e_cntr_b,n_cntr_b,language_w) %>%
           rename(source_file=`source file`,
                  Cluster_ID=`Cluster ID`,
                  Link_Score=`Link Score`) %>%
           as_tibble()%>% 
           arrange(Link_Score,Cluster_ID))
}

ReadinGroundTruthFiles <- function(file_name,separator=",",path_rel){
  # This function reads in the ground truth files. 
  return(fread(paste0(path_rel,file_name), 
               sep = separator, encoding = "UTF-8") %>%
           mutate(id_0=as.character(ID),
                  id_1=as.character(personenid)) %>%
           select(id_0,id_1) %>%
           as_tibble())
}

CreateWideDataset <- function(data,vars=c("name","firstname","sex",
                                          "e_cntr_w","n_cntr_w","birthyear","e_cntr_b",
                                          "n_cntr_b","language_w","id"),
                              generate_differences=F){
  
  Cluster_ID_max <- max(data$Cluster_ID,na.rm=T)
  
  if (generate_differences==T){
    return(data %>% 
             mutate(Cluster_ID = ifelse(is.na(Cluster_ID),Cluster_ID_max+row_number(),Cluster_ID)) %>%
             pivot_wider(id_cols=c(Cluster_ID,Link_Score),
                         values_from=vars,
                         names_from=source_file) %>%
             mutate(birthyear_missing=is.na(birthyear_1),
                    birthyear_diff=birthyear_1-birthyear_0,
                    w_dist=sqrt((e_cntr_w_0-e_cntr_w_1)^2 + (n_cntr_w_0-n_cntr_w_1)^2),
                    b_dist=sqrt((e_cntr_b_0-e_cntr_b_1)^2 + (n_cntr_b_0-n_cntr_b_1)^2))%>%
             select(id_0,id_1,Cluster_ID,Link_Score, name_1 ,name_0,  
                    firstname_1, firstname_0, sex_1, sex_0,w_dist,b_dist,birthyear_diff,
                    everything())
    )} else{
      return(data %>% 
               mutate(Cluster_ID = ifelse(is.na(Cluster_ID),Cluster_ID_max+row_number(),Cluster_ID)) %>%
               pivot_wider(id_cols=c(Cluster_ID,Link_Score),
                           values_from=vars,
                           names_from=source_file)) %>%
        select(id_0,id_1,Cluster_ID,Link_Score, name_1 ,name_0,  
               firstname_1, firstname_0, sex_1, sex_0,
               everything())                  
    }
}

KeepOnlyMatches <- function(data){
  return(
    data %>%
      group_by(Cluster_ID) %>%
      mutate(nobs=n()) %>%
      ungroup() %>%
      filter(nobs>1 & is.na(Link_Score)==F) 
  )
}

PreparationPrecisionRecall <- function(recordlinkagedata,
                                       groundtruthdata,
                                       idvars=c("id_0","id_1"),
                                       tvvars=c("tset","vset"),
                                       source="rl",
                                       linkscore_cutoff=0.00000000000000001){
  # This function prepares the data for the calculation of precision and recall 
  # for the dataset "data" based on the dataset "groundtruth". It labels every
  # overvation as "true positive", "true negative", "false positive", or 
  # "false negative". 
  # source="rl" is for rl output by Sandro, source="std" is for our output using
  
  id_0_groundtruthdata <- groundtruthdata %>% 
    distinct(id_0) %>%
    pull()
  
  if (source=="rl"){
    
    
    df_cut <- recordlinkagedata %>% 
      filter(!is.na(id_1)) %>% # filter out only existing matches
      filter(is.element(id_0,id_0_groundtruthdata)) %>% # Note: focus only on politician ids in ground truth data
      mutate(id_1=as.character(id_1)) %>%
      group_by(id_0,id_1) %>%
      dplyr::summarize(Link_Score_max=max(Link_Score,na.rm=T),
                       nobs=n()) %>%
      ungroup() %>%
      dplyr::mutate(Link_Score_max=ifelse(Link_Score_max==-Inf,NA,Link_Score_max))%>%
      group_by(id_0) %>%
      dplyr::mutate(nobs_all=sum(!is.na(Link_Score_max)))
      
    
    df <- recordlinkagedata %>% 
      filter(!is.na(id_1)) %>% # filter out only existing matches
      filter(is.element(id_0,id_0_groundtruthdata)) %>% # Note: focus only on politician ids in ground truth data
      mutate(rl_data=1) %>%
      mutate(id_1=as.character(id_1)) %>%
      left_join(df_cut,by=c("id_0","id_1")) %>%
      filter(Link_Score==Link_Score_max | (is.na(Link_Score_max) & nobs_all==0)) %>% # keep only maximum of Link value per dyad (id_0-id_1) except for those only with NA match
      distinct(id_0,id_1,.keep_all = T) %>% # keep only one observations for those observations where two have the same maximum Link_Score     
      full_join(groundtruthdata %>% mutate(gt_data=1,id_1=as.character(id_1)),
                by=c(idvars)) %>%
      mutate(gt_data=ifelse(is.na(gt_data),0,1),
             rl_data=ifelse(is.na(rl_data),0,1),
             Link_Score=ifelse(is.na(Link_Score),0,Link_Score)) %>% # Note: replace Link_Score=0 for all those (i) alone in cluster or (ii) only in gt data
      arrange(id_1,id_0) %>%
      #filter(tset==1 | (rl_data==1 & is.na(vset))) %>% # Note: focus only on test set observations or on observations from rl output not in validation set
      mutate(true_positives=ifelse(rl_data==1 & gt_data==1 & Link_Score>=linkscore_cutoff,1,0),
             true_negatives=ifelse(rl_data==0 & gt_data==1 & is.na(id_1)
                                   #|(rl_data==1 & gt_data==1 & Link_Score<linkscore_cutoff),1,0), # changed on December 8, 2021
                                   |(rl_data==1 & gt_data==0 & Link_Score<linkscore_cutoff),1,0),  
             false_positives=ifelse(rl_data==1 & gt_data==0 & 
                                      Link_Score>=linkscore_cutoff & 
                                      !is.na(id_0) & !is.na(id_1)  # Note: Focus only connections (all other with NaN in id_1 or id_0 are not connections)
                                    ,1,0),
             false_negatives=ifelse((rl_data==0 & gt_data==1 & !is.na(id_0) &  !is.na(id_1))  |
                                      (rl_data==1 & gt_data==1 & !is.na(id_0) &  !is.na(id_1) 
                                    & Link_Score<linkscore_cutoff) ,1,0),
             category=case_when(
               true_positives == 1 ~ "true positive",
               true_negatives == 1 ~ "true negatives",
               false_positives == 1 ~ "false positive",
               false_negatives == 1 ~ "false negative")) %>%
      select(Cluster_ID,Link_Score,id_0,id_1,true_positives,true_negatives,
             false_positives,false_negatives,category,rl_data,gt_data,vset,tset) 
    
    
      #filter(rl_data==0 | gt_data==1 | Link_Score>=linkscore_cutoff | is.na(id_0)  ) # filter out those in rl data but not in ground truth with low link score
  }
  
  if (source=="std"){
    df <- recordlinkagedata %>% 
      filter(is.element(id_0,id_0_groundtruthdata)) %>% # Note: focus only on politician ids in ground truth data
      mutate(rl_data=1) %>%
      mutate(id_1=as.character(id_1)) %>%
      full_join(groundtruthdata %>% mutate(gt_data=1,id_1=as.character(id_1)),
                by=c(idvars)) %>%
      mutate(gt_data=ifelse(is.na(gt_data),0,1),
             rl_data=ifelse(is.na(rl_data),0,1)) %>% 
      arrange(id_1,id_0) %>% 
      mutate(true_positives=ifelse(rl_data==1 & gt_data==1 ,1,0),
             true_negatives=ifelse((rl_data==0 & gt_data==1 & is.na(id_1)) |
                                     (rl_data==1 & gt_data==1 ),1,0),
             false_positives=ifelse(rl_data==1 & gt_data==0  & !is.na(id_0) & !is.na(id_1)  # Note: Focus only on connections (all other with NaN in id_1 or id_0 are not connections)
                                    ,1,0),
             false_negatives=ifelse((rl_data==0 & gt_data==1 & !is.na(id_1))  |
                                      (rl_data==1 & gt_data==1 & !is.na(id_1)) 
                                    ,1,0),
             category=case_when(
               true_positives == 1 ~ "true positive",
               true_negatives == 1 ~ "true negatives",
               false_positives == 1 ~ "false positive",
               false_negatives == 1 ~ "false negative")) 
  }
  
  return(df)
}

GetPrecisionRecall <- function(data){
  # This function calculates precision and recall for the dataset "data" based on 
  # the dataset "groundtruth"
  
  df <- data %>% 
    summarize(true_positives=sum(true_positives,na.rm=T),
              true_negatives=sum(true_negatives,na.rm=T),
              false_positives=sum(false_positives,na.rm=T),
              false_negatives=sum(false_negatives,na.rm=T))
  
  # ############
  # df_cut <- data %>% 
  #   dplyr::group_by(id_0,id_1) %>%
  #   dplyr::summarize(nobs=sum(!is.na(id_0)),
  #                    ls_min=min(Link_Score,na.rm=T),
  #                    ls_max=max(Link_Score,na.rm=T)) %>%
  #   arrange(id_0,id_1)
  # 
  # df_dups <- data %>% left_join(df_cut, by=c("id_0","id_1"))%>%
  #   filter(nobs>1 & !is.na(category))  %>%
  #   mutate(ls_diff=ls_max-ls_min)%>%
  #   arrange(ls_diff,id_0,id_1) 
  # View(df_dups)
  # ############
  # 
  # ############
  # data %>% group_by(id_0) %>% mutate(tn_max=max(true_negatives,na.rm=T),
  #                                    tp_max=max(true_positives,na.rm=T)) %>%
  #   ungroup() %>%
  #   filter(tn_max==1&tp_max==1) %>%
  #   arrange(id_0) %>%
  #   select(id_0,id_1)
  # ############
  
  df_all <- df %>%  
    summarize(n_all=true_positives+true_negatives+false_positives+false_negatives)
  
  precision_intensive <- df$true_positives/(df$true_positives+df$false_positives)
  recall_intensive <- df$true_positives/(df$true_positives+df$false_negatives) # true positive rate (TPR)
  specificity_intensive <- df$true_negatives/(df$true_negatives+df$false_positives) # 1- false positive rate (FPR)
  
  df_true_positives <- data %>% filter(true_positives==1) %>% select(id_0) %>%
    mutate(true_pos=1)
  df_false_negatives_minustruepos <- data %>%
    filter(false_negatives==1) %>%
    left_join(df_true_positives,by=c("id_0")) %>%
    filter(is.na(true_pos))
  
  precision_extensive <- df$true_positives/(df$true_positives+df$false_positives)
  recall_extensive <- df$true_positives/(df$true_positives+dim(df_false_negatives_minustruepos)[1]) # true positive rate (TPR)
  specificity_extensive <- df$true_negatives/(df$true_negatives+df$false_positives) # 1- false positive rate (FPR)
  
  return(data.frame(precision_intensive,recall_intensive,specificity_intensive,
                    precision_extensive,recall_extensive,specificity_extensive))
}


LookUpConnection <- function(data,data2="",entry,type="pair",idvars=c("id_0","id_1"),
                             idvars_obs=""){
  # This function looks entries from the dataset "entry" in the target datset
  # "data". The type "pair" looks up specific entry combinations, the type
  # "observations" looks up all observations related to both ids (id_0 and id_1)
  # in the entry dataset. 
  
  if (type=="pair"){
    out <- data %>% left_join(entry %>% distinct(id_0,id_1,.keep_all=T)%>% 
                                mutate(indi=1),by=c(idvars))%>%
      filter(indi==1) %>%
      distinct(id_0,id_1,.keep_all=T)
  }
  if (type=="observations"){
    entry_id_1 <- unique(entry[[idvars[2]]])
    entry_id_0 <- unique(entry[[idvars[1]]])
    out <- data %>% left_join(entry %>% distinct(id_0,id_1,.keep_all=T) 
                              ,by=c(idvars))%>%
      filter(id_0  %in% entry_id_0 | id_1  %in% entry_id_1 ) %>%
      distinct(id_0,id_1,.keep_all=T)
  }
  
  if (type=="not_matched"){
    entry_id_1 <- unique(entry[[idvars[2]]])
    entry_id_0 <- unique(entry[[idvars[1]]])
    
    out1 <- entry %>% left_join(data ,by=c("id_0"))%>% 
      distinct(id_0,.keep_all=T)%>%
      select(id_0,ends_with("_0"))
    
    out2 <- entry %>% left_join(data ,by=c("id_1"))%>% 
      distinct(id_1,.keep_all=T)%>%
      select(id_1,ends_with("_1"))
    
    out <- entry %>% left_join(out1 ,by=c("id_0")) %>% 
      left_join(out2 ,by=c("id_1")) %>% 
      arrange(id_0,id_1) %>% 
      select(id_0,id_1,name_1,name_0,firstname_1,firstname_0,sex_1,sex_0,
             e_cntr_w_1,e_cntr_w_0,n_cntr_w_1,n_cntr_w_0,
             birthyear_1,birthyear_0,e_cntr_b_1,e_cntr_b_0,n_cntr_b_1,
             tset,  vset, gt_data,rl_data,true_positives,true_negatives,
             false_positives,false_negatives) %>%
      mutate(Link_Score=NA,Cluster_ID=NA)
  }
  return(out %>% select(id_0,id_1,everything())) 
}



CalculateROCCurve <- function(recordlinkagedata,
                              groundtruthdata,
                              idvars=c("id_0","id_1"),
                              tvvars=c("tset","vset"),
                              step_size=0.02,
                              return_option=F,
                              path="",
                              path_rel="",
                              file_name="",
                              source="rl"){
  # This function calculates the ROC (Receiver operating characteristic) curve
  # for a given recordlinkagedata based on a groundtruthdata
  # source="rl" is for rl output by Sandro, source="std" is for our output using
  # simple name comparisons
  #metrics_tset <- NULL
  metrics_vset <- NULL
  if (source=="rl"){
    for (i in seq(0,1,step_size)){
      dataprep <-   PreparationPrecisionRecall(recordlinkagedata=recordlinkagedata,
                                               groundtruthdata = ground_truth_sep,
                                               idvars=idvars,
                                               tvvars=tvvars,
                                               linkscore_cutoff = i,
                                               source=source) 
      
      dataprep_vs <- dataprep %>%
        filter(vset==1 | gt_data==0) %>% # comment in after Sandro has corrected assignment of test and validation set based on id_0 and not based on dyad
        group_by(id_0) %>%
        dplyr::summarize(true_positives=mean(true_positives),
                         true_negatives=mean(true_negatives), 
                         false_positives=mean(false_positives),
                         false_negatives=mean(false_negatives)) %>%
        mutate(sum_all=true_positives+true_negatives+false_positives+false_negatives)
      
      #dataprep_vs %>% filter(sum_all!=1)
      # 
      # 
      # ##############
      # 
      # dataprep_vs %>% filter(Link_Score<0.02 & !is.na(id_1)) %>% 
      #   select(id_0,id_1,Link_Score,category) 
      # 
      # ground_truth_sep %>% filter(id_0=="ZH-2011-0214")
      # dataprep_vs %>% filter(true_negatives==1)%>% 
      #   select(id_0,id_1,Link_Score,category) 
      # 
      # dataprep_vs %>% filter(id_0=="AG-1995-0110")
      # 
      # dataprep_vs %>% janitor::tabyl(category)
      # 
      # dataprep_vs %>% filter(category=="true positive") %>% 
      #   select(id_0,id_1,Link_Score,category) 
      # 
      # ground_truth_sep %>% filter(id_0=="SG-1999-0039")
      # 
      # ggplot(dataprep_vs %>% filter(!is.na(category)),aes(x=Link_Score,color=factor(category)))+
      #   geom_histogram(alpha=0.2) +
      #   scale_color_brewer(palette = "Set1")+ 
      #   facet_wrap(~category) +
      #   theme_bw()
      #   
      # dataprep_vs %>% mutate(sum=n_fn+n_tp) %>% 
      #   summarize(avg=mean(sum), max=max(sum), min=min(sum))
      # 
      #      
      # ############
      
      df_check <- dataprep %>% mutate(sum=true_positives+true_negatives+false_positives+false_negatives) 
      if(dim(df_check %>% filter(sum>1))[1]>0){
        print("More than one classification per case for the following observations:")
        print(df_check %>% select(id_0,id_1),n=1000)}
      if(dim(df_check %>% filter(sum==0))[1]>0){
        print("One case with no classification:")
        print(df_check %>% select(id_0,id_1),n=1000)}
      
      
      
      n_fn <- sum(dataprep_vs$false_negatives)
      n_fp <- sum(dataprep_vs$false_positives)
      n_tn <- sum(dataprep_vs$true_negatives)
      n_tp <- sum(dataprep_vs$true_positives)
      
      # ###############
      # 
      # dataprep_vs %>% filter(true_negatives==1) %>%
      #   summarize(sum=sum(Link_Score))
      # 
      # dataprep_vs %>% filter(true_negatives==1 & Link_Score>0) 
      # 
      # ###############
      
      e_m <- (n_fp-n_fn)/(n_fp+n_fn+n_tn+n_tp)
      var_e <- (n_fp*(1-e_m)^2+(n_tn+n_tp)*(-e_m)^2+n_fn*(-1-e_m)^2)/(n_fp+n_fn+n_tn+n_tp-1)
    
      metrics_vset <-  bind_rows(metrics_vset,
                                 GetPrecisionRecall(data=dataprep_vs) %>% 
                                   mutate(linkscore_cutoff = i,
                                          var_e=var_e,
                                          n_fn=n_fn,
                                          n_fp=n_fp,
                                          n_tn=n_tn,
                                          n_tp=n_tp,
                                          total=n_fn+n_fp+n_tn+n_tp,
                                          e_m=e_m))
    }
    ggplot(data=metrics_vset,aes(x=(1-specificity_extensive), # False positive rate
                                 y=recall_intensive  # True positive rate
    ))+
      geom_point(size=4)+
      theme_bw(base_size=48)+
      geom_text(aes(label=linkscore_cutoff,y=recall_intensive-0.02))+
      scale_x_continuous(limits=c(0,1)) +
      scale_y_continuous(limits=c(0,1)) +
      geom_abline(intercept = 0, slope = 1) +
      #scale_color_brewer(palette = "Set1") +
      theme(legend.title = element_blank()) +
      ylab("True positive rate (TP/(TP+FN))") + 
      xlab("False positive rate (FP/(FP+TN))") +
      labs(title = "ROC Curve", 
           caption = "TP=True Positives,\n FN=False Negatives,\n FP=False Positives,\n TN=True Negatives")
  
    # ggplot(data=metrics_vset,aes(x=linkscore_cutoff,
    #                              y=var_e))+
    #   geom_point(size=4)

    
    }
  
  if (source=="std"){
    # recordlinkagedata$cutoff <- paste0(recordlinkagedata$max_distance_input_name,
    #                                    "/",
    #                                    recordlinkagedata$max_distance_input_firstname )
    # levels_all <- levels(as.factor(recordlinkagedata$cutoff))
    levels_all <- seq(0,100,5)
    for (i in levels_all){
      df <- recordlinkagedata[recordlinkagedata$w_distance<=i & !is.na(w_distance),]
      dataprep <-   PreparationPrecisionRecall(recordlinkagedata=df,
                                               groundtruthdata = ground_truth_sep,
                                               idvars=idvars,
                                               tvvars=tvvars,
                                               source=source)  
      
      # metrics_tset <- bind_rows(metrics_tset,
      #                           GetPrecisionRecall(dataprep %>% 
      #                                              filter(tset==1 | gt_data==0))%>% 
      #                             mutate(linkscore_cutoff = i))
      metrics_vset <-  bind_rows(metrics_vset,
                                 GetPrecisionRecall(dataprep %>% 
                                                      filter(vset==1 | gt_data==0))%>%
                                   mutate(cutoff = i))
    }
    ggplot(data=metrics_vset,aes(x=(1-specificity_extensive), # False positive rate
                                 y=recall_intensive)) + # True positive rate
      geom_point(size=4)+
      geom_text(aes(label=cutoff,y=recall_intensive-0.02))+
      theme_bw(base_size=48)+
      scale_x_continuous(limits=c(0,1)) +
      scale_y_continuous(limits=c(0,1)) +
      geom_abline(intercept = 0, slope = 1) +
      #scale_color_brewer(palette = "Set1") +
      theme(legend.title = element_blank()) +
      ylab("True positive rate (TP/(TP+FN))") + 
      xlab("False positive rate (FP/(FP+TN))") +
      labs(title = "ROC Curve", 
           caption = "TP=True Positives,\n FN=False Negatives,\n FP=False Positives,\n TN=True Negatives")
  }
  
  # metrics_all <- bind_rows(metrics_vset %>% mutate(category="Validation set"),
  #                          metrics_tset %>% mutate(category="Test set"))
  

  if (path_rel!=""){
    ggsave(paste0(path,"/",path_rel,"/",file_name,".pdf"),width=20,height =20) 
  }
  
  if (return_option==T){return(metrics_all)}
}  


RL_StdSearch <- function(data_target,data_input,max_distances=list(c(0,0.1),c(0.005,0.001)),
                         filename,path,append_choice=F,progress=F){
  
  # This function implements a simple name search (for firstname and lastname)
  # and saves the output in a specified csv file (filename)
  # append_choice should be FALSE when we start a new iteration but TRUE when we 
  # add new observations to an existing batch (due to error in loop)
  
  for (i in 1:length(max_distances[[1]])){
    for (j in  1:length(max_distances[[2]])){
      if (append_choice==F) {append_choice <- ifelse(i==1 & j==1,F,T)}
      # Note: replace file only for first search input and not if you restart search
      # (then append_choice=T)
      GetResultBisnode(data_target=data_target,
                       data_input=data_input,
                       filename=filename,
                       append_option=append_choice,
                       path=path,
                       progress=progress,
                       max_distance_input=c(max_distances[[1]][i],max_distances[[2]][j]),
                       return_option=F)     
    }
  }
}