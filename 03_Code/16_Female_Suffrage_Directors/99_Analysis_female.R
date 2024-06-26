
save.image(file = "./02_Processed_data/16_Female_Suffrage_Directors/Female_Suffrage_Directors.RData")

  
#------------------------
# A) Descriptive evidence
#------------------------

# (i) All Switzerland VR


sugarcube_coll_by_year_agg <- sugarcube %>%
  dplyr::group_by(year) %>%
  dplyr::summarize(mean_female_vr=1-mean(male,na.rm=T))

ggplot(sugarcube_coll_by_year_agg,aes(x=year,y=mean_female_vr)) +
  geom_point()+
  geom_line()+
  theme_bw(base_size=24) +
  xlab("Year") + ylab("Share of women among directors")

ggsave("C:/Schmidlu/Dropbox/gender_over_time.pdf",width=10,height=7)

# (ii) VR by canton

sugarcube_coll_by_year_end <- sugarcube_all %>%
  filter(year==2003 ) %>%
  mutate(xpos=year,
         ypos=mean_female_vr,
         labeler=paste0(ctn)) #,"(nobs=",nobs,")"

sugarcube_coll_by_year %>% filter(mean_female_vr>0.3)
sugarcube %>% filter(year==1933 & ctn=="ai")

library(RColorBrewer)


nb.cols <- 26
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

ggplot(sugarcube_all %>% filter(ctn!="ai"),aes(x=year,
                                  y=mean_female_vr,
                                  color=ctn)) +
  geom_point()+
  geom_line() +
  geom_text_repel(data=sugarcube_coll_by_year_end,
                  aes(x=xpos,y=ypos,label=labeler), 
                  max.overlaps = 200) +
  theme_bw(base_size=24) +
  scale_color_manual(values = mycolors) +
  theme(legend.position="none")+
  xlab("Year") + ylab("Share of women among directors")
  
ggsave("C:/Schmidlu/Dropbox/gender_over_time_by_canton.pdf",width=10,height=7)


# (iii) VR by support for female suffrage

sugarcube_all %>% 
  filter(year==2001) %>%
  arrange(yes_percent_female_1971_q) %>%
  select(yes_percent_female_1971_q,ctn,yes_percent_female_1971) %>%
  print(n=25)

sugarcube_coll_by_quintile <- sugarcube_all %>%
  dplyr::group_by(yes_percent_female_1971_q,year) %>%
  dplyr::summarize(mean_female_vr=mean(mean_female_vr,na.rm=T)) %>%
  as_tibble()

sugarcube_end <- sugarcube_coll_by_quintile %>%
  filter(year==2003 & !is.na(yes_percent_female_1971_q)) %>%
  mutate(xpos=year,
         ypos=mean_female_vr,
         labeler=paste0("Quintile ", yes_percent_female_1971_q)) #,"(nobs=",nobs,")"


mycolors <- brewer.pal(7, "Greens")[3:7]

ggplot(sugarcube_coll_by_quintile %>% filter(!is.na(yes_percent_female_1971_q)) ,
       aes(x=year,y=mean_female_vr,color=factor(yes_percent_female_1971_q))) +
  geom_point()+
  geom_line() +
  geom_text_repel(data=sugarcube_end,
                   aes(x=xpos,y=ypos,label=labeler), 
                   max.overlaps = 200) +
  theme_bw(base_size=24) +
  scale_color_manual(values = mycolors) +
  theme(legend.position="none")+
  xlab("Year") + ylab("Share of women among directors")
  

ggsave("C:/Schmidlu/Dropbox/gender_over_time_by_canton.pdf",width=10,height=7)




# (iii) NR by canton

ct_elections_coll_by_year_end <- ct_elections_coll_by_year %>%
  filter(year==LastLegislation) %>%
  mutate(xpos=year+1.2,
         ypos=mean_female_nr,
         labeler=paste0(canton))

ggplot(ct_elections_coll_by_year,
       aes(x=year,y=mean_female_nr,color=canton)) +
  geom_point()+
  geom_line() +
  geom_text(data=ct_elections_coll_by_year_end,
            aes(x=xpos,y=ypos,label=labeler)) +
  theme_bw(base_size=24)+
  scale_color_manual(values = mycolors) +
  theme(legend.position="none")+
  xlab("Year") + ylab("Share of women among candidates")

ggsave("C:/Schmidlu/Dropbox/gender_over_time_nr.pdf",width=10,height=7)



ct_elections_coll_by_year_end <- ct_elections_coll_by_year  %>%
  group_by(canton) %>%
  mutate(year_rel_max=max(year_treated,na.rm=T)) %>%
  ungroup() %>%
  mutate(xpos=year_treated+1.2,
         ypos=mean_female_nr,
         labeler=paste0(canton," (",year(suffc_date),")")) %>%
  filter(year_treated==year_rel_max)

ggplot(ct_elections_coll_by_year,
       aes(x=year_treated,
           y=mean_female_nr,
           color=canton)) +
  geom_point()+
  geom_line() +
  geom_text_repel(data=ct_elections_coll_by_year_end,aes(x=xpos,y=ypos,label=labeler)) +
  theme_bw()+
  theme(legend.position="none")+
  xlab("Year") + ylab("Share of women among candidates")


#----------------
# D) Regressions
#----------------


# (i) Female share in VR

xnam <- paste0("ctn_", ct_levels)
xnam2 <- paste0("ctn_", ct_levels,"_sq")
covs <- c("emplsec1","emplsec2","relgcath","relgprot","relgothr","infmrate")

(fmla_linear <- as.formula(paste("mean_female_vr ~ female_suffrage_ct +", paste(covs, collapse= "+"),"+",
                                 paste(xnam, collapse= "+"),"| ctn + year ")))
(fmla_squared <- as.formula(paste("mean_female_vr ~ female_suffrage_ct +", paste(covs, collapse= "+"),"+",
                                  paste(xnam, collapse= "+"),"+",paste(xnam2, collapse= "+"),"| ctn + year ")))

res1 = feols(mean_female_vr ~ female_suffrage_ct   | ctn + year, sugarcube_all)
res2 = feols(fmla_linear , sugarcube_all )
res3 = feols(fmla_squared , sugarcube_all %>% filter(ctn_year!="ai_1933"))

etable(list(res1, res2,res3), keep="female_suffrage_ct", tex=FALSE) # file ='tt.txt'

etable(list(res1, res2,res3), keep="female_suffrage_ct",
       tex=T , 
       #file =paste0(out_path,'/female_vr.tex'),
       replace = T, 
       style.tex = style.tex("aer")) # 


# drop units

sugarcube_ctn <- sugarcube_all %>% tabyl(ctn) %>% filter(ctn!="ai") %>% select(ctn) %>% pull() 

out <- data.frame()

for (ct in sugarcube_ctn){
res3_drop = feols(fmla_squared , sugarcube_all %>% filter(ctn_year!="ai_1933" & ctn!=ct))
out <- bind_rows(out,res3_drop$coeftable[1,]%>% mutate(ctn=ct)) 
}

ggplot(data=out,aes(x=factor(ctn),y=Estimate)) +
  geom_errorbar(aes(ymin = Estimate -1.96*`Std. Error`, 
                    ymax = Estimate +1.96*`Std. Error`), width = 0.2) +
  geom_point()+
  theme_bw(base_size=24)+
  geom_hline(yintercept =res3$coeftable[1,]$Estimate, linetype="dashed") +
  xlab("Canton dropped") + ylab("Estimate") +
  scale_y_continuous(limits=c(0,0.03))




# weighted regression

res1 = feols(mean_female_vr ~ female_suffrage_ct   | ctn + year, sugarcube_all, weights=~1/populcen)
res2 = feols(fmla_linear , sugarcube_all , weights=~1/populcen)
res3 = feols(fmla_squared , sugarcube_all %>% filter(ctn_year!="ai_1933"), weights=~1/populcen)

etable(list(res1, res2,res3), keep="female_suffrage_ct", tex=FALSE) # file ='tt.txt'

etable(list(res1, res2,res3), keep="female_suffrage_ct",
       tex=T , 
       #file =paste0(out_path,'/female_vr.tex'),
       replace = T, 
       style.tex = style.tex("aer")) # 




# dynamic effect

est_did = feols(mean_female_vr ~ female_suffrage_ct + i(as.factor(time_to_treatment), as.factor(female_suffrage_ct), -1) | ctn + year,
                data=sugarcube_all %>% filter(year>1960))

summary(est_did,n=100)
iplot(est_did)

sugarcube_all %>% tabyl(year)
sugarcube_all %>% tabyl(time_to_treatment)


# Sun and Abraham 2020

res_twfe = feols(mean_female_vr ~ female_suffrage_ct + i(time_to_treatment, ref = c(-1, -1000)) | ctn + year, sugarcube_all)
res_sa20 = feols(mean_female_vr ~ female_suffrage_ct + 
                   sunab(suff_year, year) | ctn + year, sugarcube_all)

# Plot the two TWFE results
iplot(list(res_twfe, res_sa20), sep = 0.5)

# Add the true results

legend("topleft", col = c(1, 4, 2), pch = c(20, 15, 17), 
       legend = c("TWFE", "Truth", "Sun & Abraham (2020)"))
       

# (v) Regression on female share in NR

xnam <- paste0("canton_", ct_levels[!ct_levels%in%c("ai","ar","gr")])
xnam2 <- paste0("canton_", ct_levels[!ct_levels%in%c("ai","ar","gr")],"_sq")

(fmla_linear <- as.formula(paste("mean_female_nr ~ female_suffrage_ct +", 
                                 paste(xnam, collapse= "+"),"| canton + year  ")))
(fmla_squared <- as.formula(paste("mean_female_nr ~ female_suffrage_ct +", 
                                  paste(xnam, collapse= "+"),"+",paste(xnam2, collapse= "+"),"| canton + year ")))


res1 = feols(mean_female_nr ~ female_suffrage_ct   | canton + year, ct_elections_coll_by_year)
res2 = feols(fmla_linear , ct_elections_coll_by_year)


# (vi) Robustness on female share in VR

xnam <- paste0("ctn_", ct_levels)
xnam2 <- paste0("ctn_", ct_levels,"_sq")

(fmla_linear <- as.formula(paste("mean_female_vr ~ female_suffrage_ct +", 
                                 paste(xnam, collapse= "+"),"| ctn + year ")))
(fmla_squared <- as.formula(paste("mean_female_vr ~ female_suffrage_ct +", 
                                  paste(xnam, collapse= "+"),"+",paste(xnam2, collapse= "+"),"| ctn + year ")))
(fmla_linear_1971 <- as.formula(paste("mean_female_vr ~ dummy_1971 +", 
                                 paste(xnam, collapse= "+"),"| ctn ")))
  
res1_exclude_small_ct = feols(fmla_linear , sugarcube_coll_by_year %>% filter(!ctn %in% c("ar","ai","nw","ow","ur")))
res1_exclude_big_ct = feols(fmla_linear , sugarcube_coll_by_year %>% filter(!ctn %in% c("ge","vd","be","zh")))
res1_1971 = feols(fmla_linear_1971 , sugarcube_coll_by_year)

summary(res1_1971,n=85)


# (vii) Robustness on female share in VR

xnam <- paste0("ctn_", ct_levels[ct_levels%in%c(c("ne","vd","bs","bl","ge"))])
(fmla_linear <- as.formula(paste("mean_female_vr ~ female_suffrage_ct +", 
                                 paste(xnam, collapse= "+"),"| ctn + year ")))

res1_exclude_small_ct = feols(fmla_linear , sugarcube_coll_by_year %>% filter(early_intro==1))
res1_exclude_small_ct = feols(fmla_linear , sugarcube_coll_by_year %>% filter(late_intro==1))


