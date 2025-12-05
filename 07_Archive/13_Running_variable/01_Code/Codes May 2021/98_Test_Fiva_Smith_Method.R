party <- "h"
data_input <- data_input_test 

if (method=="dHondt"){
  df_mult <- data.frame(mult=c(1:no_seats),
                        seat=c(1:no_seats))
}
if (method=="SainteLague"){
  df_mult <- data.frame(mult=c(1.4,seq(3,no_seats*2,2)),
                        seat=c(1:no_seats))
}

data_input <- data_input %>% mutate(seat=as.numeric(seat),
                                    No=as.numeric(No))
data_party_j <- data_input %>%
  filter(party == "h") %>%
  mutate(seat_merge= as.numeric(seat)) %>%  
  left_join(df_mult, by = c("seat"))

data_others <- data_input %>%
  mutate(No_rank = dense_rank(-as.numeric(No))) %>%
  filter(No_rank <= no_seats) %>%
  mutate(seat_merge = no_seats - No_rank + 1) %>%
  filter(party != "h") %>%
  left_join(df_mult, by = c("seat")) %>%
  select(party,No,mult)

data_party_j_long <- AppendDataFrames(data_input=data_party_j,
                 parties=levels(as.factor(data_others$party)))

data_out <- data_party_j_long  %>%
  left_join(data_others, by = c("party")) %>%
  mutate(mult=pmin(mult.x , mult.y)) %>%
  mutate(margin_all = (No.x - No.y) * mult) %>%
  filter(seat==3)
  group_by(seat) %>%
  summarize(margin_min=min(margin_all),
            margin_max=max(margin_all))%>%
  mutate(party="h",
         margin=ifelse(margin_min<0,margin_max,margin_min))%>%
  select(party,seat, margin) 
  

if (margin_type="add"){
  
  data_out <- data_out  %>%
    left_join(data_others, by = c("seat_merge")) %>%
    mutate(mult=mult.x) %>%
    mutate(margin = (as.numeric(No.x) - as.numeric(No.y)) * mult) %>%
    select(party, seat, margin)
}


