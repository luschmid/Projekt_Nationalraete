library(tidyverse)
df_vr <- tibble(year_vr=c(1934, 1943, 1960, seq(1962,1966,1), 1969, 1972, 1975, 
                          seq(1979,2017,1)))

df_el <- tibble(year_el=seq(1931, 2015,4))
df_el$nobs <- NA

for (i in 1:length(df_el$year_el)){
df_el$nobs[i] <- sum(df_vr$year_vr %in%   seq(df_el$year_el[i]+1,df_el$year_el[i]+4,1))
}

df_vr$year_el <- NA

for (i in 1:length(df_vr$year_vr)){
df_vr$year_el[i] <- df_el$year_el[df_el$year_el %in%  seq(df_vr$year_vr[i]-4,df_vr$year_vr[i]-1,1)]
}

df_final <- df_vr %>%
  left_join(df_el, by= c("year_el"))
print(df_final,n=60)
