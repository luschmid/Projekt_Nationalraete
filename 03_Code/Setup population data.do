clear
set more off
version 17

*global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "C:\Current\Dropbox\Projekt Nationalräte"
global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

**************************************************
** Setup file: Construct population comparisons **
**************************************************

** Councilor data
use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep ID ID_num year elected dir_year n_all_sum_c1 n_lrg_sum_c1 n_sml_sum_c1 i_all_c1 i_lrg_c1 i_sml_c1

tsset ID_num year
gen ellegperiod = .
forv i = 1(1)4 {
replace ellegperiod = L`i'.elected if L`i'.elected != .
}
order ID_num ID year elected ellegperiod
keep if dir_year == 1
drop elected dir_year ID_num

** Director data
append using "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_nonNRs.dta"
append using "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_nonNRs.dta"
order ID PID year 
save "$path\02_Processed_data\21_Population\NR&directors.dta", replace

gen n_directors = 1
collapse (sum) n_directors, by(year)
save "$path\02_Processed_data\21_Population\n_directors.dta", replace


** Polulation data

// import population statistics
import exc using "$path\01_Raw_data\19_Population\px-x-0102030000_101_20240912-120139.xlsx",  ///
              cellrange(C3:CQ22) firstrow clear

* The full legal age was 20 before 1996. Howwever, in some cases, officials
* could lower it to 18 for specific individuals. Thus, we use 18 throughout (see
* https://hls-dhs-dss.ch/de/articles/010367/2009-11-26/; accessed September 
* 9, 2024.)

foreach v of varlist C-CQ {
   local x : variable label `v'
   rename `v' adult_pop`x'
}
gen adult = 1 if [_n] == 1
replace adult = -1 if [_n] > 1

collapse (sum) adult_pop1931-adult_pop2023, by(adult)

foreach v of varlist adult_pop1931-adult_pop2023 {
   replace `v' = `v' * adult
}
collapse (sum) adult_pop1931-adult_pop2023
gen obs = 1

reshape long adult_pop, i(obs) j(year)
drop obs

merge 1:1 year using "$path\02_Processed_data\21_Population\n_directors.dta"
keep if _merge == 3
drop _merge

gen dummy_pop = adult_pop - n_directors
keep year dummy_pop
save "$path\02_Processed_data\21_Population\pop_res.dta", replace


// create dummy population

clear all
set obs 2931172
gen year = 1934
save "$path\02_Processed_data\21_Population\dummy_pop.dta", replace

clear all
program dummy_pop
clear
set obs `2' 
gen year = `1'
append using "$path\02_Processed_data\21_Population\dummy_pop.dta"
save "$path\02_Processed_data\21_Population\dummy_pop.dta", replace
end

*pop_res.dta
*dummy_pop	1934	2931172
dummy_pop	1943	3142891
dummy_pop	1960	3763555
dummy_pop	1962	3925800
dummy_pop	1963	4001585
dummy_pop	1964	4064868
dummy_pop	1965	4107196
dummy_pop	1966	4162451
dummy_pop	1969	4358888
dummy_pop	1972	4453449
dummy_pop	1975	4523253
dummy_pop	1979	4602626
dummy_pop	1980	4655635
dummy_pop	1981	4713568
dummy_pop	1982	4773655
dummy_pop	1983	4815380
dummy_pop	1984	4865319
dummy_pop	1985	4910168
dummy_pop	1986	4955648
dummy_pop	1987	5010450
dummy_pop	1988	5057829
dummy_pop	1989	5102752
dummy_pop	1990	5156697
dummy_pop	1991	5224260
dummy_pop	1992	5261459
dummy_pop	1993	5301173
dummy_pop	1994	5334975
dummy_pop	1995	5362253
dummy_pop	1996	5376070
dummy_pop	1997	5387290
dummy_pop	1998	5417070
dummy_pop	1999	5454490
dummy_pop	2000	5492367
dummy_pop	2001	5561346
dummy_pop	2002	5618156
dummy_pop	2003	5670923
dummy_pop	2004	5654343
dummy_pop	2005	5700483
dummy_pop	2006	5753311
dummy_pop	2007	5833888
dummy_pop	2008	5931530
dummy_pop	2009	6009032
dummy_pop	2010	6081182
dummy_pop	2011	6151881
dummy_pop	2012	6221203
dummy_pop	2013	6303138
dummy_pop	2014	6383544
dummy_pop	2015	6452139
dummy_pop	2016	6518435
dummy_pop	2017	6560173

foreach var in n_all_sum_c1 n_lrg_sum_c1 n_sml_sum_c1 i_all_c1 i_lrg_c1 i_sml_c1 {
	gen `var' = 0
}


** Create population datases: add NR & director data
append using "$path\02_Processed_data\21_Population\NR&directors.dta"
compress

// year within legislative period
gen yrlegperiod = .
replace yrlegperiod = 4 if year == 1931
replace yrlegperiod = 4 if year == 1935
replace yrlegperiod = 4 if year == 1939
replace yrlegperiod = 4 if year == 1943
replace yrlegperiod = 4 if year == 1947
replace yrlegperiod = 4 if year == 1951
replace yrlegperiod = 4 if year == 1955
replace yrlegperiod = 4 if year == 1959
replace yrlegperiod = 4 if year == 1963
replace yrlegperiod = 4 if year == 1967
replace yrlegperiod = 4 if year == 1971
replace yrlegperiod = 4 if year == 1975
replace yrlegperiod = 4 if year == 1979
replace yrlegperiod = 4 if year == 1983
replace yrlegperiod = 4 if year == 1987
replace yrlegperiod = 4 if year == 1991
replace yrlegperiod = 4 if year == 1995
replace yrlegperiod = 4 if year == 1999
replace yrlegperiod = 4 if year == 2003
replace yrlegperiod = 4 if year == 2007
replace yrlegperiod = 4 if year == 2011
replace yrlegperiod = 4 if year == 2015

replace yrlegperiod = 1 if year == 1932
replace yrlegperiod = 1 if year == 1936
replace yrlegperiod = 1 if year == 1940
replace yrlegperiod = 1 if year == 1944
replace yrlegperiod = 1 if year == 1948
replace yrlegperiod = 1 if year == 1952
replace yrlegperiod = 1 if year == 1956
replace yrlegperiod = 1 if year == 1960
replace yrlegperiod = 1 if year == 1964
replace yrlegperiod = 1 if year == 1968
replace yrlegperiod = 1 if year == 1972
replace yrlegperiod = 1 if year == 1976
replace yrlegperiod = 1 if year == 1980
replace yrlegperiod = 1 if year == 1984
replace yrlegperiod = 1 if year == 1988
replace yrlegperiod = 1 if year == 1992
replace yrlegperiod = 1 if year == 1996
replace yrlegperiod = 1 if year == 2000
replace yrlegperiod = 1 if year == 2004
replace yrlegperiod = 1 if year == 2008
replace yrlegperiod = 1 if year == 2012
replace yrlegperiod = 1 if year == 2016

replace yrlegperiod = 2 if year == 1933
replace yrlegperiod = 2 if year == 1937
replace yrlegperiod = 2 if year == 1941
replace yrlegperiod = 2 if year == 1945
replace yrlegperiod = 2 if year == 1949
replace yrlegperiod = 2 if year == 1953
replace yrlegperiod = 2 if year == 1957
replace yrlegperiod = 2 if year == 1961
replace yrlegperiod = 2 if year == 1965
replace yrlegperiod = 2 if year == 1969
replace yrlegperiod = 2 if year == 1973
replace yrlegperiod = 2 if year == 1977
replace yrlegperiod = 2 if year == 1981
replace yrlegperiod = 2 if year == 1985
replace yrlegperiod = 2 if year == 1989
replace yrlegperiod = 2 if year == 1993
replace yrlegperiod = 2 if year == 1997
replace yrlegperiod = 2 if year == 2001
replace yrlegperiod = 2 if year == 2005
replace yrlegperiod = 2 if year == 2009
replace yrlegperiod = 2 if year == 2013
replace yrlegperiod = 2 if year == 2017

replace yrlegperiod = 3 if year == 1934
replace yrlegperiod = 3 if year == 1938
replace yrlegperiod = 3 if year == 1942
replace yrlegperiod = 3 if year == 1946
replace yrlegperiod = 3 if year == 1950
replace yrlegperiod = 3 if year == 1954
replace yrlegperiod = 3 if year == 1958
replace yrlegperiod = 3 if year == 1962
replace yrlegperiod = 3 if year == 1966
replace yrlegperiod = 3 if year == 1970
replace yrlegperiod = 3 if year == 1974
replace yrlegperiod = 3 if year == 1978
replace yrlegperiod = 3 if year == 1982
replace yrlegperiod = 3 if year == 1986
replace yrlegperiod = 3 if year == 1990
replace yrlegperiod = 3 if year == 1994
replace yrlegperiod = 3 if year == 1998
replace yrlegperiod = 3 if year == 2002
replace yrlegperiod = 3 if year == 2006
replace yrlegperiod = 3 if year == 2010
replace yrlegperiod = 3 if year == 2014


sort year ID PID
order ID PID year yrlegperiod ellegperiod
save "$path\02_Processed_data\21_Population\population_data.dta", replace

erase "$path\02_Processed_data\21_Population\dummy_pop.dta"
erase "$path\02_Processed_data\21_Population\pop_res.dta"
erase "$path\02_Processed_data\21_Population\n_directors.dta"
erase "$path\02_Processed_data\21_Population\NR&directors.dta"
