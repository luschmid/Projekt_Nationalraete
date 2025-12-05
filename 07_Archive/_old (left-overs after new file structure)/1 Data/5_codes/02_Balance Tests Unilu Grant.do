
******************************************************
* (0) Set directories
******************************************************

clear
clear matrix
set mem 500m
set more off


	/*** DEFINE MAIN PATH ***/
	capture cd "C:\Dropbox\Incumbency Advantage"
	if _rc==0{
		global hauptpfad "C:\Dropbox\Projekt Nationalräte\1 Data"
		global datapath "C:\Dropbox\Projekt Nationalräte\1 Data\1_data"
		global path_original_data "C:\Dropbox\Projekt Nationalräte\1 Data\0_original_data"
		global path_fuzzy_data "C:\Dropbox\Projekt Nationalräte\1 Data\1_data\fuzzy merge datasets"
		global path_fuzzy "C:\Dropbox\Projekt Nationalräte\1 Data\3_fuzzy merge do files"

	}
	
    capture cd "D:\SchmidLu\Dropbox\Projekt Nationalräte"

	if _rc==0{ 
	global hauptpfad "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data"
		global datapath "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data"
		global path_original_data "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\0_original_data"
		global path_fuzzy_data "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data\fuzzy merge datasets"
		global path_fuzzy "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\3_fuzzy merge do files"
	}
	
*************************************		
* (A) Prepare fuzzy merge data
*************************************	
clear
cd "$datapath\6 fuzzy merge final check"
import excel "Final_Check.xlsx",  first 
//saveold Final_Check.dta, version(13)

rename Vorname firstname
rename Nachname name
rename Kanton canton
rename Listenbezeichnung list
rename Gemeinde municipality
rename Geburtsjahr birthyear 
rename Wahljahr year

keep firstname name canton list birthyear municipality year id1

gen canton_merge=canton
replace canton_merge="BEJU" if canton=="BE" | canton=="JU"

replace name="Haffter" if  name=="Hafter" & firstname=="Arthur" &	birthyear==1927 & year==1971
replace firstname="Karl" if  name=="Läuchli" & firstname=="Carl" &	birthyear==1881 & year==1931
replace birthyear=1901 if  name=="Revaclier" & firstname=="François" &	birthyear==1903 & year==1963
replace birthyear=1874 if  name=="Varesi" & firstname=="Giovanni" &	birthyear==1864 & year==1943
replace birthyear=1890 if  name=="Holliger" & firstname=="Hans" &	birthyear==1892 & year==1939
replace birthyear=1908 if  name=="Arnold" & firstname=="Max" &	birthyear==1909 & year==1951	
replace birthyear=1873 if  name=="Ryffel-Schiess" & firstname=="Albert" &	birthyear==1878 & year==1935
replace name="Wäfler" if  name=="Wäffler" & firstname=="Markus" &	birthyear==1948 & year==1983



*************************************		
* (B) Prepare electoral results data
*************************************	



// reverse all name changes made by coders in Fribourg


preserve
cd "$path_fuzzy_data"
use nationalraete_1931_2015_fuzzy.dta, clear

replace birthyear=1942 if firstname=="Brigitte" & name=="Schmid" & year==1983 & canton=="ZH" 
replace birthyear=1942 if firstname=="Brigitte" & name=="Schmid" & year==1987 & canton=="ZH" 
replace birthyear=1951 if firstname=="Christoph" & name=="Keller" & year==1979 & canton=="ZH" 
replace birthyear=1909 if firstname=="Christoph" & name=="Wolfensberger" & year==1975 & canton=="ZH" 
replace birthyear=1960 if firstname=="Daniel" & name=="Holzreuter" & year==1995 & canton=="ZH" 
replace birthyear=1945 if firstname=="Françoise" & name=="Glatz-Studer" & year==2015 & canton=="VD" 
replace birthyear=1956 if firstname=="Guy" & name=="Morin" & year==1991 & canton=="BS" 
replace birthyear=1943 if firstname=="Kurt" & name=="Schreiber" & year==1991 & canton=="ZH" 
replace birthyear=1946 if firstname=="Marianne" &  name=="Jaccard" & year==1987 & canton=="VD" 
replace birthyear=1946 if firstname=="Marianne" &  name=="Jaccard" & year==1991 & canton=="VD" 
replace birthyear=1944 if firstname=="Niklaus" &  name=="Scherr" & year==1979 & canton=="ZH" 
replace birthyear=1941 if firstname=="Rudolf" &  name=="Bautz" & year==1983 & canton=="ZH" 


replace firstname="Margrit" if firstname=="Margrith" & name=="Kolp" & year==1999 & canton=="ZH" 
replace name="Wäfler"  if firstname=="Markus" & name=="Wäffler"  & year==1983 & canton=="ZH" 
replace name="Sturny"  if firstname=="Max" & name=="Sturni"  & year==1979 & canton=="ZH" 
replace firstname="Regine"  if firstname=="Regina" & name=="Aeppli"  & year==1995 & canton=="ZH" 
replace firstname="Stefan"  if firstname=="Stephan" & name=="Dollenmeier"  & year==2011 & canton=="ZH" 

replace name="Stebler" if name=="Stehler" & year==1967 & birthyear==1924
replace name="Degen" if  name=="Degea" & year==1967 & birthyear==1933
replace name="Corbat" if name=="Corbaz"& year==1967 & birthyear==1925
replace name="Hafter" if name=="Haffter" & year==1967 & birthyear==1927
replace name="Hugli" if  name=="Huggli" & year==1971 & birthyear==1915
replace name="Günthard" if  name=="Günthart" & year==1939 & birthyear==1880
replace name="Bosshart" if  name=="Bosshardt" & year==1939 & birthyear==1901
replace name="Bänninger" if name=="Bänniger"  & year==1935 & birthyear==1898
replace name="König" if name=="Kong" & year==1963 & birthyear==1916
replace name="Sturny" if name=="Sturni" 


replace firstname="Patrick" if   name=="Mächler"   &	birthyear==1983 & year==2011
replace firstname="Massimiliano" if name=="Ay" 	 &	birthyear==1982 & year==2011
replace firstname="Sergio" if  name=="Arigoni" &		birthyear==1965 & year==2011
replace firstname="Werner 2" if name=="Müller" &	 	birthyear==1932 & year==1975 & firstname=="Werner"
replace firstname="Willi" if  name=="Walker" & 	birthyear==1922 & year==1975


replace birthyear=1893 if  name=="Hofmann" & firstname=="Paul" &	birthyear==1898 & year==1939


gen canton_merge=canton
replace canton_merge="BEJU" if canton=="BE" | canton=="JU"

cd "$path_fuzzy_data"
saveold electionresults_1931_2015.dta, version(13) replace
restore


*********************************************************
* (C) Merge fuzzy merge data to electoral results data
*********************************************************

cd "$path_fuzzy_data"
merge 1:1 firstname name canton_merge  birthyear municipality year using electionresults_1931_2015.dta, gen(merge_fuzzy_election_results)


gen id=id1
gen votes_total=votes

bysort canton year list: egen list_votes_total=sum(votes_total)
bysort canton year : egen votes_total_ct=sum(votes_total)

gen candidate_votes_perc= votes_total/list_votes_total*100 // generate individual vote share in terms of total list votes
gen candidate_votes_perc_ct= votes_total/votes_total_ct*100 // generate individual vote share in terms of total ctl votes

*(ii) ranking candidates on a list


bysort id: gen nr_participation=_n

gsort year canton  list -votes_total
bysort year canton  list: gen candidate_rank=_n // rank of candidate on list


bysort year canton  list: egen elected_perlist= sum(elected) // indicator if anybody is elected on list
bysort year canton : egen elected_ct= sum(elected) 
		
*(iii) identifying marginal guy

sort  year canton  list candidate_rank
gen elect_rank=elected*candidate_rank		// interaction rank and elected
bysort year canton  list: egen elect_rank_max=max(elect_rank)
gen candidate_dist=candidate_rank-elect_rank_max // distance to marginal guy in positions

gen marginal=0 // marginal guy from below (for all non-elected)
replace marginal=1 if candidate_dist==0

gen marginal_above=0 // marginal guy from above (for all elected)
replace marginal_above=1 if candidate_dist==1

*(iv) generating running variable

gen candidate_votes_perc_marg=. 
replace candidate_votes_perc_marg=candidate_votes_perc if marginal==1

gen candidate_votes_perc_marg_above=. 
replace candidate_votes_perc_marg_above=candidate_votes_perc if marginal_above==1

bysort year canton  list:egen  candidate_votes_perc_marg_all=max(candidate_votes_perc_marg) 
bysort year canton list:egen  candidate_votes_perc_mar_ab_all=min(candidate_votes_perc_marg_above) 

gen candidate_diff_marginal=candidate_votes_perc-candidate_votes_perc_marg_all if elected==0
replace candidate_diff_marginal=candidate_votes_perc-candidate_votes_perc_mar_ab_all if elected==1


drop job // not able to save in stata 12 format

// temporary: drop duplicates

drop if id1==""
bysort id year: gen indimax=_N
drop if indimax>1

// end of temporary


egen id_num=group(id)
egen time_id=group(year)


sort id_num year
xtset id_num time_id

gen elected_lag1=L.elected
gen elected_lag2=L2.elected


cd "$datapath"
saveold nationarat_all.dta, version(12) replace
