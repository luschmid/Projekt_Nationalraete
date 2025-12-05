* This do-file sets up the origin data etc. from Florence; see email from
* 8.12.2017

use "$path\07_Archive\01_Elections_1931_1975\nationalraete_1931_2015_mit_Heimatort_und_Job.dta", ///
	clear
replace ID="OW-1995-9001" if name_BU=="Durrer" & firstname_BU=="Rudolf O." ///
	& canton=="OW" & year==1995
replace ID="UR-1995-9001" if name_BU=="Merminod" & firstname_BU=="Yves" ///
	& canton=="UR" & year==1995
* Note: According to email from Florence from 8.12.2017 there are two missing 
* candidates in the BfS data for which we add Bundesblatt information.

* Note: we drop the following variables (among others)
* Mehrere: indicates cases with amibuous origin name(s) for which we use all 
* possible origins (e.g., Rüti ZH and Rüti BE)
* FNB: indiciates cases for which Florence tried to clarify ambiguous origin 
* names with Familiennamenbuch
* HistNummer: historical municipality number

keep ID year canton Gemeindename* BFSGdenummer* job_BU votes_BU birthyear_BU
save "$path\01_Raw_data\01_Elections_1931_1975\origins_and_jobs_tomerge.dta", ///
	replace
