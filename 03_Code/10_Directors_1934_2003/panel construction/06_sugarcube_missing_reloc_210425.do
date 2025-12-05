version 17
clear all
set more off
cap log close

global path "C:\Users\schmutzy\Dropbox\Projekt Nationalräte\" /* Office */
global dump "C:\Users\schmutzy\Desktop\Projekt Nationalräte\" /* Office Dumpster */
global ED "C:\Users\schmutzy\Dropbox\Travaux_EmilieDousse\" /* Office ED */


* Use latest version in previous matching round (same mun)
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_38.dta", clear

**************************
* 1934 & 1943

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total dup_panel maxpanel minpanel group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year
local agrp 1934 1943 /*1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003*/
foreach a of local agrp{
	preserve
		keep if year == `a'
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_`a'_remaining_2.dta", replace
	restore
	drop if year == `a'
}


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1934 & 1943
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1943_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.85
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_1943_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_1943_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1943_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_1943_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_1943_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1934_1943_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1934 to 1943
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_38.dta", clear

 replace UCID_1934_2003 = 	198331	 if UCID_1934_2003 == 	43964
 replace UCID_1934_2003 = 	10937	 if UCID_1934_2003 == 	47681
 replace UCID_1934_2003 = 	56125	 if UCID_1934_2003 == 	55970
 replace UCID_1934_2003 = 	430702	 if UCID_1934_2003 == 	56309
 replace UCID_1934_2003 = 	57620	 if UCID_1934_2003 == 	57693
 replace UCID_1934_2003 = 	411851	 if UCID_1934_2003 == 	89474
 replace UCID_1934_2003 = 	213953	 if UCID_1934_2003 == 	95049
 replace UCID_1934_2003 = 	215701	 if UCID_1934_2003 == 	95581
 replace UCID_1934_2003 = 	182050	 if UCID_1934_2003 == 	104190
 replace UCID_1934_2003 = 	60614	 if UCID_1934_2003 == 	104869
 replace UCID_1934_2003 = 	274450	 if UCID_1934_2003 == 	104982
 replace UCID_1934_2003 = 	410245	 if UCID_1934_2003 == 	106177
 replace UCID_1934_2003 = 	2624	 if UCID_1934_2003 == 	107672
 replace UCID_1934_2003 = 	219760	 if UCID_1934_2003 == 	114249
 replace UCID_1934_2003 = 	214320	 if UCID_1934_2003 == 	117127
 replace UCID_1934_2003 = 	130379	 if UCID_1934_2003 == 	125789
 replace UCID_1934_2003 = 	393085	 if UCID_1934_2003 == 	126648
 replace UCID_1934_2003 = 	413426	 if UCID_1934_2003 == 	130790
 replace UCID_1934_2003 = 	135508	 if UCID_1934_2003 == 	135515
 replace UCID_1934_2003 = 	235628	 if UCID_1934_2003 == 	146578
 replace UCID_1934_2003 = 	151942	 if UCID_1934_2003 == 	153539
 replace UCID_1934_2003 = 	263601	 if UCID_1934_2003 == 	160448
 replace UCID_1934_2003 = 	378577	 if UCID_1934_2003 == 	162260
 replace UCID_1934_2003 = 	209698	 if UCID_1934_2003 == 	209489
 replace UCID_1934_2003 = 	345198	 if UCID_1934_2003 == 	216385
 replace UCID_1934_2003 = 	226259	 if UCID_1934_2003 == 	216534
 replace UCID_1934_2003 = 	228624	 if UCID_1934_2003 == 	217271
 replace UCID_1934_2003 = 	159706	 if UCID_1934_2003 == 	217595
 replace UCID_1934_2003 = 	268856	 if UCID_1934_2003 == 	218728
 replace UCID_1934_2003 = 	48082	 if UCID_1934_2003 == 	231109
 replace UCID_1934_2003 = 	233672	 if UCID_1934_2003 == 	233550
 replace UCID_1934_2003 = 	188674	 if UCID_1934_2003 == 	241274
 replace UCID_1934_2003 = 	416369	 if UCID_1934_2003 == 	241275
 replace UCID_1934_2003 = 	132774	 if UCID_1934_2003 == 	255299
 replace UCID_1934_2003 = 	395762	 if UCID_1934_2003 == 	263842
 replace UCID_1934_2003 = 	313725	 if UCID_1934_2003 == 	307935
 replace UCID_1934_2003 = 	320583	 if UCID_1934_2003 == 	311880
 replace UCID_1934_2003 = 	320577	 if UCID_1934_2003 == 	311880
 replace UCID_1934_2003 = 	322258	 if UCID_1934_2003 == 	321842
 replace UCID_1934_2003 = 	322343	 if UCID_1934_2003 == 	321842
 replace UCID_1934_2003 = 	324327	 if UCID_1934_2003 == 	324468
 replace UCID_1934_2003 = 	324323	 if UCID_1934_2003 == 	324468
 replace UCID_1934_2003 = 	411131	 if UCID_1934_2003 == 	331335
 replace UCID_1934_2003 = 	328501	 if UCID_1934_2003 == 	332449
 replace UCID_1934_2003 = 	343357	 if UCID_1934_2003 == 	333920
 replace UCID_1934_2003 = 	394919	 if UCID_1934_2003 == 	335149
 replace UCID_1934_2003 = 	349441	 if UCID_1934_2003 == 	349668
 replace UCID_1934_2003 = 	41621	 if UCID_1934_2003 == 	353001
 replace UCID_1934_2003 = 	366565	 if UCID_1934_2003 == 	361288
 replace UCID_1934_2003 = 	408233	 if UCID_1934_2003 == 	362819
 replace UCID_1934_2003 = 	389860	 if UCID_1934_2003 == 	363320
 replace UCID_1934_2003 = 	367825	 if UCID_1934_2003 == 	366274
 replace UCID_1934_2003 = 	333630	 if UCID_1934_2003 == 	366799
 replace UCID_1934_2003 = 	367825	 if UCID_1934_2003 == 	367036
 replace UCID_1934_2003 = 	369459	 if UCID_1934_2003 == 	368158
 replace UCID_1934_2003 = 	397206	 if UCID_1934_2003 == 	369148
 replace UCID_1934_2003 = 	373319	 if UCID_1934_2003 == 	369710
 replace UCID_1934_2003 = 	386241	 if UCID_1934_2003 == 	369883
 replace UCID_1934_2003 = 	420600	 if UCID_1934_2003 == 	369885
 replace UCID_1934_2003 = 	397232	 if UCID_1934_2003 == 	370650
 replace UCID_1934_2003 = 	23261	 if UCID_1934_2003 == 	374743
 replace UCID_1934_2003 = 	387143	 if UCID_1934_2003 == 	390258
 replace UCID_1934_2003 = 	362319	 if UCID_1934_2003 == 	394257
 replace UCID_1934_2003 = 	376996	 if UCID_1934_2003 == 	402973
 replace UCID_1934_2003 = 	389644	 if UCID_1934_2003 == 	405708
 replace UCID_1934_2003 = 	409300	 if UCID_1934_2003 == 	406088
 replace UCID_1934_2003 = 	420700	 if UCID_1934_2003 == 	406417
 replace UCID_1934_2003 = 	389686	 if UCID_1934_2003 == 	409100
 replace UCID_1934_2003 = 	215600	 if UCID_1934_2003 == 	409733
 replace UCID_1934_2003 = 	430100	 if UCID_1934_2003 == 	409798
 replace UCID_1934_2003 = 	380388	 if UCID_1934_2003 == 	409812
 replace UCID_1934_2003 = 	275559	 if UCID_1934_2003 == 	410084
 replace UCID_1934_2003 = 	390848	 if UCID_1934_2003 == 	410115
 replace UCID_1934_2003 = 	105791	 if UCID_1934_2003 == 	410117
 replace UCID_1934_2003 = 	23604	 if UCID_1934_2003 == 	410252
 replace UCID_1934_2003 = 	242933	 if UCID_1934_2003 == 	410274
 replace UCID_1934_2003 = 	322258	 if UCID_1934_2003 == 	410381
 replace UCID_1934_2003 = 	322343	 if UCID_1934_2003 == 	410381
 replace UCID_1934_2003 = 	130474	 if UCID_1934_2003 == 	410485
 replace UCID_1934_2003 = 	82962	 if UCID_1934_2003 == 	410728
 replace UCID_1934_2003 = 	209663	 if UCID_1934_2003 == 	410773
 replace UCID_1934_2003 = 	352941	 if UCID_1934_2003 == 	410774
 replace UCID_1934_2003 = 	41666	 if UCID_1934_2003 == 	410908
 replace UCID_1934_2003 = 	328439	 if UCID_1934_2003 == 	411114
 replace UCID_1934_2003 = 	89744	 if UCID_1934_2003 == 	411252
 replace UCID_1934_2003 = 	119949	 if UCID_1934_2003 == 	411351
 replace UCID_1934_2003 = 	390827	 if UCID_1934_2003 == 	411544
 replace UCID_1934_2003 = 	130934	 if UCID_1934_2003 == 	411794
 replace UCID_1934_2003 = 	364821	 if UCID_1934_2003 == 	412897
 replace UCID_1934_2003 = 	323798	 if UCID_1934_2003 == 	412985
 replace UCID_1934_2003 = 	335883	 if UCID_1934_2003 == 	412995
 replace UCID_1934_2003 = 	106817	 if UCID_1934_2003 == 	413074
 replace UCID_1934_2003 = 	8541	 if UCID_1934_2003 == 	413113
 replace UCID_1934_2003 = 	346095	 if UCID_1934_2003 == 	413133
 replace UCID_1934_2003 = 	97668	 if UCID_1934_2003 == 	413256
 replace UCID_1934_2003 = 	361213	 if UCID_1934_2003 == 	413282
 replace UCID_1934_2003 = 	373132	 if UCID_1934_2003 == 	413306
 replace UCID_1934_2003 = 	195361	 if UCID_1934_2003 == 	413385
 replace UCID_1934_2003 = 	216113	 if UCID_1934_2003 == 	413427
 replace UCID_1934_2003 = 	396475	 if UCID_1934_2003 == 	413471
 replace UCID_1934_2003 = 	402271	 if UCID_1934_2003 == 	413471
 replace UCID_1934_2003 = 	388105	 if UCID_1934_2003 == 	413492
 replace UCID_1934_2003 = 	274971	 if UCID_1934_2003 == 	413513
 replace UCID_1934_2003 = 	401960	 if UCID_1934_2003 == 	413564
 replace UCID_1934_2003 = 	75666	 if UCID_1934_2003 == 	413673
 replace UCID_1934_2003 = 	397232	 if UCID_1934_2003 == 	413678
 replace UCID_1934_2003 = 	371230	 if UCID_1934_2003 == 	413715
 replace UCID_1934_2003 = 	369516	 if UCID_1934_2003 == 	413758
 replace UCID_1934_2003 = 	330303	 if UCID_1934_2003 == 	413763
 replace UCID_1934_2003 = 	275250	 if UCID_1934_2003 == 	413769
 replace UCID_1934_2003 = 	395649	 if UCID_1934_2003 == 	413788
 replace UCID_1934_2003 = 	264158	 if UCID_1934_2003 == 	413800
 replace UCID_1934_2003 = 	27321	 if UCID_1934_2003 == 	413813
 replace UCID_1934_2003 = 	349145	 if UCID_1934_2003 == 	413845
 replace UCID_1934_2003 = 	165039	 if UCID_1934_2003 == 	413873
 replace UCID_1934_2003 = 	329098	 if UCID_1934_2003 == 	413875
 replace UCID_1934_2003 = 	345198	 if UCID_1934_2003 == 	413940
 replace UCID_1934_2003 = 	288623	 if UCID_1934_2003 == 	413961
 replace UCID_1934_2003 = 	84475	 if UCID_1934_2003 == 	414035
 replace UCID_1934_2003 = 	93738	 if UCID_1934_2003 == 	414089
 replace UCID_1934_2003 = 	334166	 if UCID_1934_2003 == 	414392
 replace UCID_1934_2003 = 	345071	 if UCID_1934_2003 == 	414392
 replace UCID_1934_2003 = 	264356	 if UCID_1934_2003 == 	414628
 replace UCID_1934_2003 = 	104889	 if UCID_1934_2003 == 	414793
 replace UCID_1934_2003 = 	276919	 if UCID_1934_2003 == 	414796
 replace UCID_1934_2003 = 	88283	 if UCID_1934_2003 == 	414854
 replace UCID_1934_2003 = 	427760	 if UCID_1934_2003 == 	414988
 replace UCID_1934_2003 = 	106817	 if UCID_1934_2003 == 	415012
 replace UCID_1934_2003 = 	423851	 if UCID_1934_2003 == 	415017
 replace UCID_1934_2003 = 	343991	 if UCID_1934_2003 == 	415081
 replace UCID_1934_2003 = 	119890	 if UCID_1934_2003 == 	415104
 replace UCID_1934_2003 = 	121326	 if UCID_1934_2003 == 	415104
 replace UCID_1934_2003 = 	408158	 if UCID_1934_2003 == 	415280
 replace UCID_1934_2003 = 	367010	 if UCID_1934_2003 == 	415381
 replace UCID_1934_2003 = 	255075	 if UCID_1934_2003 == 	415388
 replace UCID_1934_2003 = 	277854	 if UCID_1934_2003 == 	415622
 replace UCID_1934_2003 = 	109742	 if UCID_1934_2003 == 	415709
 replace UCID_1934_2003 = 	223415	 if UCID_1934_2003 == 	415764
 replace UCID_1934_2003 = 	367519	 if UCID_1934_2003 == 	415828
 replace UCID_1934_2003 = 	367615	 if UCID_1934_2003 == 	415828
 replace UCID_1934_2003 = 	903	 if UCID_1934_2003 == 	415995
 replace UCID_1934_2003 = 	389941	 if UCID_1934_2003 == 	416047
 replace UCID_1934_2003 = 	390541	 if UCID_1934_2003 == 	416167
 replace UCID_1934_2003 = 	388066	 if UCID_1934_2003 == 	416167
 replace UCID_1934_2003 = 	228128	 if UCID_1934_2003 == 	416172
 replace UCID_1934_2003 = 	27383	 if UCID_1934_2003 == 	416172
 replace UCID_1934_2003 = 	870	 if UCID_1934_2003 == 	416203
 replace UCID_1934_2003 = 	388836	 if UCID_1934_2003 == 	416249
 replace UCID_1934_2003 = 	115499	 if UCID_1934_2003 == 	416340
 replace UCID_1934_2003 = 	222570	 if UCID_1934_2003 == 	416342
 replace UCID_1934_2003 = 	123678	 if UCID_1934_2003 == 	416429
 replace UCID_1934_2003 = 	109762	 if UCID_1934_2003 == 	416441
 replace UCID_1934_2003 = 	389897	 if UCID_1934_2003 == 	429383
 replace UCID_1934_2003 = 	313725	 if UCID_1934_2003 == 	429923
 replace UCID_1934_2003 = 	313725	 if UCID_1934_2003 == 	429969
 replace UCID_1934_2003 = 	56125	 if UCID_1934_2003 == 	429974
 replace UCID_1934_2003 = 	409645	 if UCID_1934_2003 == 	430024 
 
* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}

	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_39.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_38.dta"

********************************************************************************
* 1943 & 1960

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_39.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1943 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1943
		drop if dup_panel == minpanel
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1943_remaining_2.dta", replace
	restore

*Second year 1960 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1960
		drop if dup_panel == maxpanel
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1960_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1943 & 1960
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1943_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1960_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1943_1960_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1943_1960_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1943_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1960_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1943_1960_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1943_1960_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1943_1960_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1943 to 1960
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_39.dta", clear

 replace UCID_1934_2003 = 	377973	 if UCID_1934_2003 == 	22633
 replace UCID_1934_2003 = 	23559	 if UCID_1934_2003 == 	23563
 replace UCID_1934_2003 = 	365896	 if UCID_1934_2003 == 	42501
 replace UCID_1934_2003 = 	58516	 if UCID_1934_2003 == 	84475
 replace UCID_1934_2003 = 	11408	 if UCID_1934_2003 == 	121135
 replace UCID_1934_2003 = 	2614	 if UCID_1934_2003 == 	217316
 replace UCID_1934_2003 = 	402339	 if UCID_1934_2003 == 	362641
 replace UCID_1934_2003 = 	43078	 if UCID_1934_2003 == 	397206
 replace UCID_1934_2003 = 	370910	 if UCID_1934_2003 == 	397232
 replace UCID_1934_2003 = 	409824	 if UCID_1934_2003 == 	410245
 replace UCID_1934_2003 = 	405327	 if UCID_1934_2003 == 	415879
 replace UCID_1934_2003 = 	129827	 if UCID_1934_2003 == 	426130
 replace UCID_1934_2003 = 	130156	 if UCID_1934_2003 == 	426249
 replace UCID_1934_2003 = 	131018	 if UCID_1934_2003 == 	428792

 * Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_40.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_39.dta"

********************************************************************************
* 1960 & 1962

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_40.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1960 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1960
		drop if dup_panel == minpanel
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1960_remaining_2.dta", replace
	restore

*Second year 1962 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1962
		drop if dup_panel == maxpanel
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1962_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1960 & 1962
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1960_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1962_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1960_1962_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1960_1962_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1960_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1962_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1960_1962_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1960_1962_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1960_1962_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1960 to 1962
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_40.dta", clear

replace UCID_1934_2003 = 9555 if UCID_1934_2003 == 146611

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_41.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_40.dta"

********************************************************************************
* 1962 & 1963

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_41.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1962 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1962
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1962
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1962_remaining_2.dta", replace
	restore

*Second year 1963 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1963
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1963
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1963_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1962 & 1963
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1962_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1963_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1962_1963_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1962_1963_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1962_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1963_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1962_1963_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1962_1963_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1962_1963_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1962 to 1963
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_41.dta", clear

 replace UCID_1934_2003 = 	427912	 if UCID_1934_2003 == 	32456
 replace UCID_1934_2003 = 	429219	 if UCID_1934_2003 == 	47786
 replace UCID_1934_2003 = 	70593	 if UCID_1934_2003 == 	62894
 replace UCID_1934_2003 = 	266148	 if UCID_1934_2003 == 	83751
 replace UCID_1934_2003 = 	72838	 if UCID_1934_2003 == 	92300
 replace UCID_1934_2003 = 	413844	 if UCID_1934_2003 == 	94176
 replace UCID_1934_2003 = 	59790	 if UCID_1934_2003 == 	152759
 replace UCID_1934_2003 = 	410264	 if UCID_1934_2003 == 	177794
 replace UCID_1934_2003 = 	345837	 if UCID_1934_2003 == 	322277
 replace UCID_1934_2003 = 	257528	 if UCID_1934_2003 == 	385137
 
* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_42.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_41.dta"

*******************************************************************************
* 1963 & 1964

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_42.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1963 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1963
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1963
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1963_remaining_2.dta", replace
	restore

*Second year 1964 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1964
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1964
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1964_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1963 & 1964
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1963_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1964_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1963_1964_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1963_1964_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1963_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1964_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1963_1964_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1963_1964_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1963_1964_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1963 to 1964
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_42.dta", clear

 replace UCID_1934_2003 = 	27159	 if UCID_1934_2003 == 	27126
 replace UCID_1934_2003 = 	196531	 if UCID_1934_2003 == 	111877
 replace UCID_1934_2003 = 	409415	 if UCID_1934_2003 == 	119104
 replace UCID_1934_2003 = 	70103	 if UCID_1934_2003 == 	160765
 replace UCID_1934_2003 = 	230315	 if UCID_1934_2003 == 	209712
 replace UCID_1934_2003 = 	47538	 if UCID_1934_2003 == 	284498
 replace UCID_1934_2003 = 	32350	 if UCID_1934_2003 == 	296484
 replace UCID_1934_2003 = 	348029	 if UCID_1934_2003 == 	333151
 replace UCID_1934_2003 = 	8354	 if UCID_1934_2003 == 	366916
 replace UCID_1934_2003 = 	103845	 if UCID_1934_2003 == 	392660

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}

	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0	
	
drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_43.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_42.dta"

*******************************************************************************
* 1964 & 1965

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_43.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1964 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1964
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1964
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1964_remaining_2.dta", replace
	restore

*Second year 1965 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1965
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1965
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1965_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1964 & 1965
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1964_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1965_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1964_1965_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1964_1965_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1964_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1965_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1964_1965_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1964_1965_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1964_1965_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1964 to 1965
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_43.dta", clear

 replace UCID_1934_2003 = 	145933	 if UCID_1934_2003 == 	24153
 replace UCID_1934_2003 = 	357361	 if UCID_1934_2003 == 	41917
 replace UCID_1934_2003 = 	24426	 if UCID_1934_2003 == 	90018
 replace UCID_1934_2003 = 	226186	 if UCID_1934_2003 == 	125217
 replace UCID_1934_2003 = 	47329	 if UCID_1934_2003 == 	321915
 replace UCID_1934_2003 = 	346531	 if UCID_1934_2003 == 	335113
 replace UCID_1934_2003 = 	46295	 if UCID_1934_2003 == 	369332
 replace UCID_1934_2003 = 	370354	 if UCID_1934_2003 == 	381578
 replace UCID_1934_2003 = 	412420	 if UCID_1934_2003 == 	416166
 replace UCID_1934_2003 = 	45231	 if UCID_1934_2003 == 	425766

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_44.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_43.dta"

*******************************************************************************
* 1965 & 1966

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_44.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1965 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1965
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1965
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1965_remaining_2.dta", replace
	restore

*Second year 1966 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1966
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1966
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1966_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1965 & 1966
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1965_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1966_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1965_1966_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1965_1966_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1965_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1966_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1965_1966_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1965_1966_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1965_1966_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1965 to 1966
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_44.dta", clear

replace gdenr = 3851 if CID == 56678 & year == 1966
replace gdenr_2012 = 3851 if CID == 56678 & year == 1966
replace gdenr_2018 = 3851 if CID == 56678 & year == 1966

 replace UCID_1934_2003 = 	13543	 if UCID_1934_2003 == 	19793
 replace UCID_1934_2003 = 	42953	 if UCID_1934_2003 == 	44146
 replace UCID_1934_2003 = 	191895	 if UCID_1934_2003 == 	87615
 replace UCID_1934_2003 = 	199261	 if UCID_1934_2003 == 	185098
 replace UCID_1934_2003 = 	238619	 if UCID_1934_2003 == 	233685
 replace UCID_1934_2003 = 	258514	 if UCID_1934_2003 == 	259102
 replace UCID_1934_2003 = 	94671	 if UCID_1934_2003 == 	268642
 replace UCID_1934_2003 = 	32923	 if UCID_1934_2003 == 	277921
 replace UCID_1934_2003 = 	319917	 if UCID_1934_2003 == 	314513
 replace UCID_1934_2003 = 	181191	 if UCID_1934_2003 == 	325552
 replace UCID_1934_2003 = 	346654	 if UCID_1934_2003 == 	346877
 replace UCID_1934_2003 = 	171419	 if UCID_1934_2003 == 	383918
 replace UCID_1934_2003 = 	376797	 if UCID_1934_2003 == 	385049
 replace UCID_1934_2003 = 	373923	 if UCID_1934_2003 == 	400359
 replace UCID_1934_2003 = 	322176	 if UCID_1934_2003 == 	414400
 replace UCID_1934_2003 = 	150371	 if UCID_1934_2003 == 	414653
 replace UCID_1934_2003 = 	189257	 if UCID_1934_2003 == 	415705

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_45.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_44.dta"

*******************************************************************************
* 1966 & 1969

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_45.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1966 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1966
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1966
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1966_remaining_2.dta", replace
	restore

*Second year 1969 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1969
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1969
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1969_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1966 & 1969
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1966_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1969_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1966_1969_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1966_1969_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1966_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1969_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1966_1969_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1966_1969_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1966_1969_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1966 to 1969
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_45.dta", clear

 replace UCID_1934_2003 = 	16805	 if UCID_1934_2003 == 	3412
 replace UCID_1934_2003 = 	1388	 if UCID_1934_2003 == 	5970
 replace UCID_1934_2003 = 	11765	 if UCID_1934_2003 == 	25341
 replace UCID_1934_2003 = 	14094	 if UCID_1934_2003 == 	25922
 replace UCID_1934_2003 = 	27364	 if UCID_1934_2003 == 	26008
 replace UCID_1934_2003 = 	432970	 if UCID_1934_2003 == 	27967
 replace UCID_1934_2003 = 	289913	 if UCID_1934_2003 == 	34178
 replace UCID_1934_2003 = 	42596	 if UCID_1934_2003 == 	42280
 replace UCID_1934_2003 = 	337665	 if UCID_1934_2003 == 	46580
 replace UCID_1934_2003 = 	272066	 if UCID_1934_2003 == 	74697
 replace UCID_1934_2003 = 	54207	 if UCID_1934_2003 == 	81085
 replace UCID_1934_2003 = 	72281	 if UCID_1934_2003 == 	84858
 replace UCID_1934_2003 = 	179730	 if UCID_1934_2003 == 	87145
 replace UCID_1934_2003 = 	62976	 if UCID_1934_2003 == 	88640
 replace UCID_1934_2003 = 	66892	 if UCID_1934_2003 == 	96266
 replace UCID_1934_2003 = 	64481	 if UCID_1934_2003 == 	103582
 replace UCID_1934_2003 = 	125270	 if UCID_1934_2003 == 	114233
 replace UCID_1934_2003 = 	108429	 if UCID_1934_2003 == 	116612
 replace UCID_1934_2003 = 	12269	 if UCID_1934_2003 == 	121147
 replace UCID_1934_2003 = 	196404	 if UCID_1934_2003 == 	131221
 replace UCID_1934_2003 = 	137135	 if UCID_1934_2003 == 	140388
 replace UCID_1934_2003 = 	411820	 if UCID_1934_2003 == 	145283
 replace UCID_1934_2003 = 	318847	 if UCID_1934_2003 == 	175956
 replace UCID_1934_2003 = 	392619	 if UCID_1934_2003 == 	195742
 replace UCID_1934_2003 = 	3026	 if UCID_1934_2003 == 	213085
 replace UCID_1934_2003 = 	225623	 if UCID_1934_2003 == 	213249
 replace UCID_1934_2003 = 	183512	 if UCID_1934_2003 == 	216790
 replace UCID_1934_2003 = 	31084	 if UCID_1934_2003 == 	243266
 replace UCID_1934_2003 = 	246235	 if UCID_1934_2003 == 	246703
 replace UCID_1934_2003 = 	267078	 if UCID_1934_2003 == 	266059
 replace UCID_1934_2003 = 	411903	 if UCID_1934_2003 == 	266302
 replace UCID_1934_2003 = 	82610	 if UCID_1934_2003 == 	269113
 replace UCID_1934_2003 = 	263452	 if UCID_1934_2003 == 	307646
 replace UCID_1934_2003 = 	352544	 if UCID_1934_2003 == 	355719
 replace UCID_1934_2003 = 	405740	 if UCID_1934_2003 == 	374224
 replace UCID_1934_2003 = 	172729	 if UCID_1934_2003 == 	375034
 replace UCID_1934_2003 = 	42931	 if UCID_1934_2003 == 	382644
 replace UCID_1934_2003 = 	82444	 if UCID_1934_2003 == 	396367
 replace UCID_1934_2003 = 	386551	 if UCID_1934_2003 == 	405696
 replace UCID_1934_2003 = 	410939	 if UCID_1934_2003 == 	409119
 replace UCID_1934_2003 = 	419708	 if UCID_1934_2003 == 	409428
 replace UCID_1934_2003 = 	358081	 if UCID_1934_2003 == 	411493
 replace UCID_1934_2003 = 	386122	 if UCID_1934_2003 == 	420614
 replace UCID_1934_2003 = 	416865	 if UCID_1934_2003 == 	424729

 * Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_46.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_45.dta"

*******************************************************************************
* 1969 & 1972

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_46.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1969 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1969
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1969
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1969_remaining_2.dta", replace
	restore

*Second year 1972 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1972
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1972
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1972_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1969 & 1972
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1969_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1972_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1969_1972_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1969_1972_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1969_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1972_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1969_1972_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1969_1972_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1969_1972_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1969 to 1972
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_46.dta", clear

 replace UCID_1934_2003 = 	120638	 if UCID_1934_2003 == 	12254
 replace UCID_1934_2003 = 	27086	 if UCID_1934_2003 == 	25987
 replace UCID_1934_2003 = 	55227	 if UCID_1934_2003 == 	53560
 replace UCID_1934_2003 = 	55992	 if UCID_1934_2003 == 	57037
 replace UCID_1934_2003 = 	9248	 if UCID_1934_2003 == 	58094
 replace UCID_1934_2003 = 	295006	 if UCID_1934_2003 == 	62577
 replace UCID_1934_2003 = 	55857	 if UCID_1934_2003 == 	64793
 replace UCID_1934_2003 = 	87528	 if UCID_1934_2003 == 	67061
 replace UCID_1934_2003 = 	405470	 if UCID_1934_2003 == 	82665
 replace UCID_1934_2003 = 	25511	 if UCID_1934_2003 == 	97938
 replace UCID_1934_2003 = 	265652	 if UCID_1934_2003 == 	105722
 replace UCID_1934_2003 = 	131783	 if UCID_1934_2003 == 	114966
 replace UCID_1934_2003 = 	223022	 if UCID_1934_2003 == 	115178
 replace UCID_1934_2003 = 	305776	 if UCID_1934_2003 == 	124988
 replace UCID_1934_2003 = 	153095	 if UCID_1934_2003 == 	166500
 replace UCID_1934_2003 = 	193785	 if UCID_1934_2003 == 	190767
 replace UCID_1934_2003 = 	292046	 if UCID_1934_2003 == 	198342
 replace UCID_1934_2003 = 	7376	 if UCID_1934_2003 == 	202214
 replace UCID_1934_2003 = 	179072	 if UCID_1934_2003 == 	213308
 replace UCID_1934_2003 = 	224096	 if UCID_1934_2003 == 	213535
 replace UCID_1934_2003 = 	270903	 if UCID_1934_2003 == 	217013
 replace UCID_1934_2003 = 	82325	 if UCID_1934_2003 == 	217684
 replace UCID_1934_2003 = 	363911	 if UCID_1934_2003 == 	250083
 replace UCID_1934_2003 = 	161861	 if UCID_1934_2003 == 	262505
 replace UCID_1934_2003 = 	428024	 if UCID_1934_2003 == 	267760
 replace UCID_1934_2003 = 	52777	 if UCID_1934_2003 == 	268204
 replace UCID_1934_2003 = 	271777	 if UCID_1934_2003 == 	272270
 replace UCID_1934_2003 = 	292048	 if UCID_1934_2003 == 	293171
 replace UCID_1934_2003 = 	178148	 if UCID_1934_2003 == 	302953
 replace UCID_1934_2003 = 	320030	 if UCID_1934_2003 == 	315800
 replace UCID_1934_2003 = 	341312	 if UCID_1934_2003 == 	332051
 replace UCID_1934_2003 = 	409439	 if UCID_1934_2003 == 	361188
 replace UCID_1934_2003 = 	413042	 if UCID_1934_2003 == 	361921
 replace UCID_1934_2003 = 	326987	 if UCID_1934_2003 == 	379614
 replace UCID_1934_2003 = 	297964	 if UCID_1934_2003 == 	393713
 replace UCID_1934_2003 = 	395418	 if UCID_1934_2003 == 	405646
 replace UCID_1934_2003 = 	15242	 if UCID_1934_2003 == 	414941
 replace UCID_1934_2003 = 	47512	 if UCID_1934_2003 == 	425788
 replace UCID_1934_2003 = 	432184	 if UCID_1934_2003 == 	430252

 * Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_47.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_46.dta"

*******************************************************************************
* 1972 & 1975

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_47.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1972 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1972
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1972
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1972_remaining_2.dta", replace
	restore

*Second year 1975 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1975
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1975
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1975_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1972 & 1975
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1972_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1975_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1972_1975_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1972_1975_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1972_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1975_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1972_1975_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1972_1975_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1972_1975_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1972 to 1975
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_47.dta", clear

 replace UCID_1934_2003 = 	171357	 if UCID_1934_2003 == 	2596
 replace UCID_1934_2003 = 	309	 if UCID_1934_2003 == 	8203
 replace UCID_1934_2003 = 	316977	 if UCID_1934_2003 == 	9930
 replace UCID_1934_2003 = 	40240	 if UCID_1934_2003 == 	10106
 replace UCID_1934_2003 = 	11913	 if UCID_1934_2003 == 	11176
 replace UCID_1934_2003 = 	12954	 if UCID_1934_2003 == 	13131
 replace UCID_1934_2003 = 	137328	 if UCID_1934_2003 == 	13400
 replace UCID_1934_2003 = 	23585	 if UCID_1934_2003 == 	23573
 replace UCID_1934_2003 = 	47537	 if UCID_1934_2003 == 	23972
 replace UCID_1934_2003 = 	268354	 if UCID_1934_2003 == 	32517
 replace UCID_1934_2003 = 	21880	 if UCID_1934_2003 == 	40725
 replace UCID_1934_2003 = 	47684	 if UCID_1934_2003 == 	47682
 replace UCID_1934_2003 = 	21253	 if UCID_1934_2003 == 	48057
 replace UCID_1934_2003 = 	72799	 if UCID_1934_2003 == 	62030
 replace UCID_1934_2003 = 	64790	 if UCID_1934_2003 == 	73906
 replace UCID_1934_2003 = 	72338	 if UCID_1934_2003 == 	88228
 replace UCID_1934_2003 = 	169449	 if UCID_1934_2003 == 	95100
 replace UCID_1934_2003 = 	119511	 if UCID_1934_2003 == 	111959
 replace UCID_1934_2003 = 	52184	 if UCID_1934_2003 == 	140695
 replace UCID_1934_2003 = 	14256	 if UCID_1934_2003 == 	150446
 replace UCID_1934_2003 = 	75458	 if UCID_1934_2003 == 	197028
 replace UCID_1934_2003 = 	190481	 if UCID_1934_2003 == 	203296
 replace UCID_1934_2003 = 	207349	 if UCID_1934_2003 == 	206184
 replace UCID_1934_2003 = 	208740	 if UCID_1934_2003 == 	207361
 replace UCID_1934_2003 = 	28041	 if UCID_1934_2003 == 	234356
 replace UCID_1934_2003 = 	294372	 if UCID_1934_2003 == 	256201
 replace UCID_1934_2003 = 	262685	 if UCID_1934_2003 == 	264738
 replace UCID_1934_2003 = 	274972	 if UCID_1934_2003 == 	275070
 replace UCID_1934_2003 = 	182966	 if UCID_1934_2003 == 	275624
 replace UCID_1934_2003 = 	388051	 if UCID_1934_2003 == 	302226
 replace UCID_1934_2003 = 	221847	 if UCID_1934_2003 == 	325427
 replace UCID_1934_2003 = 	83046	 if UCID_1934_2003 == 	325861
 replace UCID_1934_2003 = 	214748	 if UCID_1934_2003 == 	326247
 replace UCID_1934_2003 = 	338907	 if UCID_1934_2003 == 	329575
 replace UCID_1934_2003 = 	200689	 if UCID_1934_2003 == 	337691
 replace UCID_1934_2003 = 	327299	 if UCID_1934_2003 == 	339491
 replace UCID_1934_2003 = 	192907	 if UCID_1934_2003 == 	377108
 replace UCID_1934_2003 = 	3450	 if UCID_1934_2003 == 	425671
 replace UCID_1934_2003 = 	47632	 if UCID_1934_2003 == 	426311

 * Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_48.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_47.dta"

*******************************************************************************
* 1975 & 1979

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_48.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1975 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1975
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1975
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1975_remaining_2.dta", replace
	restore

*Second year 1979 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1979
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1979
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1979_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1975 & 1979
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1975_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1979_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1975_1979_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1975_1979_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1975_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1979_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1975_1979_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1975_1979_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1975_1979_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1975 to 1979
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_48.dta", clear

 replace UCID_1934_2003 = 	354371	 if UCID_1934_2003 == 	2120
 replace UCID_1934_2003 = 	236821	 if UCID_1934_2003 == 	2941
 replace UCID_1934_2003 = 	1023	 if UCID_1934_2003 == 	2996
 replace UCID_1934_2003 = 	61490	 if UCID_1934_2003 == 	3005
 replace UCID_1934_2003 = 	180728	 if UCID_1934_2003 == 	3544
 replace UCID_1934_2003 = 	65869	 if UCID_1934_2003 == 	3880
 replace UCID_1934_2003 = 	44089	 if UCID_1934_2003 == 	4652
 replace UCID_1934_2003 = 	67153	 if UCID_1934_2003 == 	4813
 replace UCID_1934_2003 = 	319818	 if UCID_1934_2003 == 	5663
 replace UCID_1934_2003 = 	171674	 if UCID_1934_2003 == 	6542
 replace UCID_1934_2003 = 	576	 if UCID_1934_2003 == 	7221
 replace UCID_1934_2003 = 	63800	 if UCID_1934_2003 == 	7595
 replace UCID_1934_2003 = 	412189	 if UCID_1934_2003 == 	9022
 replace UCID_1934_2003 = 	279309	 if UCID_1934_2003 == 	9351
 replace UCID_1934_2003 = 	171772	 if UCID_1934_2003 == 	11330
 replace UCID_1934_2003 = 	322017	 if UCID_1934_2003 == 	11539
 replace UCID_1934_2003 = 	412640	 if UCID_1934_2003 == 	12010
 replace UCID_1934_2003 = 	144371	 if UCID_1934_2003 == 	12920
 replace UCID_1934_2003 = 	142387	 if UCID_1934_2003 == 	13754
 replace UCID_1934_2003 = 	344796	 if UCID_1934_2003 == 	15806
 replace UCID_1934_2003 = 	164689	 if UCID_1934_2003 == 	17088
 replace UCID_1934_2003 = 	162295	 if UCID_1934_2003 == 	17900
 replace UCID_1934_2003 = 	374867	 if UCID_1934_2003 == 	20026
 replace UCID_1934_2003 = 	200408	 if UCID_1934_2003 == 	21974
 replace UCID_1934_2003 = 	411877	 if UCID_1934_2003 == 	23637
 replace UCID_1934_2003 = 	228518	 if UCID_1934_2003 == 	26840
 replace UCID_1934_2003 = 	240105	 if UCID_1934_2003 == 	28118
 replace UCID_1934_2003 = 	244161	 if UCID_1934_2003 == 	28444
 replace UCID_1934_2003 = 	242918	 if UCID_1934_2003 == 	29212
 replace UCID_1934_2003 = 	84074	 if UCID_1934_2003 == 	29327
 replace UCID_1934_2003 = 	241128	 if UCID_1934_2003 == 	29478
 replace UCID_1934_2003 = 	29894	 if UCID_1934_2003 == 	29890
 replace UCID_1934_2003 = 	256095	 if UCID_1934_2003 == 	30461
 replace UCID_1934_2003 = 	183644	 if UCID_1934_2003 == 	30645
 replace UCID_1934_2003 = 	15979	 if UCID_1934_2003 == 	31292
 replace UCID_1934_2003 = 	413366	 if UCID_1934_2003 == 	32236
 replace UCID_1934_2003 = 	272220	 if UCID_1934_2003 == 	32725
 replace UCID_1934_2003 = 	279015	 if UCID_1934_2003 == 	32878
 replace UCID_1934_2003 = 	300446	 if UCID_1934_2003 == 	35262
 replace UCID_1934_2003 = 	96823	 if UCID_1934_2003 == 	36463
 replace UCID_1934_2003 = 	342742	 if UCID_1934_2003 == 	40299
 replace UCID_1934_2003 = 	40080	 if UCID_1934_2003 == 	41430
 replace UCID_1934_2003 = 	355825	 if UCID_1934_2003 == 	41683
 replace UCID_1934_2003 = 	356225	 if UCID_1934_2003 == 	41953
 replace UCID_1934_2003 = 	356381	 if UCID_1934_2003 == 	41966
 replace UCID_1934_2003 = 	356295	 if UCID_1934_2003 == 	42014
 replace UCID_1934_2003 = 	402983	 if UCID_1934_2003 == 	44594
 replace UCID_1934_2003 = 	155303	 if UCID_1934_2003 == 	47438
 replace UCID_1934_2003 = 	47915	 if UCID_1934_2003 == 	47451
 replace UCID_1934_2003 = 	21481	 if UCID_1934_2003 == 	47473
 replace UCID_1934_2003 = 	415437	 if UCID_1934_2003 == 	47475
 replace UCID_1934_2003 = 	156080	 if UCID_1934_2003 == 	47476
 replace UCID_1934_2003 = 	410290	 if UCID_1934_2003 == 	47482
 replace UCID_1934_2003 = 	409946	 if UCID_1934_2003 == 	47582
 replace UCID_1934_2003 = 	429658	 if UCID_1934_2003 == 	47619
 replace UCID_1934_2003 = 	273797	 if UCID_1934_2003 == 	47620
 replace UCID_1934_2003 = 	47778	 if UCID_1934_2003 == 	47628
 replace UCID_1934_2003 = 	11819	 if UCID_1934_2003 == 	47630
 replace UCID_1934_2003 = 	414058	 if UCID_1934_2003 == 	47640
 replace UCID_1934_2003 = 	415257	 if UCID_1934_2003 == 	47714
 replace UCID_1934_2003 = 	415902	 if UCID_1934_2003 == 	47715
 replace UCID_1934_2003 = 	30095	 if UCID_1934_2003 == 	47721
 replace UCID_1934_2003 = 	409385	 if UCID_1934_2003 == 	47753
 replace UCID_1934_2003 = 	410355	 if UCID_1934_2003 == 	47831
 replace UCID_1934_2003 = 	411881	 if UCID_1934_2003 == 	47832
 replace UCID_1934_2003 = 	410020	 if UCID_1934_2003 == 	47833
 replace UCID_1934_2003 = 	414744	 if UCID_1934_2003 == 	47897
 replace UCID_1934_2003 = 	409960	 if UCID_1934_2003 == 	47938
 replace UCID_1934_2003 = 	53184	 if UCID_1934_2003 == 	48461
 replace UCID_1934_2003 = 	70737	 if UCID_1934_2003 == 	50127
 replace UCID_1934_2003 = 	53890	 if UCID_1934_2003 == 	50716
 replace UCID_1934_2003 = 	181543	 if UCID_1934_2003 == 	51231
 replace UCID_1934_2003 = 	259773	 if UCID_1934_2003 == 	52071
 replace UCID_1934_2003 = 	104643	 if UCID_1934_2003 == 	53971
 replace UCID_1934_2003 = 	90092	 if UCID_1934_2003 == 	54524
 replace UCID_1934_2003 = 	62448	 if UCID_1934_2003 == 	55058
 replace UCID_1934_2003 = 	415756	 if UCID_1934_2003 == 	57655
 replace UCID_1934_2003 = 	132058	 if UCID_1934_2003 == 	57979
 replace UCID_1934_2003 = 	181445	 if UCID_1934_2003 == 	58084
 replace UCID_1934_2003 = 	55443	 if UCID_1934_2003 == 	64537
 replace UCID_1934_2003 = 	93552	 if UCID_1934_2003 == 	65486
 replace UCID_1934_2003 = 	56476	 if UCID_1934_2003 == 	66332
 replace UCID_1934_2003 = 	60996	 if UCID_1934_2003 == 	66485
 replace UCID_1934_2003 = 	65229	 if UCID_1934_2003 == 	66646
 replace UCID_1934_2003 = 	63432	 if UCID_1934_2003 == 	70303
 replace UCID_1934_2003 = 	67741	 if UCID_1934_2003 == 	70309
 replace UCID_1934_2003 = 	269665	 if UCID_1934_2003 == 	70598
 replace UCID_1934_2003 = 	269647	 if UCID_1934_2003 == 	70631
 replace UCID_1934_2003 = 	272829	 if UCID_1934_2003 == 	71309
 replace UCID_1934_2003 = 	414684	 if UCID_1934_2003 == 	71544
 replace UCID_1934_2003 = 	1258	 if UCID_1934_2003 == 	72319
 replace UCID_1934_2003 = 	61277	 if UCID_1934_2003 == 	73757
 replace UCID_1934_2003 = 	158172	 if UCID_1934_2003 == 	74189
 replace UCID_1934_2003 = 	295440	 if UCID_1934_2003 == 	75095
 replace UCID_1934_2003 = 	214332	 if UCID_1934_2003 == 	76459
 replace UCID_1934_2003 = 	412164	 if UCID_1934_2003 == 	77827
 replace UCID_1934_2003 = 	410258	 if UCID_1934_2003 == 	78371
 replace UCID_1934_2003 = 	249099	 if UCID_1934_2003 == 	79034
 replace UCID_1934_2003 = 	56263	 if UCID_1934_2003 == 	79411
 replace UCID_1934_2003 = 	283975	 if UCID_1934_2003 == 	79565
 replace UCID_1934_2003 = 	415842	 if UCID_1934_2003 == 	81143
 replace UCID_1934_2003 = 	410498	 if UCID_1934_2003 == 	81712
 replace UCID_1934_2003 = 	62681	 if UCID_1934_2003 == 	83464
 replace UCID_1934_2003 = 	2382	 if UCID_1934_2003 == 	83471
 replace UCID_1934_2003 = 	72265	 if UCID_1934_2003 == 	84679
 replace UCID_1934_2003 = 	413965	 if UCID_1934_2003 == 	84750
 replace UCID_1934_2003 = 	300132	 if UCID_1934_2003 == 	87249
 replace UCID_1934_2003 = 	52940	 if UCID_1934_2003 == 	87353
 replace UCID_1934_2003 = 	412404	 if UCID_1934_2003 == 	88697
 replace UCID_1934_2003 = 	53052	 if UCID_1934_2003 == 	90964
 replace UCID_1934_2003 = 	57421	 if UCID_1934_2003 == 	91803
 replace UCID_1934_2003 = 	234680	 if UCID_1934_2003 == 	92714
 replace UCID_1934_2003 = 	58508	 if UCID_1934_2003 == 	93155
 replace UCID_1934_2003 = 	72234	 if UCID_1934_2003 == 	96259
 replace UCID_1934_2003 = 	50693	 if UCID_1934_2003 == 	97185
 replace UCID_1934_2003 = 	57281	 if UCID_1934_2003 == 	99995
 replace UCID_1934_2003 = 	250494	 if UCID_1934_2003 == 	101534
 replace UCID_1934_2003 = 	15280	 if UCID_1934_2003 == 	101541
 replace UCID_1934_2003 = 	62280	 if UCID_1934_2003 == 	103519
 replace UCID_1934_2003 = 	414936	 if UCID_1934_2003 == 	104015
 replace UCID_1934_2003 = 	83037	 if UCID_1934_2003 == 	104645
 replace UCID_1934_2003 = 	69923	 if UCID_1934_2003 == 	104744
 replace UCID_1934_2003 = 	107968	 if UCID_1934_2003 == 	104975
 replace UCID_1934_2003 = 	58615	 if UCID_1934_2003 == 	105142
 replace UCID_1934_2003 = 	50120	 if UCID_1934_2003 == 	106421
 replace UCID_1934_2003 = 	413823	 if UCID_1934_2003 == 	108220
 replace UCID_1934_2003 = 	108452	 if UCID_1934_2003 == 	108562
 replace UCID_1934_2003 = 	211478	 if UCID_1934_2003 == 	111726
 replace UCID_1934_2003 = 	131431	 if UCID_1934_2003 == 	114617
 replace UCID_1934_2003 = 	128407	 if UCID_1934_2003 == 	117452
 replace UCID_1934_2003 = 	411066	 if UCID_1934_2003 == 	118175
 replace UCID_1934_2003 = 	118718	 if UCID_1934_2003 == 	118731
 replace UCID_1934_2003 = 	118736	 if UCID_1934_2003 == 	118747
 replace UCID_1934_2003 = 	129372	 if UCID_1934_2003 == 	120951
 replace UCID_1934_2003 = 	129943	 if UCID_1934_2003 == 	121144
 replace UCID_1934_2003 = 	364330	 if UCID_1934_2003 == 	122068
 replace UCID_1934_2003 = 	190143	 if UCID_1934_2003 == 	122800
 replace UCID_1934_2003 = 	414469	 if UCID_1934_2003 == 	124005
 replace UCID_1934_2003 = 	127723	 if UCID_1934_2003 == 	127728
 replace UCID_1934_2003 = 	415935	 if UCID_1934_2003 == 	129459
 replace UCID_1934_2003 = 	120073	 if UCID_1934_2003 == 	129869
 replace UCID_1934_2003 = 	115822	 if UCID_1934_2003 == 	131632
 replace UCID_1934_2003 = 	111301	 if UCID_1934_2003 == 	132080
 replace UCID_1934_2003 = 	133116	 if UCID_1934_2003 == 	132887
 replace UCID_1934_2003 = 	134423	 if UCID_1934_2003 == 	134462
 replace UCID_1934_2003 = 	123759	 if UCID_1934_2003 == 	135589
 replace UCID_1934_2003 = 	137118	 if UCID_1934_2003 == 	137137
 replace UCID_1934_2003 = 	135966	 if UCID_1934_2003 == 	137484
 replace UCID_1934_2003 = 	144005	 if UCID_1934_2003 == 	138099
 replace UCID_1934_2003 = 	137891	 if UCID_1934_2003 == 	138372
 replace UCID_1934_2003 = 	137413	 if UCID_1934_2003 == 	140704
 replace UCID_1934_2003 = 	139362	 if UCID_1934_2003 == 	144116
 replace UCID_1934_2003 = 	149101	 if UCID_1934_2003 == 	148897
 replace UCID_1934_2003 = 	204956	 if UCID_1934_2003 == 	161709
 replace UCID_1934_2003 = 	184494	 if UCID_1934_2003 == 	165034
 replace UCID_1934_2003 = 	164755	 if UCID_1934_2003 == 	175854
 replace UCID_1934_2003 = 	169354	 if UCID_1934_2003 == 	176410
 replace UCID_1934_2003 = 	163158	 if UCID_1934_2003 == 	184242
 replace UCID_1934_2003 = 	189898	 if UCID_1934_2003 == 	190182
 replace UCID_1934_2003 = 	182736	 if UCID_1934_2003 == 	194814
 replace UCID_1934_2003 = 	199641	 if UCID_1934_2003 == 	196321
 replace UCID_1934_2003 = 	411562	 if UCID_1934_2003 == 	204493
 replace UCID_1934_2003 = 	172828	 if UCID_1934_2003 == 	205837
 replace UCID_1934_2003 = 	206274	 if UCID_1934_2003 == 	209215
 replace UCID_1934_2003 = 	25846	 if UCID_1934_2003 == 	210604
 replace UCID_1934_2003 = 	415172	 if UCID_1934_2003 == 	210708
 replace UCID_1934_2003 = 	5119	 if UCID_1934_2003 == 	214070
 replace UCID_1934_2003 = 	427379	 if UCID_1934_2003 == 	215375
 replace UCID_1934_2003 = 	33117	 if UCID_1934_2003 == 	215552
 replace UCID_1934_2003 = 	109525	 if UCID_1934_2003 == 	217161
 replace UCID_1934_2003 = 	240729	 if UCID_1934_2003 == 	219110
 replace UCID_1934_2003 = 	214180	 if UCID_1934_2003 == 	223058
 replace UCID_1934_2003 = 	416352	 if UCID_1934_2003 == 	223356
 replace UCID_1934_2003 = 	229754	 if UCID_1934_2003 == 	228664
 replace UCID_1934_2003 = 	135214	 if UCID_1934_2003 == 	230116
 replace UCID_1934_2003 = 	229775	 if UCID_1934_2003 == 	230532
 replace UCID_1934_2003 = 	33479	 if UCID_1934_2003 == 	235766
 replace UCID_1934_2003 = 	158454	 if UCID_1934_2003 == 	237532
 replace UCID_1934_2003 = 	252287	 if UCID_1934_2003 == 	238917
 replace UCID_1934_2003 = 	87619	 if UCID_1934_2003 == 	241437
 replace UCID_1934_2003 = 	242295	 if UCID_1934_2003 == 	242615
 replace UCID_1934_2003 = 	411924	 if UCID_1934_2003 == 	242883
 replace UCID_1934_2003 = 	243015	 if UCID_1934_2003 == 	242923
 replace UCID_1934_2003 = 	242805	 if UCID_1934_2003 == 	243197
 replace UCID_1934_2003 = 	243001	 if UCID_1934_2003 == 	243224
 replace UCID_1934_2003 = 	242961	 if UCID_1934_2003 == 	243238
 replace UCID_1934_2003 = 	242891	 if UCID_1934_2003 == 	243315
 replace UCID_1934_2003 = 	244063	 if UCID_1934_2003 == 	243747
 replace UCID_1934_2003 = 	244449	 if UCID_1934_2003 == 	244665
 replace UCID_1934_2003 = 	243091	 if UCID_1934_2003 == 	245263
 replace UCID_1934_2003 = 	245311	 if UCID_1934_2003 == 	245302
 replace UCID_1934_2003 = 	58361	 if UCID_1934_2003 == 	245641
 replace UCID_1934_2003 = 	411596	 if UCID_1934_2003 == 	245719
 replace UCID_1934_2003 = 	415750	 if UCID_1934_2003 == 	245772
 replace UCID_1934_2003 = 	410698	 if UCID_1934_2003 == 	248715
 replace UCID_1934_2003 = 	249500	 if UCID_1934_2003 == 	249427
 replace UCID_1934_2003 = 	263676	 if UCID_1934_2003 == 	252844
 replace UCID_1934_2003 = 	263702	 if UCID_1934_2003 == 	254005
 replace UCID_1934_2003 = 	263092	 if UCID_1934_2003 == 	254309
 replace UCID_1934_2003 = 	47892	 if UCID_1934_2003 == 	258876
 replace UCID_1934_2003 = 	259475	 if UCID_1934_2003 == 	259502
 replace UCID_1934_2003 = 	171965	 if UCID_1934_2003 == 	260588
 replace UCID_1934_2003 = 	254193	 if UCID_1934_2003 == 	261059
 replace UCID_1934_2003 = 	31689	 if UCID_1934_2003 == 	264600
 replace UCID_1934_2003 = 	51941	 if UCID_1934_2003 == 	265919
 replace UCID_1934_2003 = 	267238	 if UCID_1934_2003 == 	267646
 replace UCID_1934_2003 = 	74554	 if UCID_1934_2003 == 	268695
 replace UCID_1934_2003 = 	70603	 if UCID_1934_2003 == 	269684
 replace UCID_1934_2003 = 	188962	 if UCID_1934_2003 == 	270385
 replace UCID_1934_2003 = 	71880	 if UCID_1934_2003 == 	271123
 replace UCID_1934_2003 = 	409729	 if UCID_1934_2003 == 	271849
 replace UCID_1934_2003 = 	151530	 if UCID_1934_2003 == 	272854
 replace UCID_1934_2003 = 	410507	 if UCID_1934_2003 == 	273084
 replace UCID_1934_2003 = 	250328	 if UCID_1934_2003 == 	277012
 replace UCID_1934_2003 = 	30051	 if UCID_1934_2003 == 	277671
 replace UCID_1934_2003 = 	277306	 if UCID_1934_2003 == 	277703
 replace UCID_1934_2003 = 	279267	 if UCID_1934_2003 == 	278940
 replace UCID_1934_2003 = 	270874	 if UCID_1934_2003 == 	280511
 replace UCID_1934_2003 = 	269835	 if UCID_1934_2003 == 	281345
 replace UCID_1934_2003 = 	147394	 if UCID_1934_2003 == 	281471
 replace UCID_1934_2003 = 	230142	 if UCID_1934_2003 == 	284716
 replace UCID_1934_2003 = 	412367	 if UCID_1934_2003 == 	287090
 replace UCID_1934_2003 = 	287960	 if UCID_1934_2003 == 	287917
 replace UCID_1934_2003 = 	223608	 if UCID_1934_2003 == 	289195
 replace UCID_1934_2003 = 	291440	 if UCID_1934_2003 == 	290463
 replace UCID_1934_2003 = 	294734	 if UCID_1934_2003 == 	294730
 replace UCID_1934_2003 = 	294629	 if UCID_1934_2003 == 	300754
 replace UCID_1934_2003 = 	210204	 if UCID_1934_2003 == 	301482
 replace UCID_1934_2003 = 	293740	 if UCID_1934_2003 == 	301615
 replace UCID_1934_2003 = 	111015	 if UCID_1934_2003 == 	305286
 replace UCID_1934_2003 = 	265136	 if UCID_1934_2003 == 	306930
 replace UCID_1934_2003 = 	1089	 if UCID_1934_2003 == 	306983
 replace UCID_1934_2003 = 	196141	 if UCID_1934_2003 == 	307058
 replace UCID_1934_2003 = 	182765	 if UCID_1934_2003 == 	307331
 replace UCID_1934_2003 = 	293542	 if UCID_1934_2003 == 	320004
 replace UCID_1934_2003 = 	323688	 if UCID_1934_2003 == 	323462
 replace UCID_1934_2003 = 	323707	 if UCID_1934_2003 == 	323496
 replace UCID_1934_2003 = 	339791	 if UCID_1934_2003 == 	325551
 replace UCID_1934_2003 = 	347385	 if UCID_1934_2003 == 	326448
 replace UCID_1934_2003 = 	339201	 if UCID_1934_2003 == 	326735
 replace UCID_1934_2003 = 	340020	 if UCID_1934_2003 == 	327112
 replace UCID_1934_2003 = 	365282	 if UCID_1934_2003 == 	327556
 replace UCID_1934_2003 = 	413604	 if UCID_1934_2003 == 	327702
 replace UCID_1934_2003 = 	341951	 if UCID_1934_2003 == 	327938
 replace UCID_1934_2003 = 	413825	 if UCID_1934_2003 == 	330501
 replace UCID_1934_2003 = 	410719	 if UCID_1934_2003 == 	330663
 replace UCID_1934_2003 = 	411115	 if UCID_1934_2003 == 	332389
 replace UCID_1934_2003 = 	169908	 if UCID_1934_2003 == 	333193
 replace UCID_1934_2003 = 	408553	 if UCID_1934_2003 == 	333279
 replace UCID_1934_2003 = 	335886	 if UCID_1934_2003 == 	337411
 replace UCID_1934_2003 = 	342119	 if UCID_1934_2003 == 	342301
 replace UCID_1934_2003 = 	331420	 if UCID_1934_2003 == 	346433
 replace UCID_1934_2003 = 	346887	 if UCID_1934_2003 == 	346872
 replace UCID_1934_2003 = 	170255	 if UCID_1934_2003 == 	348136
 replace UCID_1934_2003 = 	346947	 if UCID_1934_2003 == 	349338
 replace UCID_1934_2003 = 	350098	 if UCID_1934_2003 == 	349835
 replace UCID_1934_2003 = 	194932	 if UCID_1934_2003 == 	350565
 replace UCID_1934_2003 = 	352670	 if UCID_1934_2003 == 	352742
 replace UCID_1934_2003 = 	253564	 if UCID_1934_2003 == 	352805
 replace UCID_1934_2003 = 	355683	 if UCID_1934_2003 == 	355660
 replace UCID_1934_2003 = 	356193	 if UCID_1934_2003 == 	356231
 replace UCID_1934_2003 = 	356389	 if UCID_1934_2003 == 	356250
 replace UCID_1934_2003 = 	419150	 if UCID_1934_2003 == 	356269
 replace UCID_1934_2003 = 	356194	 if UCID_1934_2003 == 	356337
 replace UCID_1934_2003 = 	356214	 if UCID_1934_2003 == 	356355
 replace UCID_1934_2003 = 	373166	 if UCID_1934_2003 == 	357120
 replace UCID_1934_2003 = 	356412	 if UCID_1934_2003 == 	357370
 replace UCID_1934_2003 = 	423695	 if UCID_1934_2003 == 	357404
 replace UCID_1934_2003 = 	356420	 if UCID_1934_2003 == 	357416
 replace UCID_1934_2003 = 	356178	 if UCID_1934_2003 == 	357461
 replace UCID_1934_2003 = 	356346	 if UCID_1934_2003 == 	357463
 replace UCID_1934_2003 = 	356245	 if UCID_1934_2003 == 	357471
 replace UCID_1934_2003 = 	356316	 if UCID_1934_2003 == 	357526
 replace UCID_1934_2003 = 	356274	 if UCID_1934_2003 == 	357541
 replace UCID_1934_2003 = 	356161	 if UCID_1934_2003 == 	357567
 replace UCID_1934_2003 = 	356130	 if UCID_1934_2003 == 	357616
 replace UCID_1934_2003 = 	356212	 if UCID_1934_2003 == 	357627
 replace UCID_1934_2003 = 	356219	 if UCID_1934_2003 == 	357657
 replace UCID_1934_2003 = 	356299	 if UCID_1934_2003 == 	357674
 replace UCID_1934_2003 = 	356144	 if UCID_1934_2003 == 	357733
 replace UCID_1934_2003 = 	356163	 if UCID_1934_2003 == 	357771
 replace UCID_1934_2003 = 	356210	 if UCID_1934_2003 == 	357814
 replace UCID_1934_2003 = 	356419	 if UCID_1934_2003 == 	357878
 replace UCID_1934_2003 = 	356228	 if UCID_1934_2003 == 	357902
 replace UCID_1934_2003 = 	356410	 if UCID_1934_2003 == 	357991
 replace UCID_1934_2003 = 	356423	 if UCID_1934_2003 == 	358014
 replace UCID_1934_2003 = 	356233	 if UCID_1934_2003 == 	358053
 replace UCID_1934_2003 = 	356279	 if UCID_1934_2003 == 	358153
 replace UCID_1934_2003 = 	423230	 if UCID_1934_2003 == 	358156
 replace UCID_1934_2003 = 	356198	 if UCID_1934_2003 == 	358267
 replace UCID_1934_2003 = 	356302	 if UCID_1934_2003 == 	358361
 replace UCID_1934_2003 = 	351750	 if UCID_1934_2003 == 	359379
 replace UCID_1934_2003 = 	352815	 if UCID_1934_2003 == 	360050
 replace UCID_1934_2003 = 	361508	 if UCID_1934_2003 == 	364884
 replace UCID_1934_2003 = 	338305	 if UCID_1934_2003 == 	369976
 replace UCID_1934_2003 = 	412988	 if UCID_1934_2003 == 	370617
 replace UCID_1934_2003 = 	372192	 if UCID_1934_2003 == 	372784
 replace UCID_1934_2003 = 	215605	 if UCID_1934_2003 == 	373191
 replace UCID_1934_2003 = 	391631	 if UCID_1934_2003 == 	375576
 replace UCID_1934_2003 = 	379122	 if UCID_1934_2003 == 	376081
 replace UCID_1934_2003 = 	170980	 if UCID_1934_2003 == 	376807
 replace UCID_1934_2003 = 	78857	 if UCID_1934_2003 == 	380178
 replace UCID_1934_2003 = 	402898	 if UCID_1934_2003 == 	384735
 replace UCID_1934_2003 = 	418942	 if UCID_1934_2003 == 	386239
 replace UCID_1934_2003 = 	368552	 if UCID_1934_2003 == 	387208
 replace UCID_1934_2003 = 	192034	 if UCID_1934_2003 == 	389088
 replace UCID_1934_2003 = 	174959	 if UCID_1934_2003 == 	395013
 replace UCID_1934_2003 = 	404898	 if UCID_1934_2003 == 	396306
 replace UCID_1934_2003 = 	398117	 if UCID_1934_2003 == 	397306
 replace UCID_1934_2003 = 	41828	 if UCID_1934_2003 == 	398964
 replace UCID_1934_2003 = 	383851	 if UCID_1934_2003 == 	401676
 replace UCID_1934_2003 = 	330482	 if UCID_1934_2003 == 	409068
 replace UCID_1934_2003 = 	411243	 if UCID_1934_2003 == 	409078
 replace UCID_1934_2003 = 	201225	 if UCID_1934_2003 == 	409095
 replace UCID_1934_2003 = 	110068	 if UCID_1934_2003 == 	409099
 replace UCID_1934_2003 = 	412090	 if UCID_1934_2003 == 	409146
 replace UCID_1934_2003 = 	47692	 if UCID_1934_2003 == 	409161
 replace UCID_1934_2003 = 	409690	 if UCID_1934_2003 == 	409245
 replace UCID_1934_2003 = 	32432	 if UCID_1934_2003 == 	409247
 replace UCID_1934_2003 = 	414010	 if UCID_1934_2003 == 	409274
 replace UCID_1934_2003 = 	410277	 if UCID_1934_2003 == 	409294
 replace UCID_1934_2003 = 	411329	 if UCID_1934_2003 == 	409296
 replace UCID_1934_2003 = 	409288	 if UCID_1934_2003 == 	409328
 replace UCID_1934_2003 = 	411374	 if UCID_1934_2003 == 	409333
 replace UCID_1934_2003 = 	409067	 if UCID_1934_2003 == 	409344
 replace UCID_1934_2003 = 	409117	 if UCID_1934_2003 == 	409459
 replace UCID_1934_2003 = 	414792	 if UCID_1934_2003 == 	409475
 replace UCID_1934_2003 = 	414464	 if UCID_1934_2003 == 	409485
 replace UCID_1934_2003 = 	410819	 if UCID_1934_2003 == 	409504
 replace UCID_1934_2003 = 	413081	 if UCID_1934_2003 == 	409669
 replace UCID_1934_2003 = 	415401	 if UCID_1934_2003 == 	409682
 replace UCID_1934_2003 = 	412732	 if UCID_1934_2003 == 	409712
 replace UCID_1934_2003 = 	409088	 if UCID_1934_2003 == 	409720
 replace UCID_1934_2003 = 	322871	 if UCID_1934_2003 == 	409731
 replace UCID_1934_2003 = 	267083	 if UCID_1934_2003 == 	409750
 replace UCID_1934_2003 = 	126006	 if UCID_1934_2003 == 	409816
 replace UCID_1934_2003 = 	411511	 if UCID_1934_2003 == 	409976
 replace UCID_1934_2003 = 	410448	 if UCID_1934_2003 == 	409978
 replace UCID_1934_2003 = 	410922	 if UCID_1934_2003 == 	409982
 replace UCID_1934_2003 = 	415243	 if UCID_1934_2003 == 	410002
 replace UCID_1934_2003 = 	135383	 if UCID_1934_2003 == 	410081
 replace UCID_1934_2003 = 	410095	 if UCID_1934_2003 == 	410114
 replace UCID_1934_2003 = 	252401	 if UCID_1934_2003 == 	410159
 replace UCID_1934_2003 = 	409870	 if UCID_1934_2003 == 	410296
 replace UCID_1934_2003 = 	414254	 if UCID_1934_2003 == 	410393
 replace UCID_1934_2003 = 	48000	 if UCID_1934_2003 == 	410450
 replace UCID_1934_2003 = 	322582	 if UCID_1934_2003 == 	410573
 replace UCID_1934_2003 = 	158781	 if UCID_1934_2003 == 	410595
 replace UCID_1934_2003 = 	48096	 if UCID_1934_2003 == 	410685
 replace UCID_1934_2003 = 	32688	 if UCID_1934_2003 == 	410808
 replace UCID_1934_2003 = 	416005	 if UCID_1934_2003 == 	410980
 replace UCID_1934_2003 = 	409980	 if UCID_1934_2003 == 	411061
 replace UCID_1934_2003 = 	409239	 if UCID_1934_2003 == 	411109
 replace UCID_1934_2003 = 	411223	 if UCID_1934_2003 == 	411266
 replace UCID_1934_2003 = 	415365	 if UCID_1934_2003 == 	411293
 replace UCID_1934_2003 = 	273691	 if UCID_1934_2003 == 	411394
 replace UCID_1934_2003 = 	410461	 if UCID_1934_2003 == 	411468
 replace UCID_1934_2003 = 	416433	 if UCID_1934_2003 == 	411482
 replace UCID_1934_2003 = 	32430	 if UCID_1934_2003 == 	411594
 replace UCID_1934_2003 = 	414360	 if UCID_1934_2003 == 	411655
 replace UCID_1934_2003 = 	411029	 if UCID_1934_2003 == 	411903
 replace UCID_1934_2003 = 	410460	 if UCID_1934_2003 == 	411907
 replace UCID_1934_2003 = 	201246	 if UCID_1934_2003 == 	412080
 replace UCID_1934_2003 = 	252011	 if UCID_1934_2003 == 	412180
 replace UCID_1934_2003 = 	245807	 if UCID_1934_2003 == 	412241
 replace UCID_1934_2003 = 	56188	 if UCID_1934_2003 == 	412284
 replace UCID_1934_2003 = 	415006	 if UCID_1934_2003 == 	412327
 replace UCID_1934_2003 = 	415288	 if UCID_1934_2003 == 	412382
 replace UCID_1934_2003 = 	227392	 if UCID_1934_2003 == 	412427
 replace UCID_1934_2003 = 	409238	 if UCID_1934_2003 == 	412447
 replace UCID_1934_2003 = 	189808	 if UCID_1934_2003 == 	412460
 replace UCID_1934_2003 = 	409354	 if UCID_1934_2003 == 	412706
 replace UCID_1934_2003 = 	55797	 if UCID_1934_2003 == 	412713
 replace UCID_1934_2003 = 	414538	 if UCID_1934_2003 == 	412856
 replace UCID_1934_2003 = 	410645	 if UCID_1934_2003 == 	412892
 replace UCID_1934_2003 = 	411051	 if UCID_1934_2003 == 	412904
 replace UCID_1934_2003 = 	275438	 if UCID_1934_2003 == 	412909
 replace UCID_1934_2003 = 	414046	 if UCID_1934_2003 == 	412912
 replace UCID_1934_2003 = 	409680	 if UCID_1934_2003 == 	413068
 replace UCID_1934_2003 = 	411859	 if UCID_1934_2003 == 	413232
 replace UCID_1934_2003 = 	409491	 if UCID_1934_2003 == 	413280
 replace UCID_1934_2003 = 	248795	 if UCID_1934_2003 == 	413656
 replace UCID_1934_2003 = 	413700	 if UCID_1934_2003 == 	413776
 replace UCID_1934_2003 = 	29457	 if UCID_1934_2003 == 	413847
 replace UCID_1934_2003 = 	410334	 if UCID_1934_2003 == 	413898
 replace UCID_1934_2003 = 	281598	 if UCID_1934_2003 == 	413960
 replace UCID_1934_2003 = 	201228	 if UCID_1934_2003 == 	414179
 replace UCID_1934_2003 = 	415940	 if UCID_1934_2003 == 	414239
 replace UCID_1934_2003 = 	413017	 if UCID_1934_2003 == 	414285
 replace UCID_1934_2003 = 	273732	 if UCID_1934_2003 == 	414313
 replace UCID_1934_2003 = 	415805	 if UCID_1934_2003 == 	414366
 replace UCID_1934_2003 = 	410851	 if UCID_1934_2003 == 	414371
 replace UCID_1934_2003 = 	412568	 if UCID_1934_2003 == 	414462
 replace UCID_1934_2003 = 	414483	 if UCID_1934_2003 == 	414473
 replace UCID_1934_2003 = 	414900	 if UCID_1934_2003 == 	414645
 replace UCID_1934_2003 = 	125684	 if UCID_1934_2003 == 	414682
 replace UCID_1934_2003 = 	412294	 if UCID_1934_2003 == 	414707
 replace UCID_1934_2003 = 	413571	 if UCID_1934_2003 == 	414722
 replace UCID_1934_2003 = 	413240	 if UCID_1934_2003 == 	414967
 replace UCID_1934_2003 = 	409524	 if UCID_1934_2003 == 	415014
 replace UCID_1934_2003 = 	270707	 if UCID_1934_2003 == 	415055
 replace UCID_1934_2003 = 	412021	 if UCID_1934_2003 == 	415055
 replace UCID_1934_2003 = 	224015	 if UCID_1934_2003 == 	415087
 replace UCID_1934_2003 = 	410013	 if UCID_1934_2003 == 	415107
 replace UCID_1934_2003 = 	320607	 if UCID_1934_2003 == 	415128
 replace UCID_1934_2003 = 	409791	 if UCID_1934_2003 == 	415137
 replace UCID_1934_2003 = 	412739	 if UCID_1934_2003 == 	415169
 replace UCID_1934_2003 = 	252639	 if UCID_1934_2003 == 	415253
 replace UCID_1934_2003 = 	273714	 if UCID_1934_2003 == 	415408
 replace UCID_1934_2003 = 	34090	 if UCID_1934_2003 == 	415445
 replace UCID_1934_2003 = 	255873	 if UCID_1934_2003 == 	415535
 replace UCID_1934_2003 = 	409292	 if UCID_1934_2003 == 	415542
 replace UCID_1934_2003 = 	349690	 if UCID_1934_2003 == 	415625
 replace UCID_1934_2003 = 	199948	 if UCID_1934_2003 == 	415657
 replace UCID_1934_2003 = 	133663	 if UCID_1934_2003 == 	415784
 replace UCID_1934_2003 = 	413092	 if UCID_1934_2003 == 	415897
 replace UCID_1934_2003 = 	47573	 if UCID_1934_2003 == 	415920
 replace UCID_1934_2003 = 	275779	 if UCID_1934_2003 == 	415973
 replace UCID_1934_2003 = 	412566	 if UCID_1934_2003 == 	416038
 replace UCID_1934_2003 = 	27609	 if UCID_1934_2003 == 	416079
 replace UCID_1934_2003 = 	251834	 if UCID_1934_2003 == 	416082
 replace UCID_1934_2003 = 	288556	 if UCID_1934_2003 == 	416157
 replace UCID_1934_2003 = 	251821	 if UCID_1934_2003 == 	416207
 replace UCID_1934_2003 = 	410597	 if UCID_1934_2003 == 	416283
 replace UCID_1934_2003 = 	410993	 if UCID_1934_2003 == 	416358
 replace UCID_1934_2003 = 	27340	 if UCID_1934_2003 == 	416388
 replace UCID_1934_2003 = 	411529	 if UCID_1934_2003 == 	416406
 replace UCID_1934_2003 = 	356142	 if UCID_1934_2003 == 	417951
 replace UCID_1934_2003 = 	423696	 if UCID_1934_2003 == 	421285
 replace UCID_1934_2003 = 	418926	 if UCID_1934_2003 == 	421289
 replace UCID_1934_2003 = 	153730	 if UCID_1934_2003 == 	426399
 replace UCID_1934_2003 = 	37610	 if UCID_1934_2003 == 	428364
 replace UCID_1934_2003 = 	411292	 if UCID_1934_2003 == 	429419
 replace UCID_1934_2003 = 	361999	 if UCID_1934_2003 == 	429430

 * Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_49.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_48.dta"

*******************************************************************************
* 1979 & 1980

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_49.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1979 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1979
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1979
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1979_remaining_2.dta", replace
	restore

*Second year 1980 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1980
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1980
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1980_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1979 & 1980
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1979_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1980_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1979_1980_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1979_1980_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1979_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1980_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1979_1980_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1979_1980_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1979_1980_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1979 to 1980
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_49.dta", clear

 replace UCID_1934_2003 = 	410665	 if UCID_1934_2003 == 	41097
 replace UCID_1934_2003 = 	5878	 if UCID_1934_2003 == 	42704
 replace UCID_1934_2003 = 	409111	 if UCID_1934_2003 == 	56376
 replace UCID_1934_2003 = 	58727	 if UCID_1934_2003 == 	62570
 replace UCID_1934_2003 = 	53369	 if UCID_1934_2003 == 	76075
 replace UCID_1934_2003 = 	181956	 if UCID_1934_2003 == 	88350
 replace UCID_1934_2003 = 	405870	 if UCID_1934_2003 == 	88535
 replace UCID_1934_2003 = 	112284	 if UCID_1934_2003 == 	131487
 replace UCID_1934_2003 = 	144237	 if UCID_1934_2003 == 	143182
 replace UCID_1934_2003 = 	57985	 if UCID_1934_2003 == 	143735
 replace UCID_1934_2003 = 	409270	 if UCID_1934_2003 == 	145436
 replace UCID_1934_2003 = 	200140	 if UCID_1934_2003 == 	193181
 replace UCID_1934_2003 = 	345795	 if UCID_1934_2003 == 	201402
 replace UCID_1934_2003 = 	227362	 if UCID_1934_2003 == 	217992
 replace UCID_1934_2003 = 	250951	 if UCID_1934_2003 == 	250557
 replace UCID_1934_2003 = 	84277	 if UCID_1934_2003 == 	328533
 replace UCID_1934_2003 = 	337207	 if UCID_1934_2003 == 	331729
 replace UCID_1934_2003 = 	67052	 if UCID_1934_2003 == 	339863
 replace UCID_1934_2003 = 	338803	 if UCID_1934_2003 == 	340329
 replace UCID_1934_2003 = 	333720	 if UCID_1934_2003 == 	346371
 replace UCID_1934_2003 = 	324349	 if UCID_1934_2003 == 	350308
 replace UCID_1934_2003 = 	413626	 if UCID_1934_2003 == 	361259
 replace UCID_1934_2003 = 	405424	 if UCID_1934_2003 == 	370307
 replace UCID_1934_2003 = 	189806	 if UCID_1934_2003 == 	411634
 replace UCID_1934_2003 = 	411731	 if UCID_1934_2003 == 	412808
 replace UCID_1934_2003 = 	412264	 if UCID_1934_2003 == 	412839
 replace UCID_1934_2003 = 	189805	 if UCID_1934_2003 == 	413121
 replace UCID_1934_2003 = 	414011	 if UCID_1934_2003 == 	416257
 replace UCID_1934_2003 = 	146518	 if UCID_1934_2003 == 	430770

 * Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_50.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_49.dta"

*******************************************************************************
* 1980 & 1981

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_50.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1980 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1980
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1980
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1980_remaining_2.dta", replace
	restore

*Second year 1981 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1981
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1981
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1981_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1980 & 1981
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1980_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1981_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1980_1981_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1980_1981_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1980_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1981_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1980_1981_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1980_1981_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1980_1981_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1980 to 1981
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_50.dta", clear

 replace UCID_1934_2003 = 	118569	 if UCID_1934_2003 == 	9659
 replace UCID_1934_2003 = 	162560	 if UCID_1934_2003 == 	18693
 replace UCID_1934_2003 = 	403956	 if UCID_1934_2003 == 	20863
 replace UCID_1934_2003 = 	27813	 if UCID_1934_2003 == 	47742
 replace UCID_1934_2003 = 	1107082	 if UCID_1934_2003 == 	47918
 replace UCID_1934_2003 = 	108094	 if UCID_1934_2003 == 	63924
 replace UCID_1934_2003 = 	50163	 if UCID_1934_2003 == 	71529
 replace UCID_1934_2003 = 	63628	 if UCID_1934_2003 == 	76337
 replace UCID_1934_2003 = 	52513	 if UCID_1934_2003 == 	83086
 replace UCID_1934_2003 = 	276629	 if UCID_1934_2003 == 	87545
 replace UCID_1934_2003 = 	61309	 if UCID_1934_2003 == 	91567
 replace UCID_1934_2003 = 	296010	 if UCID_1934_2003 == 	93705
 replace UCID_1934_2003 = 	107955	 if UCID_1934_2003 == 	104179
 replace UCID_1934_2003 = 	222705	 if UCID_1934_2003 == 	105876
 replace UCID_1934_2003 = 	118044	 if UCID_1934_2003 == 	113284
 replace UCID_1934_2003 = 	119433	 if UCID_1934_2003 == 	117685
 replace UCID_1934_2003 = 	114587	 if UCID_1934_2003 == 	118029
 replace UCID_1934_2003 = 	110760	 if UCID_1934_2003 == 	119351
 replace UCID_1934_2003 = 	132465	 if UCID_1934_2003 == 	133083
 replace UCID_1934_2003 = 	130402	 if UCID_1934_2003 == 	134025
 replace UCID_1934_2003 = 	187021	 if UCID_1934_2003 == 	137088
 replace UCID_1934_2003 = 	9150	 if UCID_1934_2003 == 	137641
 replace UCID_1934_2003 = 	183303	 if UCID_1934_2003 == 	143082
 replace UCID_1934_2003 = 	416528	 if UCID_1934_2003 == 	149414
 replace UCID_1934_2003 = 	155238	 if UCID_1934_2003 == 	154617
 replace UCID_1934_2003 = 	158085	 if UCID_1934_2003 == 	160968
 replace UCID_1934_2003 = 	162403	 if UCID_1934_2003 == 	170604
 replace UCID_1934_2003 = 	277279	 if UCID_1934_2003 == 	171118
 replace UCID_1934_2003 = 	169209	 if UCID_1934_2003 == 	171704
 replace UCID_1934_2003 = 	334743	 if UCID_1934_2003 == 	181929
 replace UCID_1934_2003 = 	404350	 if UCID_1934_2003 == 	186442
 replace UCID_1934_2003 = 	386446	 if UCID_1934_2003 == 	188298
 replace UCID_1934_2003 = 	163866	 if UCID_1934_2003 == 	188479
 replace UCID_1934_2003 = 	391134	 if UCID_1934_2003 == 	191993
 replace UCID_1934_2003 = 	205297	 if UCID_1934_2003 == 	209005
 replace UCID_1934_2003 = 	81188	 if UCID_1934_2003 == 	217838
 replace UCID_1934_2003 = 	232409	 if UCID_1934_2003 == 	231778
 replace UCID_1934_2003 = 	413297	 if UCID_1934_2003 == 	246241
 replace UCID_1934_2003 = 	47772	 if UCID_1934_2003 == 	258813
 replace UCID_1934_2003 = 	253969	 if UCID_1934_2003 == 	263548
 replace UCID_1934_2003 = 	77559	 if UCID_1934_2003 == 	274829
 replace UCID_1934_2003 = 	285266	 if UCID_1934_2003 == 	284468
 replace UCID_1934_2003 = 	162756	 if UCID_1934_2003 == 	304489
 replace UCID_1934_2003 = 	274087	 if UCID_1934_2003 == 	305452
 replace UCID_1934_2003 = 	313596	 if UCID_1934_2003 == 	307975
 replace UCID_1934_2003 = 	366029	 if UCID_1934_2003 == 	328850
 replace UCID_1934_2003 = 	7506	 if UCID_1934_2003 == 	344322
 replace UCID_1934_2003 = 	147586	 if UCID_1934_2003 == 	350871
 replace UCID_1934_2003 = 	357375	 if UCID_1934_2003 == 	356309
 replace UCID_1934_2003 = 	153115	 if UCID_1934_2003 == 	357011
 replace UCID_1934_2003 = 	22548	 if UCID_1934_2003 == 	375498
 replace UCID_1934_2003 = 	383599	 if UCID_1934_2003 == 	394900
 replace UCID_1934_2003 = 	47837	 if UCID_1934_2003 == 	409730
 replace UCID_1934_2003 = 	240444	 if UCID_1934_2003 == 	411255
 replace UCID_1934_2003 = 	410855	 if UCID_1934_2003 == 	413818
 replace UCID_1934_2003 = 	56725	 if UCID_1934_2003 == 	414060
 replace UCID_1934_2003 = 	250006	 if UCID_1934_2003 == 	414135
 replace UCID_1934_2003 = 	47937	 if UCID_1934_2003 == 	416393
 replace UCID_1934_2003 = 	409172	 if UCID_1934_2003 == 	416410

 * Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0
*list if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_51.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_50.dta"

*******************************************************************************
* 1981 & 1982

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_51.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1981 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1981
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1981
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1981_remaining_2.dta", replace
	restore

*Second year 1982 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1982
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1982
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1982_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1981 & 1982
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1981_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1982_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1981_1982_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1981_1982_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1981_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1982_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1981_1982_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1981_1982_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1981_1982_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1981 to 1982
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_51.dta", clear

 replace UCID_1934_2003 = 	317	 if UCID_1934_2003 == 	383
 replace UCID_1934_2003 = 	409726	 if UCID_1934_2003 == 	9088
 replace UCID_1934_2003 = 	116453	 if UCID_1934_2003 == 	11186
 replace UCID_1934_2003 = 	167590	 if UCID_1934_2003 == 	19677
 replace UCID_1934_2003 = 	409890	 if UCID_1934_2003 == 	29457
 replace UCID_1934_2003 = 	43079	 if UCID_1934_2003 == 	46765
 replace UCID_1934_2003 = 	169422	 if UCID_1934_2003 == 	47635
 replace UCID_1934_2003 = 	136382	 if UCID_1934_2003 == 	54024
 replace UCID_1934_2003 = 	246609	 if UCID_1934_2003 == 	64234
 replace UCID_1934_2003 = 	65237	 if UCID_1934_2003 == 	64837
 replace UCID_1934_2003 = 	412452	 if UCID_1934_2003 == 	84683
 replace UCID_1934_2003 = 	377355	 if UCID_1934_2003 == 	87066
 replace UCID_1934_2003 = 	312827	 if UCID_1934_2003 == 	105991
 replace UCID_1934_2003 = 	74537	 if UCID_1934_2003 == 	115018
 replace UCID_1934_2003 = 	400448	 if UCID_1934_2003 == 	121612
 replace UCID_1934_2003 = 	132738	 if UCID_1934_2003 == 	132722
 replace UCID_1934_2003 = 	13017	 if UCID_1934_2003 == 	142461
 replace UCID_1934_2003 = 	209778	 if UCID_1934_2003 == 	143941
 replace UCID_1934_2003 = 	143944	 if UCID_1934_2003 == 	144430
 replace UCID_1934_2003 = 	145261	 if UCID_1934_2003 == 	145617
 replace UCID_1934_2003 = 	14144	 if UCID_1934_2003 == 	146750
 replace UCID_1934_2003 = 	141625	 if UCID_1934_2003 == 	172852
 replace UCID_1934_2003 = 	150336	 if UCID_1934_2003 == 	175230
 replace UCID_1934_2003 = 	284163	 if UCID_1934_2003 == 	185839
 replace UCID_1934_2003 = 	346367	 if UCID_1934_2003 == 	212168
 replace UCID_1934_2003 = 	341984	 if UCID_1934_2003 == 	213613
 replace UCID_1934_2003 = 	287757	 if UCID_1934_2003 == 	240855
 replace UCID_1934_2003 = 	140667	 if UCID_1934_2003 == 	245204
 replace UCID_1934_2003 = 	250275	 if UCID_1934_2003 == 	249783
 replace UCID_1934_2003 = 	302570	 if UCID_1934_2003 == 	251312
 replace UCID_1934_2003 = 	301903	 if UCID_1934_2003 == 	262329
 replace UCID_1934_2003 = 	275214	 if UCID_1934_2003 == 	277915
 replace UCID_1934_2003 = 	222202	 if UCID_1934_2003 == 	280296
 replace UCID_1934_2003 = 	794	 if UCID_1934_2003 == 	344612
 replace UCID_1934_2003 = 	372871	 if UCID_1934_2003 == 	368711
 replace UCID_1934_2003 = 	193159	 if UCID_1934_2003 == 	369054
 replace UCID_1934_2003 = 	127705	 if UCID_1934_2003 == 	398703
 replace UCID_1934_2003 = 	189822	 if UCID_1934_2003 == 	409254
 replace UCID_1934_2003 = 	413562	 if UCID_1934_2003 == 	411822
 replace UCID_1934_2003 = 	409858	 if UCID_1934_2003 == 	413496
 replace UCID_1934_2003 = 	292759	 if UCID_1934_2003 == 	414383
 replace UCID_1934_2003 = 	145303	 if UCID_1934_2003 == 	415448
 replace UCID_1934_2003 = 	1079532	 if UCID_1934_2003 == 	418568

 * Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0
*list if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_52.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_51.dta"

*******************************************************************************
* 1982 & 1983

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_52.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1982 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1982
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1982
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1982_remaining_2.dta", replace
	restore

*Second year 1983 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1983
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1983
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1983_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1982 & 1983
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1982_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1983_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1982_1983_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1982_1983_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1982_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1983_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1982_1983_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1982_1983_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1982_1983_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1982 to 1983
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_52.dta", clear

 replace UCID_1934_2003 = 	9329	 if UCID_1934_2003 == 	1352
 replace UCID_1934_2003 = 	62568	 if UCID_1934_2003 == 	1359
 replace UCID_1934_2003 = 	174743	 if UCID_1934_2003 == 	6676
 replace UCID_1934_2003 = 	138167	 if UCID_1934_2003 == 	13940
 replace UCID_1934_2003 = 	1197	 if UCID_1934_2003 == 	18225
 replace UCID_1934_2003 = 	27541	 if UCID_1934_2003 == 	27447
 replace UCID_1934_2003 = 	410116	 if UCID_1934_2003 == 	29870
 replace UCID_1934_2003 = 	48120	 if UCID_1934_2003 == 	48304
 replace UCID_1934_2003 = 	268337	 if UCID_1934_2003 == 	54752
 replace UCID_1934_2003 = 	3583	 if UCID_1934_2003 == 	57270
 replace UCID_1934_2003 = 	108289	 if UCID_1934_2003 == 	63971
 replace UCID_1934_2003 = 	67490	 if UCID_1934_2003 == 	69419
 replace UCID_1934_2003 = 	60510	 if UCID_1934_2003 == 	80351
 replace UCID_1934_2003 = 	162337	 if UCID_1934_2003 == 	88501
 replace UCID_1934_2003 = 	233315	 if UCID_1934_2003 == 	90830
 replace UCID_1934_2003 = 	52707	 if UCID_1934_2003 == 	97833
 replace UCID_1934_2003 = 	53133	 if UCID_1934_2003 == 	97943
 replace UCID_1934_2003 = 	175353	 if UCID_1934_2003 == 	98079
 replace UCID_1934_2003 = 	117218	 if UCID_1934_2003 == 	111891
 replace UCID_1934_2003 = 	81614	 if UCID_1934_2003 == 	115103
 replace UCID_1934_2003 = 	239290	 if UCID_1934_2003 == 	116390
 replace UCID_1934_2003 = 	119134	 if UCID_1934_2003 == 	119358
 replace UCID_1934_2003 = 	125252	 if UCID_1934_2003 == 	119433
 replace UCID_1934_2003 = 	119280	 if UCID_1934_2003 == 	119451
 replace UCID_1934_2003 = 	126223	 if UCID_1934_2003 == 	126647
 replace UCID_1934_2003 = 	414512	 if UCID_1934_2003 == 	139631
 replace UCID_1934_2003 = 	151429	 if UCID_1934_2003 == 	150787
 replace UCID_1934_2003 = 	58284	 if UCID_1934_2003 == 	151043
 replace UCID_1934_2003 = 	151326	 if UCID_1934_2003 == 	151590
 replace UCID_1934_2003 = 	151477	 if UCID_1934_2003 == 	151734
 replace UCID_1934_2003 = 	204551	 if UCID_1934_2003 == 	152769
 replace UCID_1934_2003 = 	14860	 if UCID_1934_2003 == 	152839
 replace UCID_1934_2003 = 	185441	 if UCID_1934_2003 == 	163459
 replace UCID_1934_2003 = 	154815	 if UCID_1934_2003 == 	173297
 replace UCID_1934_2003 = 	131161	 if UCID_1934_2003 == 	175256
 replace UCID_1934_2003 = 	157011	 if UCID_1934_2003 == 	178776
 replace UCID_1934_2003 = 	409093	 if UCID_1934_2003 == 	184107
 replace UCID_1934_2003 = 	191236	 if UCID_1934_2003 == 	190750
 replace UCID_1934_2003 = 	202194	 if UCID_1934_2003 == 	193748
 replace UCID_1934_2003 = 	393577	 if UCID_1934_2003 == 	214775
 replace UCID_1934_2003 = 	79185	 if UCID_1934_2003 == 	221566
 replace UCID_1934_2003 = 	210574	 if UCID_1934_2003 == 	229446
 replace UCID_1934_2003 = 	236739	 if UCID_1934_2003 == 	235553
 replace UCID_1934_2003 = 	254801	 if UCID_1934_2003 == 	238332
 replace UCID_1934_2003 = 	252319	 if UCID_1934_2003 == 	238908
 replace UCID_1934_2003 = 	234886	 if UCID_1934_2003 == 	241135
 replace UCID_1934_2003 = 	102298	 if UCID_1934_2003 == 	241298
 replace UCID_1934_2003 = 	159474	 if UCID_1934_2003 == 	242452
 replace UCID_1934_2003 = 	362023	 if UCID_1934_2003 == 	248473
 replace UCID_1934_2003 = 	28162	 if UCID_1934_2003 == 	252346
 replace UCID_1934_2003 = 	274946	 if UCID_1934_2003 == 	268321
 replace UCID_1934_2003 = 	32772	 if UCID_1934_2003 == 	271825
 replace UCID_1934_2003 = 	283964	 if UCID_1934_2003 == 	287754
 replace UCID_1934_2003 = 	312301	 if UCID_1934_2003 == 	307303
 replace UCID_1934_2003 = 	40844	 if UCID_1934_2003 == 	323940
 replace UCID_1934_2003 = 	20243	 if UCID_1934_2003 == 	331904
 replace UCID_1934_2003 = 	386088	 if UCID_1934_2003 == 	347807
 replace UCID_1934_2003 = 	366829	 if UCID_1934_2003 == 	367390
 replace UCID_1934_2003 = 	195097	 if UCID_1934_2003 == 	367965
 replace UCID_1934_2003 = 	370128	 if UCID_1934_2003 == 	368993
 replace UCID_1934_2003 = 	359024	 if UCID_1934_2003 == 	374031
 replace UCID_1934_2003 = 	80341	 if UCID_1934_2003 == 	394825
 replace UCID_1934_2003 = 	67247	 if UCID_1934_2003 == 	400703
 replace UCID_1934_2003 = 	404563	 if UCID_1934_2003 == 	404860
 replace UCID_1934_2003 = 	414969	 if UCID_1934_2003 == 	410005
 replace UCID_1934_2003 = 	414358	 if UCID_1934_2003 == 	410112
 replace UCID_1934_2003 = 	63944	 if UCID_1934_2003 == 	410125
 replace UCID_1934_2003 = 	409842	 if UCID_1934_2003 == 	411065
 replace UCID_1934_2003 = 	414386	 if UCID_1934_2003 == 	414393
 replace UCID_1934_2003 = 	428656	 if UCID_1934_2003 == 	428657
 replace UCID_1934_2003 = 	323919	 if UCID_1934_2003 == 	429282

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0
*list if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_53.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_52.dta"

*******************************************************************************
* 1983 & 1984

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_53.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1983 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1983
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1983
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1983_remaining_2.dta", replace
	restore

*Second year 1984 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1984
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1984
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1984_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1983 & 1984
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1983_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1984_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1983_1984_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1983_1984_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1983_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1984_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1983_1984_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1983_1984_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1983_1984_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1983 to 1984
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_53.dta", clear

 replace UCID_1934_2003 = 	53368	 if UCID_1934_2003 == 	595
 replace UCID_1934_2003 = 	54930	 if UCID_1934_2003 == 	5277
 replace UCID_1934_2003 = 	260759	 if UCID_1934_2003 == 	37432
 replace UCID_1934_2003 = 	413556	 if UCID_1934_2003 == 	47594
 replace UCID_1934_2003 = 	158532	 if UCID_1934_2003 == 	47637
 replace UCID_1934_2003 = 	309290	 if UCID_1934_2003 == 	59593
 replace UCID_1934_2003 = 	148216	 if UCID_1934_2003 == 	73628
 replace UCID_1934_2003 = 	332650	 if UCID_1934_2003 == 	104072
 replace UCID_1934_2003 = 	63893	 if UCID_1934_2003 == 	107955
 replace UCID_1934_2003 = 	346849	 if UCID_1934_2003 == 	131162
 replace UCID_1934_2003 = 	132410	 if UCID_1934_2003 == 	146607
 replace UCID_1934_2003 = 	414029	 if UCID_1934_2003 == 	146956
 replace UCID_1934_2003 = 	57836	 if UCID_1934_2003 == 	149828
 replace UCID_1934_2003 = 	83520	 if UCID_1934_2003 == 	152090
 replace UCID_1934_2003 = 	414380	 if UCID_1934_2003 == 	155831
 replace UCID_1934_2003 = 	139946	 if UCID_1934_2003 == 	156845
 replace UCID_1934_2003 = 	384194	 if UCID_1934_2003 == 	172364
 replace UCID_1934_2003 = 	341584	 if UCID_1934_2003 == 	181337
 replace UCID_1934_2003 = 	365269	 if UCID_1934_2003 == 	184159
 replace UCID_1934_2003 = 	166668	 if UCID_1934_2003 == 	187829
 replace UCID_1934_2003 = 	374576	 if UCID_1934_2003 == 	195340
 replace UCID_1934_2003 = 	200369	 if UCID_1934_2003 == 	198893
 replace UCID_1934_2003 = 	119692	 if UCID_1934_2003 == 	204832
 replace UCID_1934_2003 = 	225048	 if UCID_1934_2003 == 	214783
 replace UCID_1934_2003 = 	407534	 if UCID_1934_2003 == 	216163
 replace UCID_1934_2003 = 	152286	 if UCID_1934_2003 == 	239329
 replace UCID_1934_2003 = 	245398	 if UCID_1934_2003 == 	244931
 replace UCID_1934_2003 = 	65015	 if UCID_1934_2003 == 	246566
 replace UCID_1934_2003 = 	280407	 if UCID_1934_2003 == 	266810
 replace UCID_1934_2003 = 	33665	 if UCID_1934_2003 == 	285086
 replace UCID_1934_2003 = 	36813	 if UCID_1934_2003 == 	291487
 replace UCID_1934_2003 = 	77613	 if UCID_1934_2003 == 	309365
 replace UCID_1934_2003 = 	323682	 if UCID_1934_2003 == 	323464
 replace UCID_1934_2003 = 	160005	 if UCID_1934_2003 == 	329527
 replace UCID_1934_2003 = 	373860	 if UCID_1934_2003 == 	332078
 replace UCID_1934_2003 = 	339329	 if UCID_1934_2003 == 	333958
 replace UCID_1934_2003 = 	413119	 if UCID_1934_2003 == 	364792
 replace UCID_1934_2003 = 	42307	 if UCID_1934_2003 == 	365695
 replace UCID_1934_2003 = 	238049	 if UCID_1934_2003 == 	410205
 replace UCID_1934_2003 = 	15258	 if UCID_1934_2003 == 	411698
 replace UCID_1934_2003 = 	154372	 if UCID_1934_2003 == 	412090
 replace UCID_1934_2003 = 	410196	 if UCID_1934_2003 == 	413643
 replace UCID_1934_2003 = 	426290	 if UCID_1934_2003 == 	426289

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0
*list if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_54.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_53.dta"

*******************************************************************************
* 1984 & 1985

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_54.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1984 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1984
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1984
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1984_remaining_2.dta", replace
	restore

*Second year 1985 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1985
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1985
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1985_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1984 & 1985
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1984_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1985_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1984_1985_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1984_1985_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1984_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1985_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1984_1985_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1984_1985_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1984_1985_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1984 to 1985
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_54.dta", clear

 replace UCID_1934_2003 = 	165279	 if UCID_1934_2003 == 	3473
 replace UCID_1934_2003 = 	149826	 if UCID_1934_2003 == 	6161
 replace UCID_1934_2003 = 	117305	 if UCID_1934_2003 == 	10007
 replace UCID_1934_2003 = 	235428	 if UCID_1934_2003 == 	28699
 replace UCID_1934_2003 = 	33860	 if UCID_1934_2003 == 	34496
 replace UCID_1934_2003 = 	406724	 if UCID_1934_2003 == 	45328
 replace UCID_1934_2003 = 	88707	 if UCID_1934_2003 == 	54144
 replace UCID_1934_2003 = 	71957	 if UCID_1934_2003 == 	71658
 replace UCID_1934_2003 = 	268573	 if UCID_1934_2003 == 	71963
 replace UCID_1934_2003 = 	61664	 if UCID_1934_2003 == 	83371
 replace UCID_1934_2003 = 	183508	 if UCID_1934_2003 == 	84024
 replace UCID_1934_2003 = 	68654	 if UCID_1934_2003 == 	85899
 replace UCID_1934_2003 = 	404141	 if UCID_1934_2003 == 	88247
 replace UCID_1934_2003 = 	268458	 if UCID_1934_2003 == 	88749
 replace UCID_1934_2003 = 	261471	 if UCID_1934_2003 == 	94112
 replace UCID_1934_2003 = 	311597	 if UCID_1934_2003 == 	97789
 replace UCID_1934_2003 = 	54029	 if UCID_1934_2003 == 	104643
 replace UCID_1934_2003 = 	117381	 if UCID_1934_2003 == 	111546
 replace UCID_1934_2003 = 	125297	 if UCID_1934_2003 == 	119434
 replace UCID_1934_2003 = 	146623	 if UCID_1934_2003 == 	123642
 replace UCID_1934_2003 = 	116611	 if UCID_1934_2003 == 	130218
 replace UCID_1934_2003 = 	415447	 if UCID_1934_2003 == 	137615
 replace UCID_1934_2003 = 	5289	 if UCID_1934_2003 == 	145695
 replace UCID_1934_2003 = 	89989	 if UCID_1934_2003 == 	154030
 replace UCID_1934_2003 = 	76167	 if UCID_1934_2003 == 	154154
 replace UCID_1934_2003 = 	179228	 if UCID_1934_2003 == 	161538
 replace UCID_1934_2003 = 	431690	 if UCID_1934_2003 == 	175127
 replace UCID_1934_2003 = 	47541	 if UCID_1934_2003 == 	176920
 replace UCID_1934_2003 = 	164754	 if UCID_1934_2003 == 	177383
 replace UCID_1934_2003 = 	415425	 if UCID_1934_2003 == 	178282
 replace UCID_1934_2003 = 	46998	 if UCID_1934_2003 == 	195259
 replace UCID_1934_2003 = 	199563	 if UCID_1934_2003 == 	195624
 replace UCID_1934_2003 = 	164541	 if UCID_1934_2003 == 	199050
 replace UCID_1934_2003 = 	201136	 if UCID_1934_2003 == 	201303
 replace UCID_1934_2003 = 	114028	 if UCID_1934_2003 == 	202553
 replace UCID_1934_2003 = 	205599	 if UCID_1934_2003 == 	205212
 replace UCID_1934_2003 = 	103542	 if UCID_1934_2003 == 	213289
 replace UCID_1934_2003 = 	266245	 if UCID_1934_2003 == 	213510
 replace UCID_1934_2003 = 	236448	 if UCID_1934_2003 == 	237239
 replace UCID_1934_2003 = 	237368	 if UCID_1934_2003 == 	241695
 replace UCID_1934_2003 = 	246098	 if UCID_1934_2003 == 	246471
 replace UCID_1934_2003 = 	347214	 if UCID_1934_2003 == 	249705
 replace UCID_1934_2003 = 	35507	 if UCID_1934_2003 == 	254254
 replace UCID_1934_2003 = 	254303	 if UCID_1934_2003 == 	254443
 replace UCID_1934_2003 = 	54136	 if UCID_1934_2003 == 	266463
 replace UCID_1934_2003 = 	295776	 if UCID_1934_2003 == 	302669
 replace UCID_1934_2003 = 	291435	 if UCID_1934_2003 == 	312872
 replace UCID_1934_2003 = 	164964	 if UCID_1934_2003 == 	371642
 replace UCID_1934_2003 = 	29496	 if UCID_1934_2003 == 	409248
 replace UCID_1934_2003 = 	409386	 if UCID_1934_2003 == 	409345
 replace UCID_1934_2003 = 	15228	 if UCID_1934_2003 == 	410144
 replace UCID_1934_2003 = 	251780	 if UCID_1934_2003 == 	411422
 replace UCID_1934_2003 = 	341130	 if UCID_1934_2003 == 	411574
 replace UCID_1934_2003 = 	415720	 if UCID_1934_2003 == 	414686
 replace UCID_1934_2003 = 	137790	 if UCID_1934_2003 == 	426435
 
* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0
*list if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_55.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_54.dta"

*******************************************************************************
* 1985 & 1986

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_55.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1985 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1985
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1985
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1985_remaining_2.dta", replace
	restore

*Second year 1986 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1986
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1986
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1986_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1985 & 1986
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1985_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1986_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1985_1986_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1985_1986_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1985_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1986_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1985_1986_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1985_1986_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1985_1986_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1985 to 1986
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_55.dta", clear

 replace UCID_1934_2003 = 	62125	 if UCID_1934_2003 == 	9181
 replace UCID_1934_2003 = 	268524	 if UCID_1934_2003 == 	9987
 replace UCID_1934_2003 = 	132335	 if UCID_1934_2003 == 	12471
 replace UCID_1934_2003 = 	154212	 if UCID_1934_2003 == 	14885
 replace UCID_1934_2003 = 	138288	 if UCID_1934_2003 == 	15057
 replace UCID_1934_2003 = 	27200	 if UCID_1934_2003 == 	25133
 replace UCID_1934_2003 = 	428384	 if UCID_1934_2003 == 	30631
 replace UCID_1934_2003 = 	122124	 if UCID_1934_2003 == 	39990
 replace UCID_1934_2003 = 	391966	 if UCID_1934_2003 == 	47097
 replace UCID_1934_2003 = 	32674	 if UCID_1934_2003 == 	47658
 replace UCID_1934_2003 = 	97268	 if UCID_1934_2003 == 	57459
 replace UCID_1934_2003 = 	61991	 if UCID_1934_2003 == 	61256
 replace UCID_1934_2003 = 	60058	 if UCID_1934_2003 == 	61985
 replace UCID_1934_2003 = 	60099	 if UCID_1934_2003 == 	62022
 replace UCID_1934_2003 = 	52136	 if UCID_1934_2003 == 	82248
 replace UCID_1934_2003 = 	161758	 if UCID_1934_2003 == 	100510
 replace UCID_1934_2003 = 	370071	 if UCID_1934_2003 == 	104387
 replace UCID_1934_2003 = 	366474	 if UCID_1934_2003 == 	110251
 replace UCID_1934_2003 = 	119628	 if UCID_1934_2003 == 	113163
 replace UCID_1934_2003 = 	119649	 if UCID_1934_2003 == 	117180
 replace UCID_1934_2003 = 	119639	 if UCID_1934_2003 == 	119556
 replace UCID_1934_2003 = 	119369	 if UCID_1934_2003 == 	119697
 replace UCID_1934_2003 = 	119567	 if UCID_1934_2003 == 	119704
 replace UCID_1934_2003 = 	219120	 if UCID_1934_2003 == 	120398
 replace UCID_1934_2003 = 	108929	 if UCID_1934_2003 == 	128260
 replace UCID_1934_2003 = 	120977	 if UCID_1934_2003 == 	128930
 replace UCID_1934_2003 = 	274076	 if UCID_1934_2003 == 	129431
 replace UCID_1934_2003 = 	132631	 if UCID_1934_2003 == 	134011
 replace UCID_1934_2003 = 	12591	 if UCID_1934_2003 == 	134166
 replace UCID_1934_2003 = 	140887	 if UCID_1934_2003 == 	146423
 replace UCID_1934_2003 = 	148643	 if UCID_1934_2003 == 	148431
 replace UCID_1934_2003 = 	151703	 if UCID_1934_2003 == 	151472
 replace UCID_1934_2003 = 	29471	 if UCID_1934_2003 == 	152029
 replace UCID_1934_2003 = 	157499	 if UCID_1934_2003 == 	155917
 replace UCID_1934_2003 = 	181149	 if UCID_1934_2003 == 	163678
 replace UCID_1934_2003 = 	164948	 if UCID_1934_2003 == 	166568
 replace UCID_1934_2003 = 	180617	 if UCID_1934_2003 == 	168572
 replace UCID_1934_2003 = 	165339	 if UCID_1934_2003 == 	180505
 replace UCID_1934_2003 = 	218066	 if UCID_1934_2003 == 	183275
 replace UCID_1934_2003 = 	190753	 if UCID_1934_2003 == 	183917
 replace UCID_1934_2003 = 	410721	 if UCID_1934_2003 == 	189273
 replace UCID_1934_2003 = 	190844	 if UCID_1934_2003 == 	190401
 replace UCID_1934_2003 = 	191358	 if UCID_1934_2003 == 	199114
 replace UCID_1934_2003 = 	308694	 if UCID_1934_2003 == 	199168
 replace UCID_1934_2003 = 	202672	 if UCID_1934_2003 == 	201939
 replace UCID_1934_2003 = 	203834	 if UCID_1934_2003 == 	206275
 replace UCID_1934_2003 = 	364672	 if UCID_1934_2003 == 	208250
 replace UCID_1934_2003 = 	416516	 if UCID_1934_2003 == 	221883
 replace UCID_1934_2003 = 	224031	 if UCID_1934_2003 == 	221900
 replace UCID_1934_2003 = 	224480	 if UCID_1934_2003 == 	225959
 replace UCID_1934_2003 = 	238913	 if UCID_1934_2003 == 	234927
 replace UCID_1934_2003 = 	244271	 if UCID_1934_2003 == 	238328
 replace UCID_1934_2003 = 	409216	 if UCID_1934_2003 == 	240690
 replace UCID_1934_2003 = 	244897	 if UCID_1934_2003 == 	242381
 replace UCID_1934_2003 = 	243236	 if UCID_1934_2003 == 	242846
 replace UCID_1934_2003 = 	243299	 if UCID_1934_2003 == 	242917
 replace UCID_1934_2003 = 	243397	 if UCID_1934_2003 == 	242961
 replace UCID_1934_2003 = 	243289	 if UCID_1934_2003 == 	242996
 replace UCID_1934_2003 = 	100767	 if UCID_1934_2003 == 	250859
 replace UCID_1934_2003 = 	234024	 if UCID_1934_2003 == 	252305
 replace UCID_1934_2003 = 	253985	 if UCID_1934_2003 == 	253982
 replace UCID_1934_2003 = 	259214	 if UCID_1934_2003 == 	259207
 replace UCID_1934_2003 = 	339432	 if UCID_1934_2003 == 	265002
 replace UCID_1934_2003 = 	286736	 if UCID_1934_2003 == 	265235
 replace UCID_1934_2003 = 	270317	 if UCID_1934_2003 == 	271159
 replace UCID_1934_2003 = 	251311	 if UCID_1934_2003 == 	273078
 replace UCID_1934_2003 = 	251292	 if UCID_1934_2003 == 	273092
 replace UCID_1934_2003 = 	280761	 if UCID_1934_2003 == 	279759
 replace UCID_1934_2003 = 	291467	 if UCID_1934_2003 == 	289739
 replace UCID_1934_2003 = 	34405	 if UCID_1934_2003 == 	290965
 replace UCID_1934_2003 = 	319732	 if UCID_1934_2003 == 	304490
 replace UCID_1934_2003 = 	314372	 if UCID_1934_2003 == 	313796
 replace UCID_1934_2003 = 	304269	 if UCID_1934_2003 == 	314149
 replace UCID_1934_2003 = 	10195	 if UCID_1934_2003 == 	317516
 replace UCID_1934_2003 = 	348291	 if UCID_1934_2003 == 	322799
 replace UCID_1934_2003 = 	338353	 if UCID_1934_2003 == 	323443
 replace UCID_1934_2003 = 	338838	 if UCID_1934_2003 == 	330215
 replace UCID_1934_2003 = 	338648	 if UCID_1934_2003 == 	331719
 replace UCID_1934_2003 = 	324813	 if UCID_1934_2003 == 	331763
 replace UCID_1934_2003 = 	58864	 if UCID_1934_2003 == 	335893
 replace UCID_1934_2003 = 	351826	 if UCID_1934_2003 == 	351789
 replace UCID_1934_2003 = 	416531	 if UCID_1934_2003 == 	356131
 replace UCID_1934_2003 = 	357382	 if UCID_1934_2003 == 	356161
 replace UCID_1934_2003 = 	421290	 if UCID_1934_2003 == 	356203
 replace UCID_1934_2003 = 	357522	 if UCID_1934_2003 == 	356205
 replace UCID_1934_2003 = 	358038	 if UCID_1934_2003 == 	356223
 replace UCID_1934_2003 = 	357876	 if UCID_1934_2003 == 	356241
 replace UCID_1934_2003 = 	357707	 if UCID_1934_2003 == 	356330
 replace UCID_1934_2003 = 	356249	 if UCID_1934_2003 == 	356406
 replace UCID_1934_2003 = 	357855	 if UCID_1934_2003 == 	356425
 replace UCID_1934_2003 = 	356362	 if UCID_1934_2003 == 	358302
 replace UCID_1934_2003 = 	377147	 if UCID_1934_2003 == 	368489
 replace UCID_1934_2003 = 	40920	 if UCID_1934_2003 == 	377596
 replace UCID_1934_2003 = 	370482	 if UCID_1934_2003 == 	377627
 replace UCID_1934_2003 = 	368782	 if UCID_1934_2003 == 	378996
 replace UCID_1934_2003 = 	314705	 if UCID_1934_2003 == 	385906
 replace UCID_1934_2003 = 	405388	 if UCID_1934_2003 == 	391713
 replace UCID_1934_2003 = 	369173	 if UCID_1934_2003 == 	397313
 replace UCID_1934_2003 = 	404228	 if UCID_1934_2003 == 	400706
 replace UCID_1934_2003 = 	118408	 if UCID_1934_2003 == 	407631
 replace UCID_1934_2003 = 	360828	 if UCID_1934_2003 == 	407985
 replace UCID_1934_2003 = 	414481	 if UCID_1934_2003 == 	409157
 replace UCID_1934_2003 = 	246077	 if UCID_1934_2003 == 	409295
 replace UCID_1934_2003 = 	56367	 if UCID_1934_2003 == 	409480
 replace UCID_1934_2003 = 	132694	 if UCID_1934_2003 == 	409532
 replace UCID_1934_2003 = 	227684	 if UCID_1934_2003 == 	409765
 replace UCID_1934_2003 = 	415647	 if UCID_1934_2003 == 	409829
 replace UCID_1934_2003 = 	227396	 if UCID_1934_2003 == 	409847
 replace UCID_1934_2003 = 	431	 if UCID_1934_2003 == 	409921
 replace UCID_1934_2003 = 	409868	 if UCID_1934_2003 == 	409925
 replace UCID_1934_2003 = 	412096	 if UCID_1934_2003 == 	409939
 replace UCID_1934_2003 = 	416193	 if UCID_1934_2003 == 	409994
 replace UCID_1934_2003 = 	410309	 if UCID_1934_2003 == 	410102
 replace UCID_1934_2003 = 	411350	 if UCID_1934_2003 == 	410116
 replace UCID_1934_2003 = 	21480	 if UCID_1934_2003 == 	410118
 replace UCID_1934_2003 = 	411589	 if UCID_1934_2003 == 	410135
 replace UCID_1934_2003 = 	246053	 if UCID_1934_2003 == 	410198
 replace UCID_1934_2003 = 	248419	 if UCID_1934_2003 == 	410242
 replace UCID_1934_2003 = 	409173	 if UCID_1934_2003 == 	410353
 replace UCID_1934_2003 = 	410649	 if UCID_1934_2003 == 	410359
 replace UCID_1934_2003 = 	413618	 if UCID_1934_2003 == 	410363
 replace UCID_1934_2003 = 	412446	 if UCID_1934_2003 == 	410505
 replace UCID_1934_2003 = 	411961	 if UCID_1934_2003 == 	410524
 replace UCID_1934_2003 = 	414327	 if UCID_1934_2003 == 	410899
 replace UCID_1934_2003 = 	204296	 if UCID_1934_2003 == 	410941
 replace UCID_1934_2003 = 	118121	 if UCID_1934_2003 == 	411066
 replace UCID_1934_2003 = 	292116	 if UCID_1934_2003 == 	411067
 replace UCID_1934_2003 = 	147013	 if UCID_1934_2003 == 	411180
 replace UCID_1934_2003 = 	409365	 if UCID_1934_2003 == 	411382
 replace UCID_1934_2003 = 	412602	 if UCID_1934_2003 == 	411385
 replace UCID_1934_2003 = 	273618	 if UCID_1934_2003 == 	411575
 replace UCID_1934_2003 = 	414243	 if UCID_1934_2003 == 	411577
 replace UCID_1934_2003 = 	122660	 if UCID_1934_2003 == 	411678
 replace UCID_1934_2003 = 	412086	 if UCID_1934_2003 == 	411844
 replace UCID_1934_2003 = 	416420	 if UCID_1934_2003 == 	411896
 replace UCID_1934_2003 = 	412138	 if UCID_1934_2003 == 	411936
 replace UCID_1934_2003 = 	247281	 if UCID_1934_2003 == 	411973
 replace UCID_1934_2003 = 	412494	 if UCID_1934_2003 == 	412478
 replace UCID_1934_2003 = 	413775	 if UCID_1934_2003 == 	412727
 replace UCID_1934_2003 = 	414446	 if UCID_1934_2003 == 	413257
 replace UCID_1934_2003 = 	415281	 if UCID_1934_2003 == 	413415
 replace UCID_1934_2003 = 	411458	 if UCID_1934_2003 == 	413828
 replace UCID_1934_2003 = 	341854	 if UCID_1934_2003 == 	414021
 replace UCID_1934_2003 = 	415658	 if UCID_1934_2003 == 	414123
 replace UCID_1934_2003 = 	411945	 if UCID_1934_2003 == 	414332
 replace UCID_1934_2003 = 	2472	 if UCID_1934_2003 == 	414610
 replace UCID_1934_2003 = 	414244	 if UCID_1934_2003 == 	415117
 replace UCID_1934_2003 = 	415048	 if UCID_1934_2003 == 	415125
 replace UCID_1934_2003 = 	414866	 if UCID_1934_2003 == 	415172
 replace UCID_1934_2003 = 	70497	 if UCID_1934_2003 == 	415307
 replace UCID_1934_2003 = 	411353	 if UCID_1934_2003 == 	415419
 replace UCID_1934_2003 = 	410210	 if UCID_1934_2003 == 	415591
 replace UCID_1934_2003 = 	132183	 if UCID_1934_2003 == 	415639
 replace UCID_1934_2003 = 	146954	 if UCID_1934_2003 == 	415694
 replace UCID_1934_2003 = 	409952	 if UCID_1934_2003 == 	415723
 replace UCID_1934_2003 = 	258419	 if UCID_1934_2003 == 	415742
 replace UCID_1934_2003 = 	227880	 if UCID_1934_2003 == 	415787
 replace UCID_1934_2003 = 	413582	 if UCID_1934_2003 == 	415796
 replace UCID_1934_2003 = 	426183	 if UCID_1934_2003 == 	415808
 replace UCID_1934_2003 = 	415297	 if UCID_1934_2003 == 	415823
 replace UCID_1934_2003 = 	422732	 if UCID_1934_2003 == 	418766
 replace UCID_1934_2003 = 	426824	 if UCID_1934_2003 == 	426823
 replace UCID_1934_2003 = 	189172	 if UCID_1934_2003 == 	430968
 replace UCID_1934_2003 = 	431921	 if UCID_1934_2003 == 	1025096

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0
*list if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_56.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_55.dta"

*******************************************************************************
* 1986 & 1987

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_56.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1986 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1986
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1986
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1986_remaining_2.dta", replace
	restore

*Second year 1987 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1987
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1987
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1987_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1986 & 1987
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1986_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1987_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1986_1987_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1986_1987_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1986_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1987_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1986_1987_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1986_1987_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1986_1987_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1986 to 1987
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_56.dta", clear

 replace UCID_1934_2003 = 	49618	 if UCID_1934_2003 == 	4812
 replace UCID_1934_2003 = 	1051	 if UCID_1934_2003 == 	5515
 replace UCID_1934_2003 = 	156085	 if UCID_1934_2003 == 	13664
 replace UCID_1934_2003 = 	17798	 if UCID_1934_2003 == 	24403
 replace UCID_1934_2003 = 	27526	 if UCID_1934_2003 == 	24617
 replace UCID_1934_2003 = 	905391	 if UCID_1934_2003 == 	25352
 replace UCID_1934_2003 = 	28087	 if UCID_1934_2003 == 	28110
 replace UCID_1934_2003 = 	15143	 if UCID_1934_2003 == 	28601
 replace UCID_1934_2003 = 	425493	 if UCID_1934_2003 == 	32666
 replace UCID_1934_2003 = 	362471	 if UCID_1934_2003 == 	42584
 replace UCID_1934_2003 = 	32428	 if UCID_1934_2003 == 	47486
 replace UCID_1934_2003 = 	34072	 if UCID_1934_2003 == 	47530
 replace UCID_1934_2003 = 	29984	 if UCID_1934_2003 == 	47550
 replace UCID_1934_2003 = 	433022	 if UCID_1934_2003 == 	47649
 replace UCID_1934_2003 = 	245953	 if UCID_1934_2003 == 	47667
 replace UCID_1934_2003 = 	708	 if UCID_1934_2003 == 	47841
 replace UCID_1934_2003 = 	15181	 if UCID_1934_2003 == 	47947
 replace UCID_1934_2003 = 	29762	 if UCID_1934_2003 == 	48103
 replace UCID_1934_2003 = 	32835	 if UCID_1934_2003 == 	48108
 replace UCID_1934_2003 = 	250204	 if UCID_1934_2003 == 	58330
 replace UCID_1934_2003 = 	64443	 if UCID_1934_2003 == 	63021
 replace UCID_1934_2003 = 	143003	 if UCID_1934_2003 == 	63786
 replace UCID_1934_2003 = 	270500	 if UCID_1934_2003 == 	70660
 replace UCID_1934_2003 = 	279771	 if UCID_1934_2003 == 	73574
 replace UCID_1934_2003 = 	52760	 if UCID_1934_2003 == 	74130
 replace UCID_1934_2003 = 	91783	 if UCID_1934_2003 == 	78429
 replace UCID_1934_2003 = 	62504	 if UCID_1934_2003 == 	80736
 replace UCID_1934_2003 = 	51099	 if UCID_1934_2003 == 	90178
 replace UCID_1934_2003 = 	178326	 if UCID_1934_2003 == 	143376
 replace UCID_1934_2003 = 	145254	 if UCID_1934_2003 == 	145953
 replace UCID_1934_2003 = 	150929	 if UCID_1934_2003 == 	148408
 replace UCID_1934_2003 = 	271824	 if UCID_1934_2003 == 	153862
 replace UCID_1934_2003 = 	156922	 if UCID_1934_2003 == 	156489
 replace UCID_1934_2003 = 	156473	 if UCID_1934_2003 == 	157049
 replace UCID_1934_2003 = 	191568	 if UCID_1934_2003 == 	191695
 replace UCID_1934_2003 = 	189171	 if UCID_1934_2003 == 	195255
 replace UCID_1934_2003 = 	199845	 if UCID_1934_2003 == 	198048
 replace UCID_1934_2003 = 	202369	 if UCID_1934_2003 == 	202025
 replace UCID_1934_2003 = 	351312	 if UCID_1934_2003 == 	204071
 replace UCID_1934_2003 = 	346728	 if UCID_1934_2003 == 	206149
 replace UCID_1934_2003 = 	269759	 if UCID_1934_2003 == 	213678
 replace UCID_1934_2003 = 	136600	 if UCID_1934_2003 == 	222194
 replace UCID_1934_2003 = 	27017	 if UCID_1934_2003 == 	223816
 replace UCID_1934_2003 = 	235180	 if UCID_1934_2003 == 	234263
 replace UCID_1934_2003 = 	235268	 if UCID_1934_2003 == 	234491
 replace UCID_1934_2003 = 	242779	 if UCID_1934_2003 == 	244528
 replace UCID_1934_2003 = 	243334	 if UCID_1934_2003 == 	245943
 replace UCID_1934_2003 = 	158567	 if UCID_1934_2003 == 	246541
 replace UCID_1934_2003 = 	314388	 if UCID_1934_2003 == 	247247
 replace UCID_1934_2003 = 	262523	 if UCID_1934_2003 == 	255259
 replace UCID_1934_2003 = 	407032	 if UCID_1934_2003 == 	264414
 replace UCID_1934_2003 = 	296503	 if UCID_1934_2003 == 	265697
 replace UCID_1934_2003 = 	24594	 if UCID_1934_2003 == 	266885
 replace UCID_1934_2003 = 	155893	 if UCID_1934_2003 == 	269696
 replace UCID_1934_2003 = 	271880	 if UCID_1934_2003 == 	269735
 replace UCID_1934_2003 = 	64697	 if UCID_1934_2003 == 	270908
 replace UCID_1934_2003 = 	281248	 if UCID_1934_2003 == 	274308
 replace UCID_1934_2003 = 	62862	 if UCID_1934_2003 == 	274381
 replace UCID_1934_2003 = 	274258	 if UCID_1934_2003 == 	274885
 replace UCID_1934_2003 = 	32928	 if UCID_1934_2003 == 	275760
 replace UCID_1934_2003 = 	270924	 if UCID_1934_2003 == 	280820
 replace UCID_1934_2003 = 	313534	 if UCID_1934_2003 == 	283060
 replace UCID_1934_2003 = 	314370	 if UCID_1934_2003 == 	291457
 replace UCID_1934_2003 = 	290477	 if UCID_1934_2003 == 	293233
 replace UCID_1934_2003 = 	257275	 if UCID_1934_2003 == 	309874
 replace UCID_1934_2003 = 	320043	 if UCID_1934_2003 == 	319588
 replace UCID_1934_2003 = 	291470	 if UCID_1934_2003 == 	321229
 replace UCID_1934_2003 = 	345639	 if UCID_1934_2003 == 	323273
 replace UCID_1934_2003 = 	132974	 if UCID_1934_2003 == 	333724
 replace UCID_1934_2003 = 	156573	 if UCID_1934_2003 == 	335762
 replace UCID_1934_2003 = 	333483	 if UCID_1934_2003 == 	341029
 replace UCID_1934_2003 = 	343815	 if UCID_1934_2003 == 	342895
 replace UCID_1934_2003 = 	364423	 if UCID_1934_2003 == 	362187
 replace UCID_1934_2003 = 	287761	 if UCID_1934_2003 == 	409849
 replace UCID_1934_2003 = 	156144	 if UCID_1934_2003 == 	410082
 replace UCID_1934_2003 = 	248391	 if UCID_1934_2003 == 	410339
 replace UCID_1934_2003 = 	29238	 if UCID_1934_2003 == 	412754
 replace UCID_1934_2003 = 	155602	 if UCID_1934_2003 == 	412793
 replace UCID_1934_2003 = 	135759	 if UCID_1934_2003 == 	414101
 replace UCID_1934_2003 = 	249039	 if UCID_1934_2003 == 	414907
 replace UCID_1934_2003 = 	56185	 if UCID_1934_2003 == 	415283
 replace UCID_1934_2003 = 	248959	 if UCID_1934_2003 == 	415334
 replace UCID_1934_2003 = 	355705	 if UCID_1934_2003 == 	415739
 replace UCID_1934_2003 = 	15256	 if UCID_1934_2003 == 	426751
 replace UCID_1934_2003 = 	25849	 if UCID_1934_2003 == 	429265

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
*keep if dup4 > 0
*list if dup4 > 0

* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_57.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_56.dta"

*******************************************************************************
* 1987 & 1988

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_57.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1987 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1987
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1987
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1987_remaining_2.dta", replace
	restore

*Second year 1988 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1988
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1988
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1988_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1987 & 1988
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1987_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1988_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1987_1988_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1987_1988_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1987_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1988_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1987_1988_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1987_1988_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1987_1988_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1987 to 1988
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_57.dta", clear

 replace UCID_1934_2003 = 	89	 if UCID_1934_2003 == 	1482
 replace UCID_1934_2003 = 	69987	 if UCID_1934_2003 == 	5478
 replace UCID_1934_2003 = 	1310	 if UCID_1934_2003 == 	7090
 replace UCID_1934_2003 = 	117816	 if UCID_1934_2003 == 	9780
 replace UCID_1934_2003 = 	11385	 if UCID_1934_2003 == 	9913
 replace UCID_1934_2003 = 	229116	 if UCID_1934_2003 == 	26790
 replace UCID_1934_2003 = 	25624	 if UCID_1934_2003 == 	32879
 replace UCID_1934_2003 = 	314730	 if UCID_1934_2003 == 	38360
 replace UCID_1934_2003 = 	121689	 if UCID_1934_2003 == 	39361
 replace UCID_1934_2003 = 	404935	 if UCID_1934_2003 == 	44923
 replace UCID_1934_2003 = 	431835	 if UCID_1934_2003 == 	47508
 replace UCID_1934_2003 = 	83030	 if UCID_1934_2003 == 	51723
 replace UCID_1934_2003 = 	69224	 if UCID_1934_2003 == 	56205
 replace UCID_1934_2003 = 	62152	 if UCID_1934_2003 == 	62828
 replace UCID_1934_2003 = 	55478	 if UCID_1934_2003 == 	63593
 replace UCID_1934_2003 = 	287586	 if UCID_1934_2003 == 	66820
 replace UCID_1934_2003 = 	60623	 if UCID_1934_2003 == 	74979
 replace UCID_1934_2003 = 	97312	 if UCID_1934_2003 == 	78590
 replace UCID_1934_2003 = 	57675	 if UCID_1934_2003 == 	86139
 replace UCID_1934_2003 = 	396816	 if UCID_1934_2003 == 	86243
 replace UCID_1934_2003 = 	48999	 if UCID_1934_2003 == 	101609
 replace UCID_1934_2003 = 	271211	 if UCID_1934_2003 == 	108317
 replace UCID_1934_2003 = 	324382	 if UCID_1934_2003 == 	111090
 replace UCID_1934_2003 = 	129401	 if UCID_1934_2003 == 	120654
 replace UCID_1934_2003 = 	138402	 if UCID_1934_2003 == 	138645
 replace UCID_1934_2003 = 	157520	 if UCID_1934_2003 == 	139910
 replace UCID_1934_2003 = 	200279	 if UCID_1934_2003 == 	141984
 replace UCID_1934_2003 = 	146295	 if UCID_1934_2003 == 	155969
 replace UCID_1934_2003 = 	391949	 if UCID_1934_2003 == 	161760
 replace UCID_1934_2003 = 	252368	 if UCID_1934_2003 == 	173560
 replace UCID_1934_2003 = 	97389	 if UCID_1934_2003 == 	179396
 replace UCID_1934_2003 = 	200374	 if UCID_1934_2003 == 	191827
 replace UCID_1934_2003 = 	60130	 if UCID_1934_2003 == 	203088
 replace UCID_1934_2003 = 	195345	 if UCID_1934_2003 == 	214664
 replace UCID_1934_2003 = 	212806	 if UCID_1934_2003 == 	224401
 replace UCID_1934_2003 = 	215890	 if UCID_1934_2003 == 	226243
 replace UCID_1934_2003 = 	251464	 if UCID_1934_2003 == 	234008
 replace UCID_1934_2003 = 	28844	 if UCID_1934_2003 == 	240678
 replace UCID_1934_2003 = 	244256	 if UCID_1934_2003 == 	241836
 replace UCID_1934_2003 = 	207054	 if UCID_1934_2003 == 	260637
 replace UCID_1934_2003 = 	44017	 if UCID_1934_2003 == 	262672
 replace UCID_1934_2003 = 	267778	 if UCID_1934_2003 == 	266161
 replace UCID_1934_2003 = 	267775	 if UCID_1934_2003 == 	266692
 replace UCID_1934_2003 = 	268624	 if UCID_1934_2003 == 	269992
 replace UCID_1934_2003 = 	285153	 if UCID_1934_2003 == 	273893
 replace UCID_1934_2003 = 	276233	 if UCID_1934_2003 == 	275919
 replace UCID_1934_2003 = 	155834	 if UCID_1934_2003 == 	280331
 replace UCID_1934_2003 = 	17503	 if UCID_1934_2003 == 	281939
 replace UCID_1934_2003 = 	188028	 if UCID_1934_2003 == 	300315
 replace UCID_1934_2003 = 	311340	 if UCID_1934_2003 == 	311252
 replace UCID_1934_2003 = 	323701	 if UCID_1934_2003 == 	323688
 replace UCID_1934_2003 = 	338657	 if UCID_1934_2003 == 	325539
 replace UCID_1934_2003 = 	17923	 if UCID_1934_2003 == 	329261
 replace UCID_1934_2003 = 	336826	 if UCID_1934_2003 == 	329402
 replace UCID_1934_2003 = 	341518	 if UCID_1934_2003 == 	332281
 replace UCID_1934_2003 = 	342490	 if UCID_1934_2003 == 	344138
 replace UCID_1934_2003 = 	149304	 if UCID_1934_2003 == 	354743
 replace UCID_1934_2003 = 	44047	 if UCID_1934_2003 == 	370466
 replace UCID_1934_2003 = 	81091	 if UCID_1934_2003 == 	372503
 replace UCID_1934_2003 = 	344797	 if UCID_1934_2003 == 	382742
 replace UCID_1934_2003 = 	370474	 if UCID_1934_2003 == 	387735
 replace UCID_1934_2003 = 	369708	 if UCID_1934_2003 == 	393756
 replace UCID_1934_2003 = 	93192	 if UCID_1934_2003 == 	395655
 replace UCID_1934_2003 = 	137465	 if UCID_1934_2003 == 	411324
 replace UCID_1934_2003 = 	29935	 if UCID_1934_2003 == 	414994
 replace UCID_1934_2003 = 	23362	 if UCID_1934_2003 == 	426979

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_58.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_57.dta"

*******************************************************************************
* 1988 & 1989

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_58.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1988 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1988
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1988
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1988_remaining_2.dta", replace
	restore

*Second year 1989 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1989
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1989
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1989_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1988 & 1989
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1988_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1989_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1988_1989_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1988_1989_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1988_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1989_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1988_1989_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1988_1989_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1988_1989_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1988 to 1989
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_58.dta", clear

 replace UCID_1934_2003 = 	46959	 if UCID_1934_2003 == 	4834
 replace UCID_1934_2003 = 	129468	 if UCID_1934_2003 == 	11499
 replace UCID_1934_2003 = 	40535	 if UCID_1934_2003 == 	13242
 replace UCID_1934_2003 = 	15420	 if UCID_1934_2003 == 	18226
 replace UCID_1934_2003 = 	225542	 if UCID_1934_2003 == 	27046
 replace UCID_1934_2003 = 	232083	 if UCID_1934_2003 == 	49074
 replace UCID_1934_2003 = 	55263	 if UCID_1934_2003 == 	53613
 replace UCID_1934_2003 = 	73027	 if UCID_1934_2003 == 	72513
 replace UCID_1934_2003 = 	72359	 if UCID_1934_2003 == 	79274
 replace UCID_1934_2003 = 	244537	 if UCID_1934_2003 == 	85724
 replace UCID_1934_2003 = 	107463	 if UCID_1934_2003 == 	88828
 replace UCID_1934_2003 = 	226076	 if UCID_1934_2003 == 	93863
 replace UCID_1934_2003 = 	64415	 if UCID_1934_2003 == 	113776
 replace UCID_1934_2003 = 	130097	 if UCID_1934_2003 == 	115705
 replace UCID_1934_2003 = 	27179	 if UCID_1934_2003 == 	122240
 replace UCID_1934_2003 = 	126632	 if UCID_1934_2003 == 	127035
 replace UCID_1934_2003 = 	9401	 if UCID_1934_2003 == 	130104
 replace UCID_1934_2003 = 	142927	 if UCID_1934_2003 == 	144221
 replace UCID_1934_2003 = 	147745	 if UCID_1934_2003 == 	147363
 replace UCID_1934_2003 = 	36000	 if UCID_1934_2003 == 	147658
 replace UCID_1934_2003 = 	86149	 if UCID_1934_2003 == 	159040
 replace UCID_1934_2003 = 	190946	 if UCID_1934_2003 == 	172796
 replace UCID_1934_2003 = 	169002	 if UCID_1934_2003 == 	185585
 replace UCID_1934_2003 = 	47566	 if UCID_1934_2003 == 	195689
 replace UCID_1934_2003 = 	410268	 if UCID_1934_2003 == 	197620
 replace UCID_1934_2003 = 	413726	 if UCID_1934_2003 == 	198229
 replace UCID_1934_2003 = 	203148	 if UCID_1934_2003 == 	198470
 replace UCID_1934_2003 = 	414944	 if UCID_1934_2003 == 	198647
 replace UCID_1934_2003 = 	27427	 if UCID_1934_2003 == 	220921
 replace UCID_1934_2003 = 	252301	 if UCID_1934_2003 == 	237131
 replace UCID_1934_2003 = 	235448	 if UCID_1934_2003 == 	238713
 replace UCID_1934_2003 = 	242564	 if UCID_1934_2003 == 	242050
 replace UCID_1934_2003 = 	245545	 if UCID_1934_2003 == 	243180
 replace UCID_1934_2003 = 	158192	 if UCID_1934_2003 == 	248084
 replace UCID_1934_2003 = 	81337	 if UCID_1934_2003 == 	254161
 replace UCID_1934_2003 = 	270667	 if UCID_1934_2003 == 	268732
 replace UCID_1934_2003 = 	268911	 if UCID_1934_2003 == 	269204
 replace UCID_1934_2003 = 	269899	 if UCID_1934_2003 == 	270691
 replace UCID_1934_2003 = 	226871	 if UCID_1934_2003 == 	275899
 replace UCID_1934_2003 = 	208426	 if UCID_1934_2003 == 	281066
 replace UCID_1934_2003 = 	74098	 if UCID_1934_2003 == 	293699
 replace UCID_1934_2003 = 	314217	 if UCID_1934_2003 == 	296302
 replace UCID_1934_2003 = 	156826	 if UCID_1934_2003 == 	309377
 replace UCID_1934_2003 = 	156408	 if UCID_1934_2003 == 	311942
 replace UCID_1934_2003 = 	820399	 if UCID_1934_2003 == 	313061
 replace UCID_1934_2003 = 	320265	 if UCID_1934_2003 == 	318942
 replace UCID_1934_2003 = 	282815	 if UCID_1934_2003 == 	334760
 replace UCID_1934_2003 = 	347944	 if UCID_1934_2003 == 	349051
 replace UCID_1934_2003 = 	331818	 if UCID_1934_2003 == 	349873
 replace UCID_1934_2003 = 	361721	 if UCID_1934_2003 == 	361845
 replace UCID_1934_2003 = 	60059	 if UCID_1934_2003 == 	363759
 replace UCID_1934_2003 = 	368354	 if UCID_1934_2003 == 	370003
 replace UCID_1934_2003 = 	405138	 if UCID_1934_2003 == 	375235
 replace UCID_1934_2003 = 	321497	 if UCID_1934_2003 == 	392424
 replace UCID_1934_2003 = 	356170	 if UCID_1934_2003 == 	393207
 replace UCID_1934_2003 = 	56979	 if UCID_1934_2003 == 	415378
 replace UCID_1934_2003 = 	12246	 if UCID_1934_2003 == 	426846
 replace UCID_1934_2003 = 	292598	 if UCID_1934_2003 == 	428258
 replace UCID_1934_2003 = 	873835	 if UCID_1934_2003 == 	432601

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_59.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_58.dta"

*******************************************************************************
* 1989 & 1990

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_59.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1989 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1989
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1989
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1989_remaining_2.dta", replace
	restore

*Second year 1990 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1990
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1990
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1990_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1989 & 1990
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1989_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1990_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1989_1990_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1989_1990_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1989_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1990_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1989_1990_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1989_1990_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1989_1990_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1989 to 1990
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_59.dta", clear

replace gdenr = 5590 if UCID_1934_2003 == 353265 & year == 1990
replace gdenr_2012 = 5590 if UCID_1934_2003 == 353265 & year == 1990
replace gdenr_2018 = 5590 if UCID_1934_2003 == 353265 & year == 1990
replace gdename = "pully" if UCID_1934_2003 == 353265 & year == 1990

replace gdenr = 5590 if UCID_1934_2003 == 353277 
replace gdenr_2012 = 5590 if UCID_1934_2003 == 353277
replace gdenr_2018 = 5590 if UCID_1934_2003 == 353277

replace gdenr = 5590 if UCID_1934_2003 == 353302
replace gdenr_2012 = 5590 if UCID_1934_2003 == 353302
replace gdenr_2018 = 5590 if UCID_1934_2003 == 353302

replace gdenr = 5590 if UCID_1934_2003 == 353296
replace gdenr_2012 = 5590 if UCID_1934_2003 == 353296
replace gdenr_2018 = 5590 if UCID_1934_2003 == 353296

replace gdenr = 6133 if UCID_1934_2003 == 353282 
replace gdenr_2012 = 6133 if UCID_1934_2003 == 353282
replace gdenr_2018 = 6133 if UCID_1934_2003 == 353282

replace gdenr = 5590 if UCID_1934_2003 == 353279 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 353279 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 353279 & year == 1991
replace gdename = "pully" if UCID_1934_2003 == 353279 & year == 1991

 replace UCID_1934_2003 = 	49265	 if UCID_1934_2003 == 	41
 replace UCID_1934_2003 = 	17079	 if UCID_1934_2003 == 	2906
 replace UCID_1934_2003 = 	201308	 if UCID_1934_2003 == 	10517
 replace UCID_1934_2003 = 	324302	 if UCID_1934_2003 == 	39072
 replace UCID_1934_2003 = 	175135	 if UCID_1934_2003 == 	39572
 replace UCID_1934_2003 = 	34084	 if UCID_1934_2003 == 	47571
 replace UCID_1934_2003 = 	205131	 if UCID_1934_2003 == 	47735
 replace UCID_1934_2003 = 	340132	 if UCID_1934_2003 == 	47865
 replace UCID_1934_2003 = 	62832	 if UCID_1934_2003 == 	60501
 replace UCID_1934_2003 = 	369637	 if UCID_1934_2003 == 	78685
 replace UCID_1934_2003 = 	52968	 if UCID_1934_2003 == 	85375
 replace UCID_1934_2003 = 	69593	 if UCID_1934_2003 == 	85965
 replace UCID_1934_2003 = 	71840	 if UCID_1934_2003 == 	87161
 replace UCID_1934_2003 = 	181185	 if UCID_1934_2003 == 	87944
 replace UCID_1934_2003 = 	59870	 if UCID_1934_2003 == 	94028
 replace UCID_1934_2003 = 	54424	 if UCID_1934_2003 == 	94228
 replace UCID_1934_2003 = 	71745	 if UCID_1934_2003 == 	95067
 replace UCID_1934_2003 = 	223492	 if UCID_1934_2003 == 	107381
 replace UCID_1934_2003 = 	12745	 if UCID_1934_2003 == 	135193
 replace UCID_1934_2003 = 	144239	 if UCID_1934_2003 == 	137196
 replace UCID_1934_2003 = 	138395	 if UCID_1934_2003 == 	140274
 replace UCID_1934_2003 = 	157290	 if UCID_1934_2003 == 	140782
 replace UCID_1934_2003 = 	137062	 if UCID_1934_2003 == 	144346
 replace UCID_1934_2003 = 	299741	 if UCID_1934_2003 == 	149736
 replace UCID_1934_2003 = 	157669	 if UCID_1934_2003 == 	160682
 replace UCID_1934_2003 = 	97436	 if UCID_1934_2003 == 	172106
 replace UCID_1934_2003 = 	162520	 if UCID_1934_2003 == 	185880
 replace UCID_1934_2003 = 	152798	 if UCID_1934_2003 == 	186289
 replace UCID_1934_2003 = 	60432	 if UCID_1934_2003 == 	193525
 replace UCID_1934_2003 = 	191528	 if UCID_1934_2003 == 	201386
 replace UCID_1934_2003 = 	227628	 if UCID_1934_2003 == 	215784
 replace UCID_1934_2003 = 	394514	 if UCID_1934_2003 == 	218803
 replace UCID_1934_2003 = 	227916	 if UCID_1934_2003 == 	220892
 replace UCID_1934_2003 = 	224871	 if UCID_1934_2003 == 	229058
 replace UCID_1934_2003 = 	233376	 if UCID_1934_2003 == 	231637
 replace UCID_1934_2003 = 	367440	 if UCID_1934_2003 == 	238644
 replace UCID_1934_2003 = 	246531	 if UCID_1934_2003 == 	245974
 replace UCID_1934_2003 = 	282470	 if UCID_1934_2003 == 	250562
 replace UCID_1934_2003 = 	241689	 if UCID_1934_2003 == 	250636
 replace UCID_1934_2003 = 	288992	 if UCID_1934_2003 == 	251972
 replace UCID_1934_2003 = 	279901	 if UCID_1934_2003 == 	268769
 replace UCID_1934_2003 = 	267454	 if UCID_1934_2003 == 	269108
 replace UCID_1934_2003 = 	273377	 if UCID_1934_2003 == 	277279
 replace UCID_1934_2003 = 	17926	 if UCID_1934_2003 == 	281844
 replace UCID_1934_2003 = 	287005	 if UCID_1934_2003 == 	285844
 replace UCID_1934_2003 = 	291450	 if UCID_1934_2003 == 	289479
 replace UCID_1934_2003 = 	264456	 if UCID_1934_2003 == 	292149
 replace UCID_1934_2003 = 	264397	 if UCID_1934_2003 == 	292150
 replace UCID_1934_2003 = 	291445	 if UCID_1934_2003 == 	294028
 replace UCID_1934_2003 = 	291257	 if UCID_1934_2003 == 	294765
 replace UCID_1934_2003 = 	5761	 if UCID_1934_2003 == 	295005
 replace UCID_1934_2003 = 	320773	 if UCID_1934_2003 == 	299713
 replace UCID_1934_2003 = 	304145	 if UCID_1934_2003 == 	317620
 replace UCID_1934_2003 = 	320221	 if UCID_1934_2003 == 	319849
 replace UCID_1934_2003 = 	320241	 if UCID_1934_2003 == 	320228
 replace UCID_1934_2003 = 	319516	 if UCID_1934_2003 == 	320234
 replace UCID_1934_2003 = 	339761	 if UCID_1934_2003 == 	331759
 replace UCID_1934_2003 = 	338824	 if UCID_1934_2003 == 	337384
 replace UCID_1934_2003 = 	353265	 if UCID_1934_2003 == 	338113
 replace UCID_1934_2003 = 	321553	 if UCID_1934_2003 == 	347027
 replace UCID_1934_2003 = 	348160	 if UCID_1934_2003 == 	348942
 replace UCID_1934_2003 = 	205284	 if UCID_1934_2003 == 	361533
 replace UCID_1934_2003 = 	156557	 if UCID_1934_2003 == 	361673
 replace UCID_1934_2003 = 	70888	 if UCID_1934_2003 == 	364912
 replace UCID_1934_2003 = 	180009	 if UCID_1934_2003 == 	373352
 replace UCID_1934_2003 = 	370131	 if UCID_1934_2003 == 	385755
 replace UCID_1934_2003 = 	403345	 if UCID_1934_2003 == 	386263
 replace UCID_1934_2003 = 	101493	 if UCID_1934_2003 == 	386660
 replace UCID_1934_2003 = 	164592	 if UCID_1934_2003 == 	389598
 replace UCID_1934_2003 = 	347409	 if UCID_1934_2003 == 	390193
 replace UCID_1934_2003 = 	429761	 if UCID_1934_2003 == 	429760
 replace UCID_1934_2003 = 	402665	 if UCID_1934_2003 == 	393996
 replace UCID_1934_2003 = 	370990	 if UCID_1934_2003 == 	403253
 replace UCID_1934_2003 = 	245018	 if UCID_1934_2003 == 	415648
 replace UCID_1934_2003 = 	60758	 if UCID_1934_2003 == 	424191
 replace UCID_1934_2003 = 	52494	 if UCID_1934_2003 == 	905506
 replace UCID_1934_2003 = 	393124	 if UCID_1934_2003 == 	1381923

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_60.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_59.dta"

*******************************************************************************
* 1990 & 1991

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_60.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1990 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1990
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1990
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1990_remaining_2.dta", replace
	restore

*Second year 1991 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1991
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1991
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1991_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1990 & 1991
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1990_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1991_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1990_1991_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1990_1991_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1990_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1991_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1990_1991_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1990_1991_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1990_1991_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1990 to 1991
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_60.dta", clear

 replace UCID_1934_2003 = 	897	 if UCID_1934_2003 == 	848
 replace UCID_1934_2003 = 	1216	 if UCID_1934_2003 == 	2585
 replace UCID_1934_2003 = 	229491	 if UCID_1934_2003 == 	4500
 replace UCID_1934_2003 = 	19291	 if UCID_1934_2003 == 	5420
 replace UCID_1934_2003 = 	318	 if UCID_1934_2003 == 	5891
 replace UCID_1934_2003 = 	55489	 if UCID_1934_2003 == 	7588
 replace UCID_1934_2003 = 	343521	 if UCID_1934_2003 == 	14695
 replace UCID_1934_2003 = 	39176	 if UCID_1934_2003 == 	15291
 replace UCID_1934_2003 = 	16639	 if UCID_1934_2003 == 	17751
 replace UCID_1934_2003 = 	162599	 if UCID_1934_2003 == 	18160
 replace UCID_1934_2003 = 	136595	 if UCID_1934_2003 == 	22322
 replace UCID_1934_2003 = 	184097	 if UCID_1934_2003 == 	23342
 replace UCID_1934_2003 = 	294277	 if UCID_1934_2003 == 	34561
 replace UCID_1934_2003 = 	33727	 if UCID_1934_2003 == 	40973
 replace UCID_1934_2003 = 	113664	 if UCID_1934_2003 == 	44403
 replace UCID_1934_2003 = 	208663	 if UCID_1934_2003 == 	48338
 replace UCID_1934_2003 = 	48171	 if UCID_1934_2003 == 	48728
 replace UCID_1934_2003 = 	56179	 if UCID_1934_2003 == 	55604
 replace UCID_1934_2003 = 	195144	 if UCID_1934_2003 == 	63818
 replace UCID_1934_2003 = 	374953	 if UCID_1934_2003 == 	65493
 replace UCID_1934_2003 = 	108340	 if UCID_1934_2003 == 	66241
 replace UCID_1934_2003 = 	66350	 if UCID_1934_2003 == 	66582
 replace UCID_1934_2003 = 	288038	 if UCID_1934_2003 == 	69987
 replace UCID_1934_2003 = 	240113	 if UCID_1934_2003 == 	70170
 replace UCID_1934_2003 = 	233061	 if UCID_1934_2003 == 	70896
 replace UCID_1934_2003 = 	52215	 if UCID_1934_2003 == 	73724
 replace UCID_1934_2003 = 	61771	 if UCID_1934_2003 == 	84261
 replace UCID_1934_2003 = 	67065	 if UCID_1934_2003 == 	84719
 replace UCID_1934_2003 = 	57354	 if UCID_1934_2003 == 	86149
 replace UCID_1934_2003 = 	230861	 if UCID_1934_2003 == 	89342
 replace UCID_1934_2003 = 	21892	 if UCID_1934_2003 == 	92129
 replace UCID_1934_2003 = 	54513	 if UCID_1934_2003 == 	95875
 replace UCID_1934_2003 = 	241996	 if UCID_1934_2003 == 	96434
 replace UCID_1934_2003 = 	72486	 if UCID_1934_2003 == 	106125
 replace UCID_1934_2003 = 	120448	 if UCID_1934_2003 == 	108721
 replace UCID_1934_2003 = 	123135	 if UCID_1934_2003 == 	111004
 replace UCID_1934_2003 = 	202414	 if UCID_1934_2003 == 	114834
 replace UCID_1934_2003 = 	202380	 if UCID_1934_2003 == 	116109
 replace UCID_1934_2003 = 	206330	 if UCID_1934_2003 == 	120641
 replace UCID_1934_2003 = 	132284	 if UCID_1934_2003 == 	132337
 replace UCID_1934_2003 = 	163777	 if UCID_1934_2003 == 	134527
 replace UCID_1934_2003 = 	12985	 if UCID_1934_2003 == 	143271
 replace UCID_1934_2003 = 	197476	 if UCID_1934_2003 == 	158221
 replace UCID_1934_2003 = 	171132	 if UCID_1934_2003 == 	159189
 replace UCID_1934_2003 = 	312195	 if UCID_1934_2003 == 	159494
 replace UCID_1934_2003 = 	167647	 if UCID_1934_2003 == 	161836
 replace UCID_1934_2003 = 	171608	 if UCID_1934_2003 == 	162337
 replace UCID_1934_2003 = 	174130	 if UCID_1934_2003 == 	166436
 replace UCID_1934_2003 = 	75133	 if UCID_1934_2003 == 	167145
 replace UCID_1934_2003 = 	148488	 if UCID_1934_2003 == 	167380
 replace UCID_1934_2003 = 	177724	 if UCID_1934_2003 == 	169140
 replace UCID_1934_2003 = 	163806	 if UCID_1934_2003 == 	169890
 replace UCID_1934_2003 = 	248644	 if UCID_1934_2003 == 	180161
 replace UCID_1934_2003 = 	167500	 if UCID_1934_2003 == 	180530
 replace UCID_1934_2003 = 	196229	 if UCID_1934_2003 == 	184091
 replace UCID_1934_2003 = 	164294	 if UCID_1934_2003 == 	185640
 replace UCID_1934_2003 = 	41280	 if UCID_1934_2003 == 	185745
 replace UCID_1934_2003 = 	97205	 if UCID_1934_2003 == 	195749
 replace UCID_1934_2003 = 	341519	 if UCID_1934_2003 == 	200977
 replace UCID_1934_2003 = 	205286	 if UCID_1934_2003 == 	209053
 replace UCID_1934_2003 = 	224654	 if UCID_1934_2003 == 	210295
 replace UCID_1934_2003 = 	211496	 if UCID_1934_2003 == 	225627
 replace UCID_1934_2003 = 	223547	 if UCID_1934_2003 == 	229820
 replace UCID_1934_2003 = 	28020	 if UCID_1934_2003 == 	234019
 replace UCID_1934_2003 = 	151099	 if UCID_1934_2003 == 	247934
 replace UCID_1934_2003 = 	247960	 if UCID_1934_2003 == 	249946
 replace UCID_1934_2003 = 	248499	 if UCID_1934_2003 == 	265643
 replace UCID_1934_2003 = 	76637	 if UCID_1934_2003 == 	266353
 replace UCID_1934_2003 = 	278617	 if UCID_1934_2003 == 	279073
 replace UCID_1934_2003 = 	285934	 if UCID_1934_2003 == 	285575
 replace UCID_1934_2003 = 	291484	 if UCID_1934_2003 == 	289487
 replace UCID_1934_2003 = 	197865	 if UCID_1934_2003 == 	292328
 replace UCID_1934_2003 = 	312337	 if UCID_1934_2003 == 	300537
 replace UCID_1934_2003 = 	23350	 if UCID_1934_2003 == 	301156
 replace UCID_1934_2003 = 	320244	 if UCID_1934_2003 == 	314958
 replace UCID_1934_2003 = 	314701	 if UCID_1934_2003 == 	319827
 replace UCID_1934_2003 = 	17724	 if UCID_1934_2003 == 	320037
 replace UCID_1934_2003 = 	318937	 if UCID_1934_2003 == 	320202
 replace UCID_1934_2003 = 	314733	 if UCID_1934_2003 == 	320388
 replace UCID_1934_2003 = 	320992	 if UCID_1934_2003 == 	320654
 replace UCID_1934_2003 = 	339608	 if UCID_1934_2003 == 	324003
 replace UCID_1934_2003 = 	379006	 if UCID_1934_2003 == 	326900
 replace UCID_1934_2003 = 	226000	 if UCID_1934_2003 == 	332132
 replace UCID_1934_2003 = 	340652	 if UCID_1934_2003 == 	332326
 replace UCID_1934_2003 = 	332975	 if UCID_1934_2003 == 	337689
 replace UCID_1934_2003 = 	203175	 if UCID_1934_2003 == 	342436
 replace UCID_1934_2003 = 	347292	 if UCID_1934_2003 == 	348223
 replace UCID_1934_2003 = 	334718	 if UCID_1934_2003 == 	349167
 replace UCID_1934_2003 = 	242903	 if UCID_1934_2003 == 	349485
 replace UCID_1934_2003 = 	329892	 if UCID_1934_2003 == 	354288
 replace UCID_1934_2003 = 	193818	 if UCID_1934_2003 == 	357844
 replace UCID_1934_2003 = 	47027	 if UCID_1934_2003 == 	383877
 replace UCID_1934_2003 = 	406524	 if UCID_1934_2003 == 	386628
 replace UCID_1934_2003 = 	47223	 if UCID_1934_2003 == 	404264
 replace UCID_1934_2003 = 	250141	 if UCID_1934_2003 == 	410976
 replace UCID_1934_2003 = 	790	 if UCID_1934_2003 == 	425756
 replace UCID_1934_2003 = 	340089	 if UCID_1934_2003 == 	428803

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_61.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_60.dta"

*******************************************************************************
* 1991 & 1992

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_61.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1991 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1991
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1991
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1991_remaining_2.dta", replace
	restore

*Second year 1992 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1992
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1992
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1992_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1991 & 1992
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1991_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1992_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1991_1992_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1991_1992_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1991_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1992_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1991_1992_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1991_1992_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1991_1992_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1991 to 1992
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_61.dta", clear

 replace UCID_1934_2003 = 	390868	 if UCID_1934_2003 == 	2807
 replace UCID_1934_2003 = 	245339	 if UCID_1934_2003 == 	28837
 replace UCID_1934_2003 = 	38828	 if UCID_1934_2003 == 	40639
 replace UCID_1934_2003 = 	200433	 if UCID_1934_2003 == 	41136
 replace UCID_1934_2003 = 	93749	 if UCID_1934_2003 == 	41678
 replace UCID_1934_2003 = 	42986	 if UCID_1934_2003 == 	42985
 replace UCID_1934_2003 = 	246281	 if UCID_1934_2003 == 	50667
 replace UCID_1934_2003 = 	54807	 if UCID_1934_2003 == 	50866
 replace UCID_1934_2003 = 	137983	 if UCID_1934_2003 == 	61261
 replace UCID_1934_2003 = 	246398	 if UCID_1934_2003 == 	68875
 replace UCID_1934_2003 = 	156092	 if UCID_1934_2003 == 	72455
 replace UCID_1934_2003 = 	51242	 if UCID_1934_2003 == 	74136
 replace UCID_1934_2003 = 	62134	 if UCID_1934_2003 == 	74796
 replace UCID_1934_2003 = 	62559	 if UCID_1934_2003 == 	75137
 replace UCID_1934_2003 = 	60479	 if UCID_1934_2003 == 	76417
 replace UCID_1934_2003 = 	72604	 if UCID_1934_2003 == 	77323
 replace UCID_1934_2003 = 	62663	 if UCID_1934_2003 == 	82547
 replace UCID_1934_2003 = 	60375	 if UCID_1934_2003 == 	84018
 replace UCID_1934_2003 = 	233047	 if UCID_1934_2003 == 	84112
 replace UCID_1934_2003 = 	56271	 if UCID_1934_2003 == 	84248
 replace UCID_1934_2003 = 	148061	 if UCID_1934_2003 == 	86528
 replace UCID_1934_2003 = 	50316	 if UCID_1934_2003 == 	89185
 replace UCID_1934_2003 = 	379511	 if UCID_1934_2003 == 	92570
 replace UCID_1934_2003 = 	65543	 if UCID_1934_2003 == 	95246
 replace UCID_1934_2003 = 	60854	 if UCID_1934_2003 == 	97874
 replace UCID_1934_2003 = 	165704	 if UCID_1934_2003 == 	98467
 replace UCID_1934_2003 = 	60357	 if UCID_1934_2003 == 	104259
 replace UCID_1934_2003 = 	119073	 if UCID_1934_2003 == 	112558
 replace UCID_1934_2003 = 	119021	 if UCID_1934_2003 == 	119141
 replace UCID_1934_2003 = 	124236	 if UCID_1934_2003 == 	120055
 replace UCID_1934_2003 = 	114787	 if UCID_1934_2003 == 	123336
 replace UCID_1934_2003 = 	12046	 if UCID_1934_2003 == 	126460
 replace UCID_1934_2003 = 	121973	 if UCID_1934_2003 == 	129232
 replace UCID_1934_2003 = 	144363	 if UCID_1934_2003 == 	130469
 replace UCID_1934_2003 = 	138330	 if UCID_1934_2003 == 	137014
 replace UCID_1934_2003 = 	321307	 if UCID_1934_2003 == 	144555
 replace UCID_1934_2003 = 	141144	 if UCID_1934_2003 == 	145034
 replace UCID_1934_2003 = 	157586	 if UCID_1934_2003 == 	153422
 replace UCID_1934_2003 = 	155542	 if UCID_1934_2003 == 	157388
 replace UCID_1934_2003 = 	156447	 if UCID_1934_2003 == 	157552
 replace UCID_1934_2003 = 	173201	 if UCID_1934_2003 == 	160435
 replace UCID_1934_2003 = 	124543	 if UCID_1934_2003 == 	177451
 replace UCID_1934_2003 = 	150903	 if UCID_1934_2003 == 	179296
 replace UCID_1934_2003 = 	162386	 if UCID_1934_2003 == 	180275
 replace UCID_1934_2003 = 	192803	 if UCID_1934_2003 == 	183881
 replace UCID_1934_2003 = 	90555	 if UCID_1934_2003 == 	185174
 replace UCID_1934_2003 = 	167310	 if UCID_1934_2003 == 	185251
 replace UCID_1934_2003 = 	162402	 if UCID_1934_2003 == 	197982
 replace UCID_1934_2003 = 	197255	 if UCID_1934_2003 == 	202510
 replace UCID_1934_2003 = 	203563	 if UCID_1934_2003 == 	207141
 replace UCID_1934_2003 = 	225779	 if UCID_1934_2003 == 	225752
 replace UCID_1934_2003 = 	234768	 if UCID_1934_2003 == 	237534
 replace UCID_1934_2003 = 	243177	 if UCID_1934_2003 == 	245573
 replace UCID_1934_2003 = 	169523	 if UCID_1934_2003 == 	246009
 replace UCID_1934_2003 = 	85911	 if UCID_1934_2003 == 	248621
 replace UCID_1934_2003 = 	390549	 if UCID_1934_2003 == 	260933
 replace UCID_1934_2003 = 	58854	 if UCID_1934_2003 == 	265573
 replace UCID_1934_2003 = 	217990	 if UCID_1934_2003 == 	265957
 replace UCID_1934_2003 = 	227672	 if UCID_1934_2003 == 	266442
 replace UCID_1934_2003 = 	82333	 if UCID_1934_2003 == 	267065
 replace UCID_1934_2003 = 	279406	 if UCID_1934_2003 == 	273862
 replace UCID_1934_2003 = 	275479	 if UCID_1934_2003 == 	275293
 replace UCID_1934_2003 = 	271239	 if UCID_1934_2003 == 	277814
 replace UCID_1934_2003 = 	278279	 if UCID_1934_2003 == 	278156
 replace UCID_1934_2003 = 	235715	 if UCID_1934_2003 == 	282514
 replace UCID_1934_2003 = 	283757	 if UCID_1934_2003 == 	283956
 replace UCID_1934_2003 = 	34611	 if UCID_1934_2003 == 	289475
 replace UCID_1934_2003 = 	147737	 if UCID_1934_2003 == 	289827
 replace UCID_1934_2003 = 	150827	 if UCID_1934_2003 == 	294139
 replace UCID_1934_2003 = 	147726	 if UCID_1934_2003 == 	298341
 replace UCID_1934_2003 = 	54705	 if UCID_1934_2003 == 	308231
 replace UCID_1934_2003 = 	167605	 if UCID_1934_2003 == 	313797
 replace UCID_1934_2003 = 	156515	 if UCID_1934_2003 == 	322118
 replace UCID_1934_2003 = 	203674	 if UCID_1934_2003 == 	325316
 replace UCID_1934_2003 = 	339853	 if UCID_1934_2003 == 	328026
 replace UCID_1934_2003 = 	350444	 if UCID_1934_2003 == 	335720
 replace UCID_1934_2003 = 	322826	 if UCID_1934_2003 == 	348018
 replace UCID_1934_2003 = 	165683	 if UCID_1934_2003 == 	348398
 replace UCID_1934_2003 = 	354359	 if UCID_1934_2003 == 	351551
 replace UCID_1934_2003 = 	361058	 if UCID_1934_2003 == 	354543
 replace UCID_1934_2003 = 	369195	 if UCID_1934_2003 == 	376362
 replace UCID_1934_2003 = 	369002	 if UCID_1934_2003 == 	377779
 replace UCID_1934_2003 = 	184764	 if UCID_1934_2003 == 	386453
 replace UCID_1934_2003 = 	369354	 if UCID_1934_2003 == 	386675
 replace UCID_1934_2003 = 	246755	 if UCID_1934_2003 == 	413459
 replace UCID_1934_2003 = 	431987	 if UCID_1934_2003 == 	422778

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_62.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_61.dta"

*******************************************************************************
* 1992 & 1993

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_62.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1992 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1992
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1992
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1992_remaining_2.dta", replace
	restore

*Second year 1993 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1993
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1993
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1993_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1992 & 1993
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1992_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1993_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1992_1993_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1992_1993_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1992_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1993_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1992_1993_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1992_1993_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1992_1993_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1992 to 1993
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_62.dta", clear

 replace UCID_1934_2003 = 	1779	 if UCID_1934_2003 == 	1210
 replace UCID_1934_2003 = 	274833	 if UCID_1934_2003 == 	2193
 replace UCID_1934_2003 = 	12891	 if UCID_1934_2003 == 	13279
 replace UCID_1934_2003 = 	135989	 if UCID_1934_2003 == 	21451
 replace UCID_1934_2003 = 	343235	 if UCID_1934_2003 == 	60693
 replace UCID_1934_2003 = 	62825	 if UCID_1934_2003 == 	64455
 replace UCID_1934_2003 = 	67559	 if UCID_1934_2003 == 	69250
 replace UCID_1934_2003 = 	179778	 if UCID_1934_2003 == 	76874
 replace UCID_1934_2003 = 	413362	 if UCID_1934_2003 == 	79863
 replace UCID_1934_2003 = 	50345	 if UCID_1934_2003 == 	83834
 replace UCID_1934_2003 = 	194866	 if UCID_1934_2003 == 	85656
 replace UCID_1934_2003 = 	201174	 if UCID_1934_2003 == 	86137
 replace UCID_1934_2003 = 	185917	 if UCID_1934_2003 == 	101768
 replace UCID_1934_2003 = 	186665	 if UCID_1934_2003 == 	103861
 replace UCID_1934_2003 = 	11656	 if UCID_1934_2003 == 	116633
 replace UCID_1934_2003 = 	132994	 if UCID_1934_2003 == 	133588
 replace UCID_1934_2003 = 	60318	 if UCID_1934_2003 == 	140014
 replace UCID_1934_2003 = 	138380	 if UCID_1934_2003 == 	142095
 replace UCID_1934_2003 = 	152267	 if UCID_1934_2003 == 	152790
 replace UCID_1934_2003 = 	155949	 if UCID_1934_2003 == 	156590
 replace UCID_1934_2003 = 	193628	 if UCID_1934_2003 == 	158127
 replace UCID_1934_2003 = 	175782	 if UCID_1934_2003 == 	163256
 replace UCID_1934_2003 = 	100816	 if UCID_1934_2003 == 	164567
 replace UCID_1934_2003 = 	23339	 if UCID_1934_2003 == 	167062
 replace UCID_1934_2003 = 	341703	 if UCID_1934_2003 == 	175052
 replace UCID_1934_2003 = 	366878	 if UCID_1934_2003 == 	185825
 replace UCID_1934_2003 = 	133465	 if UCID_1934_2003 == 	207260
 replace UCID_1934_2003 = 	226678	 if UCID_1934_2003 == 	213137
 replace UCID_1934_2003 = 	229864	 if UCID_1934_2003 == 	242424
 replace UCID_1934_2003 = 	50683	 if UCID_1934_2003 == 	243609
 replace UCID_1934_2003 = 	70652	 if UCID_1934_2003 == 	244149
 replace UCID_1934_2003 = 	262375	 if UCID_1934_2003 == 	247596
 replace UCID_1934_2003 = 	249486	 if UCID_1934_2003 == 	249348
 replace UCID_1934_2003 = 	422308	 if UCID_1934_2003 == 	267240
 replace UCID_1934_2003 = 	236349	 if UCID_1934_2003 == 	269329
 replace UCID_1934_2003 = 	375630	 if UCID_1934_2003 == 	286929
 replace UCID_1934_2003 = 	33808	 if UCID_1934_2003 == 	288302
 replace UCID_1934_2003 = 	136259	 if UCID_1934_2003 == 	288511
 replace UCID_1934_2003 = 	295554	 if UCID_1934_2003 == 	297436
 replace UCID_1934_2003 = 	306639	 if UCID_1934_2003 == 	315118
 replace UCID_1934_2003 = 	320226	 if UCID_1934_2003 == 	317375
 replace UCID_1934_2003 = 	320231	 if UCID_1934_2003 == 	318393
 replace UCID_1934_2003 = 	83385	 if UCID_1934_2003 == 	324786
 replace UCID_1934_2003 = 	74140	 if UCID_1934_2003 == 	325439
 replace UCID_1934_2003 = 	38635	 if UCID_1934_2003 == 	327139
 replace UCID_1934_2003 = 	337407	 if UCID_1934_2003 == 	333552
 replace UCID_1934_2003 = 	339236	 if UCID_1934_2003 == 	335322
 replace UCID_1934_2003 = 	330415	 if UCID_1934_2003 == 	337038
 replace UCID_1934_2003 = 	3207	 if UCID_1934_2003 == 	340686
 replace UCID_1934_2003 = 	242768	 if UCID_1934_2003 == 	355554
 replace UCID_1934_2003 = 	90959	 if UCID_1934_2003 == 	365915
 replace UCID_1934_2003 = 	370246	 if UCID_1934_2003 == 	389136
 replace UCID_1934_2003 = 	403891	 if UCID_1934_2003 == 	395523
 replace UCID_1934_2003 = 	196356	 if UCID_1934_2003 == 	396496
 replace UCID_1934_2003 = 	316509	 if UCID_1934_2003 == 	397122
 replace UCID_1934_2003 = 	136952	 if UCID_1934_2003 == 	414252

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_63.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_62.dta"

*******************************************************************************
* 1993 & 1994

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_63.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1993 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1993
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1993
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1993_remaining_2.dta", replace
	restore

*Second year 1994 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1994
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1994
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1994_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1993 & 1994
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1993_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1994_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1993_1994_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1993_1994_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1993_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1994_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1993_1994_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1993_1994_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1993_1994_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1993 to 1994
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_63.dta", clear

 replace UCID_1934_2003 = 	62432	 if UCID_1934_2003 == 	6466
 replace UCID_1934_2003 = 	1595	 if UCID_1934_2003 == 	6889
 replace UCID_1934_2003 = 	128004	 if UCID_1934_2003 == 	12098
 replace UCID_1934_2003 = 	15433	 if UCID_1934_2003 == 	15018
 replace UCID_1934_2003 = 	157204	 if UCID_1934_2003 == 	15166
 replace UCID_1934_2003 = 	370388	 if UCID_1934_2003 == 	26367
 replace UCID_1934_2003 = 	416774	 if UCID_1934_2003 == 	34468
 replace UCID_1934_2003 = 	37559	 if UCID_1934_2003 == 	36333
 replace UCID_1934_2003 = 	37005	 if UCID_1934_2003 == 	38433
 replace UCID_1934_2003 = 	52609	 if UCID_1934_2003 == 	50197
 replace UCID_1934_2003 = 	53183	 if UCID_1934_2003 == 	50737
 replace UCID_1934_2003 = 	61900	 if UCID_1934_2003 == 	62969
 replace UCID_1934_2003 = 	70698	 if UCID_1934_2003 == 	63065
 replace UCID_1934_2003 = 	51452	 if UCID_1934_2003 == 	69961
 replace UCID_1934_2003 = 	49964	 if UCID_1934_2003 == 	71430
 replace UCID_1934_2003 = 	281591	 if UCID_1934_2003 == 	71946
 replace UCID_1934_2003 = 	175166	 if UCID_1934_2003 == 	74553
 replace UCID_1934_2003 = 	57717	 if UCID_1934_2003 == 	77721
 replace UCID_1934_2003 = 	236276	 if UCID_1934_2003 == 	81944
 replace UCID_1934_2003 = 	348210	 if UCID_1934_2003 == 	87651
 replace UCID_1934_2003 = 	57786	 if UCID_1934_2003 == 	90472
 replace UCID_1934_2003 = 	53859	 if UCID_1934_2003 == 	95854
 replace UCID_1934_2003 = 	73401	 if UCID_1934_2003 == 	97014
 replace UCID_1934_2003 = 	49516	 if UCID_1934_2003 == 	98033
 replace UCID_1934_2003 = 	138707	 if UCID_1934_2003 == 	103248
 replace UCID_1934_2003 = 	67295	 if UCID_1934_2003 == 	108104
 replace UCID_1934_2003 = 	323398	 if UCID_1934_2003 == 	114929
 replace UCID_1934_2003 = 	10658	 if UCID_1934_2003 == 	118152
 replace UCID_1934_2003 = 	110883	 if UCID_1934_2003 == 	118682
 replace UCID_1934_2003 = 	123245	 if UCID_1934_2003 == 	124079
 replace UCID_1934_2003 = 	118725	 if UCID_1934_2003 == 	130101
 replace UCID_1934_2003 = 	134319	 if UCID_1934_2003 == 	133824
 replace UCID_1934_2003 = 	146117	 if UCID_1934_2003 == 	135660
 replace UCID_1934_2003 = 	141712	 if UCID_1934_2003 == 	136313
 replace UCID_1934_2003 = 	136839	 if UCID_1934_2003 == 	137380
 replace UCID_1934_2003 = 	148257	 if UCID_1934_2003 == 	152764
 replace UCID_1934_2003 = 	184884	 if UCID_1934_2003 == 	152815
 replace UCID_1934_2003 = 	318936	 if UCID_1934_2003 == 	172109
 replace UCID_1934_2003 = 	165016	 if UCID_1934_2003 == 	175401
 replace UCID_1934_2003 = 	16616	 if UCID_1934_2003 == 	180821
 replace UCID_1934_2003 = 	162938	 if UCID_1934_2003 == 	182917
 replace UCID_1934_2003 = 	111759	 if UCID_1934_2003 == 	192271
 replace UCID_1934_2003 = 	166707	 if UCID_1934_2003 == 	196097
 replace UCID_1934_2003 = 	186545	 if UCID_1934_2003 == 	207958
 replace UCID_1934_2003 = 	205681	 if UCID_1934_2003 == 	209051
 replace UCID_1934_2003 = 	228029	 if UCID_1934_2003 == 	227284
 replace UCID_1934_2003 = 	228234	 if UCID_1934_2003 == 	228124
 replace UCID_1934_2003 = 	229909	 if UCID_1934_2003 == 	228730
 replace UCID_1934_2003 = 	166242	 if UCID_1934_2003 == 	231826
 replace UCID_1934_2003 = 	238407	 if UCID_1934_2003 == 	241649
 replace UCID_1934_2003 = 	244278	 if UCID_1934_2003 == 	244098
 replace UCID_1934_2003 = 	250276	 if UCID_1934_2003 == 	250152
 replace UCID_1934_2003 = 	272615	 if UCID_1934_2003 == 	273149
 replace UCID_1934_2003 = 	273930	 if UCID_1934_2003 == 	274009
 replace UCID_1934_2003 = 	199573	 if UCID_1934_2003 == 	284418
 replace UCID_1934_2003 = 	300359	 if UCID_1934_2003 == 	293987
 replace UCID_1934_2003 = 	66933	 if UCID_1934_2003 == 	294829
 replace UCID_1934_2003 = 	149896	 if UCID_1934_2003 == 	305084
 replace UCID_1934_2003 = 	290616	 if UCID_1934_2003 == 	309721
 replace UCID_1934_2003 = 	313686	 if UCID_1934_2003 == 	311501
 replace UCID_1934_2003 = 	159107	 if UCID_1934_2003 == 	313772
 replace UCID_1934_2003 = 	316570	 if UCID_1934_2003 == 	320502
 replace UCID_1934_2003 = 	150396	 if UCID_1934_2003 == 	321116
 replace UCID_1934_2003 = 	344074	 if UCID_1934_2003 == 	323165
 replace UCID_1934_2003 = 	337262	 if UCID_1934_2003 == 	329593
 replace UCID_1934_2003 = 	339133	 if UCID_1934_2003 == 	335363
 replace UCID_1934_2003 = 	322886	 if UCID_1934_2003 == 	348725
 replace UCID_1934_2003 = 	131524	 if UCID_1934_2003 == 	349472
 replace UCID_1934_2003 = 	341974	 if UCID_1934_2003 == 	350037
 replace UCID_1934_2003 = 	351378	 if UCID_1934_2003 == 	360364
 replace UCID_1934_2003 = 	348126	 if UCID_1934_2003 == 	366026
 replace UCID_1934_2003 = 	352182	 if UCID_1934_2003 == 	370578
 replace UCID_1934_2003 = 	136969	 if UCID_1934_2003 == 	380275
 replace UCID_1934_2003 = 	276177	 if UCID_1934_2003 == 	414380
 replace UCID_1934_2003 = 	288921	 if UCID_1934_2003 == 	415792

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_64.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_63.dta"

*******************************************************************************
* 1994 & 1995

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_64.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1994 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1994
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1994
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1994_remaining_2.dta", replace
	restore

*Second year 1995 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1995
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1995
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1995_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1994 & 1995
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1994_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1995_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1994_1995_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1994_1995_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1994_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1995_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1994_1995_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1994_1995_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1994_1995_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1994 to 1995
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_64.dta", clear

 replace UCID_1934_2003 = 	3795	 if UCID_1934_2003 == 	1888
 replace UCID_1934_2003 = 	71839	 if UCID_1934_2003 == 	4628
 replace UCID_1934_2003 = 	131752	 if UCID_1934_2003 == 	9752
 replace UCID_1934_2003 = 	164141	 if UCID_1934_2003 == 	17645
 replace UCID_1934_2003 = 	41532	 if UCID_1934_2003 == 	21549
 replace UCID_1934_2003 = 	427615	 if UCID_1934_2003 == 	28756
 replace UCID_1934_2003 = 	292336	 if UCID_1934_2003 == 	31197
 replace UCID_1934_2003 = 	33028	 if UCID_1934_2003 == 	33051
 replace UCID_1934_2003 = 	37815	 if UCID_1934_2003 == 	36303
 replace UCID_1934_2003 = 	53905	 if UCID_1934_2003 == 	50608
 replace UCID_1934_2003 = 	53862	 if UCID_1934_2003 == 	54467
 replace UCID_1934_2003 = 	52619	 if UCID_1934_2003 == 	56235
 replace UCID_1934_2003 = 	180277	 if UCID_1934_2003 == 	59892
 replace UCID_1934_2003 = 	65898	 if UCID_1934_2003 == 	65113
 replace UCID_1934_2003 = 	63282	 if UCID_1934_2003 == 	67041
 replace UCID_1934_2003 = 	52872	 if UCID_1934_2003 == 	67196
 replace UCID_1934_2003 = 	166777	 if UCID_1934_2003 == 	76067
 replace UCID_1934_2003 = 	203481	 if UCID_1934_2003 == 	83345
 replace UCID_1934_2003 = 	63958	 if UCID_1934_2003 == 	88191
 replace UCID_1934_2003 = 	243860	 if UCID_1934_2003 == 	92210
 replace UCID_1934_2003 = 	162793	 if UCID_1934_2003 == 	97040
 replace UCID_1934_2003 = 	402399	 if UCID_1934_2003 == 	98619
 replace UCID_1934_2003 = 	207623	 if UCID_1934_2003 == 	98709
 replace UCID_1934_2003 = 	336637	 if UCID_1934_2003 == 	98828
 replace UCID_1934_2003 = 	62214	 if UCID_1934_2003 == 	99750
 replace UCID_1934_2003 = 	73125	 if UCID_1934_2003 == 	100523
 replace UCID_1934_2003 = 	51245	 if UCID_1934_2003 == 	103780
 replace UCID_1934_2003 = 	114731	 if UCID_1934_2003 == 	117592
 replace UCID_1934_2003 = 	112759	 if UCID_1934_2003 == 	121381
 replace UCID_1934_2003 = 	277666	 if UCID_1934_2003 == 	125044
 replace UCID_1934_2003 = 	126786	 if UCID_1934_2003 == 	125936
 replace UCID_1934_2003 = 	133214	 if UCID_1934_2003 == 	130284
 replace UCID_1934_2003 = 	133018	 if UCID_1934_2003 == 	133567
 replace UCID_1934_2003 = 	157282	 if UCID_1934_2003 == 	136922
 replace UCID_1934_2003 = 	147916	 if UCID_1934_2003 == 	151025
 replace UCID_1934_2003 = 	167267	 if UCID_1934_2003 == 	153643
 replace UCID_1934_2003 = 	99637	 if UCID_1934_2003 == 	160204
 replace UCID_1934_2003 = 	167938	 if UCID_1934_2003 == 	163816
 replace UCID_1934_2003 = 	341830	 if UCID_1934_2003 == 	170337
 replace UCID_1934_2003 = 	148116	 if UCID_1934_2003 == 	171864
 replace UCID_1934_2003 = 	79751	 if UCID_1934_2003 == 	172462
 replace UCID_1934_2003 = 	351998	 if UCID_1934_2003 == 	173862
 replace UCID_1934_2003 = 	251958	 if UCID_1934_2003 == 	180890
 replace UCID_1934_2003 = 	224970	 if UCID_1934_2003 == 	204467
 replace UCID_1934_2003 = 	52016	 if UCID_1934_2003 == 	225001
 replace UCID_1934_2003 = 	228463	 if UCID_1934_2003 == 	229016
 replace UCID_1934_2003 = 	231363	 if UCID_1934_2003 == 	231130
 replace UCID_1934_2003 = 	233231	 if UCID_1934_2003 == 	231769
 replace UCID_1934_2003 = 	252748	 if UCID_1934_2003 == 	235156
 replace UCID_1934_2003 = 	83250	 if UCID_1934_2003 == 	248107
 replace UCID_1934_2003 = 	268045	 if UCID_1934_2003 == 	265865
 replace UCID_1934_2003 = 	169239	 if UCID_1934_2003 == 	269686
 replace UCID_1934_2003 = 	163461	 if UCID_1934_2003 == 	272480
 replace UCID_1934_2003 = 	274785	 if UCID_1934_2003 == 	274923
 replace UCID_1934_2003 = 	303642	 if UCID_1934_2003 == 	275217
 replace UCID_1934_2003 = 	271866	 if UCID_1934_2003 == 	279141
 replace UCID_1934_2003 = 	71630	 if UCID_1934_2003 == 	279990
 replace UCID_1934_2003 = 	267847	 if UCID_1934_2003 == 	280262
 replace UCID_1934_2003 = 	67301	 if UCID_1934_2003 == 	281382
 replace UCID_1934_2003 = 	270330	 if UCID_1934_2003 == 	281622
 replace UCID_1934_2003 = 	299676	 if UCID_1934_2003 == 	291433
 replace UCID_1934_2003 = 	48412	 if UCID_1934_2003 == 	292431
 replace UCID_1934_2003 = 	349399	 if UCID_1934_2003 == 	300443
 replace UCID_1934_2003 = 	321126	 if UCID_1934_2003 == 	306856
 replace UCID_1934_2003 = 	144555	 if UCID_1934_2003 == 	321307
 replace UCID_1934_2003 = 	339237	 if UCID_1934_2003 == 	339580
 replace UCID_1934_2003 = 	342121	 if UCID_1934_2003 == 	345371
 replace UCID_1934_2003 = 	190337	 if UCID_1934_2003 == 	348393
 replace UCID_1934_2003 = 	349681	 if UCID_1934_2003 == 	349678
 replace UCID_1934_2003 = 	42286	 if UCID_1934_2003 == 	361209
 replace UCID_1934_2003 = 	367523	 if UCID_1934_2003 == 	366141
 replace UCID_1934_2003 = 	366129	 if UCID_1934_2003 == 	367094
 replace UCID_1934_2003 = 	154834	 if UCID_1934_2003 == 	368216
 replace UCID_1934_2003 = 	242044	 if UCID_1934_2003 == 	373338
 replace UCID_1934_2003 = 	201924	 if UCID_1934_2003 == 	395443
 replace UCID_1934_2003 = 	330249	 if UCID_1934_2003 == 	397844
 replace UCID_1934_2003 = 	405489	 if UCID_1934_2003 == 	402210
 replace UCID_1934_2003 = 	251400	 if UCID_1934_2003 == 	413968
 replace UCID_1934_2003 = 	198644	 if UCID_1934_2003 == 	425989

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_65.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_64.dta"

*******************************************************************************
* 1995 & 1996

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_65.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1995 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1995
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1995
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1995_remaining_2.dta", replace
	restore

*Second year 1996 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1996
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1996
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1996_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1995 & 1996
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1995_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1996_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1995_1996_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1995_1996_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1995_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1996_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1995_1996_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1995_1996_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1995_1996_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1995 to 1996
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_65.dta", clear

 replace UCID_1934_2003 = 	4181	 if UCID_1934_2003 == 	604
 replace UCID_1934_2003 = 	29620	 if UCID_1934_2003 == 	7212
 replace UCID_1934_2003 = 	12430	 if UCID_1934_2003 == 	9558
 replace UCID_1934_2003 = 	11137	 if UCID_1934_2003 == 	9651
 replace UCID_1934_2003 = 	13976	 if UCID_1934_2003 == 	14038
 replace UCID_1934_2003 = 	37109	 if UCID_1934_2003 == 	16309
 replace UCID_1934_2003 = 	17165	 if UCID_1934_2003 == 	17302
 replace UCID_1934_2003 = 	17288	 if UCID_1934_2003 == 	20148
 replace UCID_1934_2003 = 	19292	 if UCID_1934_2003 == 	29358
 replace UCID_1934_2003 = 	37577	 if UCID_1934_2003 == 	37169
 replace UCID_1934_2003 = 	354796	 if UCID_1934_2003 == 	41809
 replace UCID_1934_2003 = 	361617	 if UCID_1934_2003 == 	42720
 replace UCID_1934_2003 = 	404897	 if UCID_1934_2003 == 	43934
 replace UCID_1934_2003 = 	43128	 if UCID_1934_2003 == 	45359
 replace UCID_1934_2003 = 	52107	 if UCID_1934_2003 == 	45599
 replace UCID_1934_2003 = 	21610	 if UCID_1934_2003 == 	47460
 replace UCID_1934_2003 = 	428001	 if UCID_1934_2003 == 	48037
 replace UCID_1934_2003 = 	55233	 if UCID_1934_2003 == 	52163
 replace UCID_1934_2003 = 	102981	 if UCID_1934_2003 == 	53077
 replace UCID_1934_2003 = 	79556	 if UCID_1934_2003 == 	59108
 replace UCID_1934_2003 = 	82547	 if UCID_1934_2003 == 	62663
 replace UCID_1934_2003 = 	64764	 if UCID_1934_2003 == 	62825
 replace UCID_1934_2003 = 	108002	 if UCID_1934_2003 == 	62847
 replace UCID_1934_2003 = 	807	 if UCID_1934_2003 == 	64007
 replace UCID_1934_2003 = 	107965	 if UCID_1934_2003 == 	66803
 replace UCID_1934_2003 = 	68189	 if UCID_1934_2003 == 	67754
 replace UCID_1934_2003 = 	80924	 if UCID_1934_2003 == 	69959
 replace UCID_1934_2003 = 	179317	 if UCID_1934_2003 == 	71275
 replace UCID_1934_2003 = 	53294	 if UCID_1934_2003 == 	73172
 replace UCID_1934_2003 = 	73092	 if UCID_1934_2003 == 	75591
 replace UCID_1934_2003 = 	66252	 if UCID_1934_2003 == 	84306
 replace UCID_1934_2003 = 	72202	 if UCID_1934_2003 == 	85252
 replace UCID_1934_2003 = 	62865	 if UCID_1934_2003 == 	86922
 replace UCID_1934_2003 = 	71245	 if UCID_1934_2003 == 	87823
 replace UCID_1934_2003 = 	388997	 if UCID_1934_2003 == 	88868
 replace UCID_1934_2003 = 	330314	 if UCID_1934_2003 == 	89932
 replace UCID_1934_2003 = 	58178	 if UCID_1934_2003 == 	91644
 replace UCID_1934_2003 = 	62875	 if UCID_1934_2003 == 	91904
 replace UCID_1934_2003 = 	286710	 if UCID_1934_2003 == 	93637
 replace UCID_1934_2003 = 	55790	 if UCID_1934_2003 == 	94543
 replace UCID_1934_2003 = 	185495	 if UCID_1934_2003 == 	101270
 replace UCID_1934_2003 = 	284216	 if UCID_1934_2003 == 	107389
 replace UCID_1934_2003 = 	109250	 if UCID_1934_2003 == 	109516
 replace UCID_1934_2003 = 	123259	 if UCID_1934_2003 == 	112665
 replace UCID_1934_2003 = 	12170	 if UCID_1934_2003 == 	112750
 replace UCID_1934_2003 = 	11275	 if UCID_1934_2003 == 	113613
 replace UCID_1934_2003 = 	118258	 if UCID_1934_2003 == 	114179
 replace UCID_1934_2003 = 	170749	 if UCID_1934_2003 == 	115241
 replace UCID_1934_2003 = 	426122	 if UCID_1934_2003 == 	117184
 replace UCID_1934_2003 = 	125611	 if UCID_1934_2003 == 	118762
 replace UCID_1934_2003 = 	11316	 if UCID_1934_2003 == 	119193
 replace UCID_1934_2003 = 	145002	 if UCID_1934_2003 == 	144659
 replace UCID_1934_2003 = 	153226	 if UCID_1934_2003 == 	153757
 replace UCID_1934_2003 = 	210774	 if UCID_1934_2003 == 	154947
 replace UCID_1934_2003 = 	156403	 if UCID_1934_2003 == 	155086
 replace UCID_1934_2003 = 	166744	 if UCID_1934_2003 == 	164022
 replace UCID_1934_2003 = 	163647	 if UCID_1934_2003 == 	164974
 replace UCID_1934_2003 = 	163030	 if UCID_1934_2003 == 	166295
 replace UCID_1934_2003 = 	168713	 if UCID_1934_2003 == 	170683
 replace UCID_1934_2003 = 	163232	 if UCID_1934_2003 == 	172290
 replace UCID_1934_2003 = 	165960	 if UCID_1934_2003 == 	174506
 replace UCID_1934_2003 = 	165261	 if UCID_1934_2003 == 	178291
 replace UCID_1934_2003 = 	167640	 if UCID_1934_2003 == 	184279
 replace UCID_1934_2003 = 	385383	 if UCID_1934_2003 == 	187265
 replace UCID_1934_2003 = 	292068	 if UCID_1934_2003 == 	187774
 replace UCID_1934_2003 = 	191732	 if UCID_1934_2003 == 	194562
 replace UCID_1934_2003 = 	204528	 if UCID_1934_2003 == 	212811
 replace UCID_1934_2003 = 	227579	 if UCID_1934_2003 == 	216544
 replace UCID_1934_2003 = 	167118	 if UCID_1934_2003 == 	217352
 replace UCID_1934_2003 = 	473368	 if UCID_1934_2003 == 	220632
 replace UCID_1934_2003 = 	170247	 if UCID_1934_2003 == 	221137
 replace UCID_1934_2003 = 	227637	 if UCID_1934_2003 == 	224791
 replace UCID_1934_2003 = 	405828	 if UCID_1934_2003 == 	235181
 replace UCID_1934_2003 = 	219985	 if UCID_1934_2003 == 	238130
 replace UCID_1934_2003 = 	241483	 if UCID_1934_2003 == 	239837
 replace UCID_1934_2003 = 	246375	 if UCID_1934_2003 == 	245977
 replace UCID_1934_2003 = 	251016	 if UCID_1934_2003 == 	249797
 replace UCID_1934_2003 = 	407320	 if UCID_1934_2003 == 	256023
 replace UCID_1934_2003 = 	64310	 if UCID_1934_2003 == 	260862
 replace UCID_1934_2003 = 	227288	 if UCID_1934_2003 == 	275845
 replace UCID_1934_2003 = 	272434	 if UCID_1934_2003 == 	281325
 replace UCID_1934_2003 = 	251637	 if UCID_1934_2003 == 	287582
 replace UCID_1934_2003 = 	314731	 if UCID_1934_2003 == 	306035
 replace UCID_1934_2003 = 	308667	 if UCID_1934_2003 == 	317230
 replace UCID_1934_2003 = 	338434	 if UCID_1934_2003 == 	327172
 replace UCID_1934_2003 = 	381568	 if UCID_1934_2003 == 	327321
 replace UCID_1934_2003 = 	342257	 if UCID_1934_2003 == 	328918
 replace UCID_1934_2003 = 	337921	 if UCID_1934_2003 == 	330381
 replace UCID_1934_2003 = 	367266	 if UCID_1934_2003 == 	337569
 replace UCID_1934_2003 = 	353159	 if UCID_1934_2003 == 	356881
 replace UCID_1934_2003 = 	356007	 if UCID_1934_2003 == 	359957
 replace UCID_1934_2003 = 	204641	 if UCID_1934_2003 == 	368576
 replace UCID_1934_2003 = 	369218	 if UCID_1934_2003 == 	394376
 replace UCID_1934_2003 = 	154796	 if UCID_1934_2003 == 	410308
 replace UCID_1934_2003 = 	327511	 if UCID_1934_2003 == 	416948
 replace UCID_1934_2003 = 	81346	 if UCID_1934_2003 == 	425469

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_66.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_65.dta"

*******************************************************************************
* 1996 & 1997

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_66.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1996 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1996
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1996
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1996_remaining_2.dta", replace
	restore

*Second year 1997 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1997
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1997
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1997_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1996 & 1997
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1996_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1997_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1996_1997_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1996_1997_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1996_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1997_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1996_1997_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1996_1997_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1996_1997_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1996 to 1997
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_66.dta", clear

 replace UCID_1934_2003 = 	1751	 if UCID_1934_2003 == 	294
 replace UCID_1934_2003 = 	14377	 if UCID_1934_2003 == 	2303
 replace UCID_1934_2003 = 	8702	 if UCID_1934_2003 == 	2320
 replace UCID_1934_2003 = 	170535	 if UCID_1934_2003 == 	2751
 replace UCID_1934_2003 = 	2321	 if UCID_1934_2003 == 	8734
 replace UCID_1934_2003 = 	128037	 if UCID_1934_2003 == 	9902
 replace UCID_1934_2003 = 	118037	 if UCID_1934_2003 == 	10984
 replace UCID_1934_2003 = 	9382	 if UCID_1934_2003 == 	12258
 replace UCID_1934_2003 = 	1381	 if UCID_1934_2003 == 	12272
 replace UCID_1934_2003 = 	70121	 if UCID_1934_2003 == 	16935
 replace UCID_1934_2003 = 	229957	 if UCID_1934_2003 == 	24419
 replace UCID_1934_2003 = 	244021	 if UCID_1934_2003 == 	29335
 replace UCID_1934_2003 = 	32121	 if UCID_1934_2003 == 	32461
 replace UCID_1934_2003 = 	32562	 if UCID_1934_2003 == 	32637
 replace UCID_1934_2003 = 	34110	 if UCID_1934_2003 == 	33645
 replace UCID_1934_2003 = 	33589	 if UCID_1934_2003 == 	33777
 replace UCID_1934_2003 = 	34418	 if UCID_1934_2003 == 	38521
 replace UCID_1934_2003 = 	339221	 if UCID_1934_2003 == 	41315
 replace UCID_1934_2003 = 	361628	 if UCID_1934_2003 == 	48073
 replace UCID_1934_2003 = 	74038	 if UCID_1934_2003 == 	50998
 replace UCID_1934_2003 = 	58729	 if UCID_1934_2003 == 	57354
 replace UCID_1934_2003 = 	79488	 if UCID_1934_2003 == 	57877
 replace UCID_1934_2003 = 	410028	 if UCID_1934_2003 == 	59253
 replace UCID_1934_2003 = 	281371	 if UCID_1934_2003 == 	67301
 replace UCID_1934_2003 = 	73049	 if UCID_1934_2003 == 	73080
 replace UCID_1934_2003 = 	72456	 if UCID_1934_2003 == 	91625
 replace UCID_1934_2003 = 	171060	 if UCID_1934_2003 == 	96181
 replace UCID_1934_2003 = 	65103	 if UCID_1934_2003 == 	97730
 replace UCID_1934_2003 = 	117716	 if UCID_1934_2003 == 	112379
 replace UCID_1934_2003 = 	118545	 if UCID_1934_2003 == 	112617
 replace UCID_1934_2003 = 	117544	 if UCID_1934_2003 == 	115222
 replace UCID_1934_2003 = 	111046	 if UCID_1934_2003 == 	117317
 replace UCID_1934_2003 = 	110062	 if UCID_1934_2003 == 	121817
 replace UCID_1934_2003 = 	125463	 if UCID_1934_2003 == 	128198
 replace UCID_1934_2003 = 	130785	 if UCID_1934_2003 == 	130721
 replace UCID_1934_2003 = 	119390	 if UCID_1934_2003 == 	133687
 replace UCID_1934_2003 = 	66927	 if UCID_1934_2003 == 	142187
 replace UCID_1934_2003 = 	145831	 if UCID_1934_2003 == 	144801
 replace UCID_1934_2003 = 	144885	 if UCID_1934_2003 == 	145568
 replace UCID_1934_2003 = 	139333	 if UCID_1934_2003 == 	145898
 replace UCID_1934_2003 = 	77285	 if UCID_1934_2003 == 	146466
 replace UCID_1934_2003 = 	298341	 if UCID_1934_2003 == 	147726
 replace UCID_1934_2003 = 	63266	 if UCID_1934_2003 == 	158147
 replace UCID_1934_2003 = 	31092	 if UCID_1934_2003 == 	159615
 replace UCID_1934_2003 = 	270597	 if UCID_1934_2003 == 	164013
 replace UCID_1934_2003 = 	157300	 if UCID_1934_2003 == 	170326
 replace UCID_1934_2003 = 	17337	 if UCID_1934_2003 == 	174372
 replace UCID_1934_2003 = 	162455	 if UCID_1934_2003 == 	185974
 replace UCID_1934_2003 = 	165227	 if UCID_1934_2003 == 	193992
 replace UCID_1934_2003 = 	23294	 if UCID_1934_2003 == 	194630
 replace UCID_1934_2003 = 	391085	 if UCID_1934_2003 == 	199164
 replace UCID_1934_2003 = 	323318	 if UCID_1934_2003 == 	201014
 replace UCID_1934_2003 = 	200839	 if UCID_1934_2003 == 	202784
 replace UCID_1934_2003 = 	204994	 if UCID_1934_2003 == 	209406
 replace UCID_1934_2003 = 	227479	 if UCID_1934_2003 == 	220397
 replace UCID_1934_2003 = 	407130	 if UCID_1934_2003 == 	224328
 replace UCID_1934_2003 = 	323947	 if UCID_1934_2003 == 	226793
 replace UCID_1934_2003 = 	226454	 if UCID_1934_2003 == 	227495
 replace UCID_1934_2003 = 	236512	 if UCID_1934_2003 == 	236646
 replace UCID_1934_2003 = 	17673	 if UCID_1934_2003 == 	243702
 replace UCID_1934_2003 = 	237623	 if UCID_1934_2003 == 	251069
 replace UCID_1934_2003 = 	147449	 if UCID_1934_2003 == 	259964
 replace UCID_1934_2003 = 	75526	 if UCID_1934_2003 == 	269062
 replace UCID_1934_2003 = 	142881	 if UCID_1934_2003 == 	277091
 replace UCID_1934_2003 = 	233865	 if UCID_1934_2003 == 	289272
 replace UCID_1934_2003 = 	54989	 if UCID_1934_2003 == 	295708
 replace UCID_1934_2003 = 	167480	 if UCID_1934_2003 == 	302027
 replace UCID_1934_2003 = 	336835	 if UCID_1934_2003 == 	335053
 replace UCID_1934_2003 = 	353252	 if UCID_1934_2003 == 	337286
 replace UCID_1934_2003 = 	355893	 if UCID_1934_2003 == 	339833
 replace UCID_1934_2003 = 	344136	 if UCID_1934_2003 == 	340438
 replace UCID_1934_2003 = 	344357	 if UCID_1934_2003 == 	342345
 replace UCID_1934_2003 = 	324362	 if UCID_1934_2003 == 	350454
 replace UCID_1934_2003 = 	370521	 if UCID_1934_2003 == 	372889
 replace UCID_1934_2003 = 	344532	 if UCID_1934_2003 == 	379550
 replace UCID_1934_2003 = 	202206	 if UCID_1934_2003 == 	408984
 replace UCID_1934_2003 = 	124031	 if UCID_1934_2003 == 	414301
 replace UCID_1934_2003 = 	278221	 if UCID_1934_2003 == 	416240
 replace UCID_1934_2003 = 	164342	 if UCID_1934_2003 == 	425891

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_67.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_66.dta"

*******************************************************************************
* 1997 & 1998

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_67.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1997 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1997
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1997
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1997_remaining_2.dta", replace
	restore

*Second year 1998 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1998
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1998
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1998_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1997 & 1998
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1997_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1998_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1997_1998_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1997_1998_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1997_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1998_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1997_1998_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1997_1998_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1997_1998_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1997 to 1998
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_67.dta", clear

replace UCID_1934_2003 = 229922 if UCID_1934_2003 == 229596

 replace UCID_1934_2003 = 	79230	 if UCID_1934_2003 == 	1087
 replace UCID_1934_2003 = 	265028	 if UCID_1934_2003 == 	2111
 replace UCID_1934_2003 = 	145548	 if UCID_1934_2003 == 	3523
 replace UCID_1934_2003 = 	332043	 if UCID_1934_2003 == 	6872
 replace UCID_1934_2003 = 	17451	 if UCID_1934_2003 == 	16571
 replace UCID_1934_2003 = 	23389	 if UCID_1934_2003 == 	23383
 replace UCID_1934_2003 = 	226467	 if UCID_1934_2003 == 	25084
 replace UCID_1934_2003 = 	242929	 if UCID_1934_2003 == 	29298
 replace UCID_1934_2003 = 	32391	 if UCID_1934_2003 == 	32507
 replace UCID_1934_2003 = 	314712	 if UCID_1934_2003 == 	35950
 replace UCID_1934_2003 = 	298622	 if UCID_1934_2003 == 	38382
 replace UCID_1934_2003 = 	39289	 if UCID_1934_2003 == 	40880
 replace UCID_1934_2003 = 	51311	 if UCID_1934_2003 == 	49562
 replace UCID_1934_2003 = 	51483	 if UCID_1934_2003 == 	49710
 replace UCID_1934_2003 = 	55254	 if UCID_1934_2003 == 	55303
 replace UCID_1934_2003 = 	244359	 if UCID_1934_2003 == 	55449
 replace UCID_1934_2003 = 	141750	 if UCID_1934_2003 == 	58730
 replace UCID_1934_2003 = 	66675	 if UCID_1934_2003 == 	67181
 replace UCID_1934_2003 = 	67079	 if UCID_1934_2003 == 	67252
 replace UCID_1934_2003 = 	51035	 if UCID_1934_2003 == 	67352
 replace UCID_1934_2003 = 	71669	 if UCID_1934_2003 == 	85285
 replace UCID_1934_2003 = 	287786	 if UCID_1934_2003 == 	86185
 replace UCID_1934_2003 = 	11211	 if UCID_1934_2003 == 	89245
 replace UCID_1934_2003 = 	65796	 if UCID_1934_2003 == 	95048
 replace UCID_1934_2003 = 	167849	 if UCID_1934_2003 == 	97393
 replace UCID_1934_2003 = 	59491	 if UCID_1934_2003 == 	99119
 replace UCID_1934_2003 = 	163379	 if UCID_1934_2003 == 	100181
 replace UCID_1934_2003 = 	294617	 if UCID_1934_2003 == 	108875
 replace UCID_1934_2003 = 	118998	 if UCID_1934_2003 == 	111384
 replace UCID_1934_2003 = 	200718	 if UCID_1934_2003 == 	116513
 replace UCID_1934_2003 = 	127252	 if UCID_1934_2003 == 	117407
 replace UCID_1934_2003 = 	129719	 if UCID_1934_2003 == 	120921
 replace UCID_1934_2003 = 	130973	 if UCID_1934_2003 == 	131108
 replace UCID_1934_2003 = 	135235	 if UCID_1934_2003 == 	134907
 replace UCID_1934_2003 = 	144520	 if UCID_1934_2003 == 	136048
 replace UCID_1934_2003 = 	145407	 if UCID_1934_2003 == 	145577
 replace UCID_1934_2003 = 	62914	 if UCID_1934_2003 == 	163548
 replace UCID_1934_2003 = 	185500	 if UCID_1934_2003 == 	164451
 replace UCID_1934_2003 = 	174980	 if UCID_1934_2003 == 	166649
 replace UCID_1934_2003 = 	200660	 if UCID_1934_2003 == 	190927
 replace UCID_1934_2003 = 	191679	 if UCID_1934_2003 == 	199712
 replace UCID_1934_2003 = 	168400	 if UCID_1934_2003 == 	222133
 replace UCID_1934_2003 = 	231690	 if UCID_1934_2003 == 	233566
 replace UCID_1934_2003 = 	242554	 if UCID_1934_2003 == 	242087
 replace UCID_1934_2003 = 	265186	 if UCID_1934_2003 == 	265238
 replace UCID_1934_2003 = 	271951	 if UCID_1934_2003 == 	268083
 replace UCID_1934_2003 = 	271663	 if UCID_1934_2003 == 	274652
 replace UCID_1934_2003 = 	275402	 if UCID_1934_2003 == 	275691
 replace UCID_1934_2003 = 	282373	 if UCID_1934_2003 == 	282006
 replace UCID_1934_2003 = 	29134	 if UCID_1934_2003 == 	282541
 replace UCID_1934_2003 = 	250145	 if UCID_1934_2003 == 	287262
 replace UCID_1934_2003 = 	379494	 if UCID_1934_2003 == 	301861
 replace UCID_1934_2003 = 	62206	 if UCID_1934_2003 == 	313814
 replace UCID_1934_2003 = 	301997	 if UCID_1934_2003 == 	317540
 replace UCID_1934_2003 = 	346343	 if UCID_1934_2003 == 	323221
 replace UCID_1934_2003 = 	337686	 if UCID_1934_2003 == 	330540
 replace UCID_1934_2003 = 	352265	 if UCID_1934_2003 == 	354033
 replace UCID_1934_2003 = 	165056	 if UCID_1934_2003 == 	354955
 replace UCID_1934_2003 = 	123037	 if UCID_1934_2003 == 	354964
 replace UCID_1934_2003 = 	45039	 if UCID_1934_2003 == 	369236
 replace UCID_1934_2003 = 	378469	 if UCID_1934_2003 == 	369352
 replace UCID_1934_2003 = 	429393	 if UCID_1934_2003 == 	429363

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_68.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_67.dta"

*******************************************************************************
* 1998 & 1999

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_68.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1998 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1998
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1998
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1998_remaining_2.dta", replace
	restore

*Second year 1999 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1999
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1999
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1999_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1998 & 1999
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1998_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1999_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1998_1999_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1998_1999_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1998_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1999_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1998_1999_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1998_1999_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1998_1999_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1998 to 1999
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_68.dta", clear

 replace UCID_1934_2003 = 	1771	 if UCID_1934_2003 == 	2148
 replace UCID_1934_2003 = 	32051	 if UCID_1934_2003 == 	2300
 replace UCID_1934_2003 = 	259	 if UCID_1934_2003 == 	5897
 replace UCID_1934_2003 = 	310	 if UCID_1934_2003 == 	8292
 replace UCID_1934_2003 = 	10039	 if UCID_1934_2003 == 	11175
 replace UCID_1934_2003 = 	11688	 if UCID_1934_2003 == 	11512
 replace UCID_1934_2003 = 	11815	 if UCID_1934_2003 == 	23720
 replace UCID_1934_2003 = 	388726	 if UCID_1934_2003 == 	26959
 replace UCID_1934_2003 = 	278568	 if UCID_1934_2003 == 	27323
 replace UCID_1934_2003 = 	11029	 if UCID_1934_2003 == 	33920
 replace UCID_1934_2003 = 	144034	 if UCID_1934_2003 == 	49373
 replace UCID_1934_2003 = 	413498	 if UCID_1934_2003 == 	54349
 replace UCID_1934_2003 = 	80064	 if UCID_1934_2003 == 	58497
 replace UCID_1934_2003 = 	65271	 if UCID_1934_2003 == 	63748
 replace UCID_1934_2003 = 	67662	 if UCID_1934_2003 == 	63835
 replace UCID_1934_2003 = 	236781	 if UCID_1934_2003 == 	65861
 replace UCID_1934_2003 = 	75006	 if UCID_1934_2003 == 	68967
 replace UCID_1934_2003 = 	51306	 if UCID_1934_2003 == 	78672
 replace UCID_1934_2003 = 	169623	 if UCID_1934_2003 == 	82761
 replace UCID_1934_2003 = 	409164	 if UCID_1934_2003 == 	85876
 replace UCID_1934_2003 = 	63886	 if UCID_1934_2003 == 	89045
 replace UCID_1934_2003 = 	66419	 if UCID_1934_2003 == 	89952
 replace UCID_1934_2003 = 	165165	 if UCID_1934_2003 == 	94123
 replace UCID_1934_2003 = 	65875	 if UCID_1934_2003 == 	98256
 replace UCID_1934_2003 = 	53334	 if UCID_1934_2003 == 	100971
 replace UCID_1934_2003 = 	12002	 if UCID_1934_2003 == 	111467
 replace UCID_1934_2003 = 	112649	 if UCID_1934_2003 == 	119589
 replace UCID_1934_2003 = 	367298	 if UCID_1934_2003 == 	122363
 replace UCID_1934_2003 = 	134884	 if UCID_1934_2003 == 	125521
 replace UCID_1934_2003 = 	354325	 if UCID_1934_2003 == 	141179
 replace UCID_1934_2003 = 	300437	 if UCID_1934_2003 == 	141761
 replace UCID_1934_2003 = 	167819	 if UCID_1934_2003 == 	166452
 replace UCID_1934_2003 = 	337908	 if UCID_1934_2003 == 	175806
 replace UCID_1934_2003 = 	163430	 if UCID_1934_2003 == 	178474
 replace UCID_1934_2003 = 	144566	 if UCID_1934_2003 == 	182395
 replace UCID_1934_2003 = 	308026	 if UCID_1934_2003 == 	182720
 replace UCID_1934_2003 = 	191724	 if UCID_1934_2003 == 	194393
 replace UCID_1934_2003 = 	226272	 if UCID_1934_2003 == 	218997
 replace UCID_1934_2003 = 	223572	 if UCID_1934_2003 == 	221889
 replace UCID_1934_2003 = 	224786	 if UCID_1934_2003 == 	226237
 replace UCID_1934_2003 = 	245857	 if UCID_1934_2003 == 	247356
 replace UCID_1934_2003 = 	249207	 if UCID_1934_2003 == 	248768
 replace UCID_1934_2003 = 	267019	 if UCID_1934_2003 == 	266417
 replace UCID_1934_2003 = 	276714	 if UCID_1934_2003 == 	277023
 replace UCID_1934_2003 = 	208741	 if UCID_1934_2003 == 	280302
 replace UCID_1934_2003 = 	225478	 if UCID_1934_2003 == 	286194
 replace UCID_1934_2003 = 	208314	 if UCID_1934_2003 == 	288038
 replace UCID_1934_2003 = 	312128	 if UCID_1934_2003 == 	304061
 replace UCID_1934_2003 = 	148394	 if UCID_1934_2003 == 	314185
 replace UCID_1934_2003 = 	324198	 if UCID_1934_2003 == 	324177
 replace UCID_1934_2003 = 	339304	 if UCID_1934_2003 == 	327396
 replace UCID_1934_2003 = 	339749	 if UCID_1934_2003 == 	341119
 replace UCID_1934_2003 = 	391755	 if UCID_1934_2003 == 	344125
 replace UCID_1934_2003 = 	324608	 if UCID_1934_2003 == 	346164
 replace UCID_1934_2003 = 	193869	 if UCID_1934_2003 == 	364333
 replace UCID_1934_2003 = 	362639	 if UCID_1934_2003 == 	365174
 replace UCID_1934_2003 = 	340277	 if UCID_1934_2003 == 	366980
 replace UCID_1934_2003 = 	402822	 if UCID_1934_2003 == 	380996
 replace UCID_1934_2003 = 	88703	 if UCID_1934_2003 == 	391154
 replace UCID_1934_2003 = 	346743	 if UCID_1934_2003 == 	393423
 replace UCID_1934_2003 = 	421283	 if UCID_1934_2003 == 	422729
 replace UCID_1934_2003 = 	426797	 if UCID_1934_2003 == 	426798
 replace UCID_1934_2003 = 	405508	 if UCID_1934_2003 == 	427499

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_69.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_68.dta"

*******************************************************************************
* 1999 & 2000

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_69.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 1999 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 1999
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 1999
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1999_remaining_2.dta", replace
	restore

*Second year 2000 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 2000
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 2000
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2000_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 1999 & 2000
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1999_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2000_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1999_2000_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1999_2000_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1999_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2000_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1999_2000_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1999_2000_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_1999_2000_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 1999 to 2000
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_69.dta", clear

 replace UCID_1934_2003 = 	108348	 if UCID_1934_2003 == 	2068
 replace UCID_1934_2003 = 	236613	 if UCID_1934_2003 == 	7612
 replace UCID_1934_2003 = 	426345	 if UCID_1934_2003 == 	13643
 replace UCID_1934_2003 = 	14668	 if UCID_1934_2003 == 	14399
 replace UCID_1934_2003 = 	2182	 if UCID_1934_2003 == 	14633
 replace UCID_1934_2003 = 	17188	 if UCID_1934_2003 == 	16837
 replace UCID_1934_2003 = 	16544	 if UCID_1934_2003 == 	17408
 replace UCID_1934_2003 = 	103757	 if UCID_1934_2003 == 	22243
 replace UCID_1934_2003 = 	122791	 if UCID_1934_2003 == 	23320
 replace UCID_1934_2003 = 	27137	 if UCID_1934_2003 == 	25876
 replace UCID_1934_2003 = 	14512	 if UCID_1934_2003 == 	27521
 replace UCID_1934_2003 = 	247034	 if UCID_1934_2003 == 	29488
 replace UCID_1934_2003 = 	31771	 if UCID_1934_2003 == 	30861
 replace UCID_1934_2003 = 	9252	 if UCID_1934_2003 == 	30874
 replace UCID_1934_2003 = 	34615	 if UCID_1934_2003 == 	34559
 replace UCID_1934_2003 = 	264931	 if UCID_1934_2003 == 	39189
 replace UCID_1934_2003 = 	48343	 if UCID_1934_2003 == 	48664
 replace UCID_1934_2003 = 	57011	 if UCID_1934_2003 == 	55861
 replace UCID_1934_2003 = 	73385	 if UCID_1934_2003 == 	56751
 replace UCID_1934_2003 = 	59574	 if UCID_1934_2003 == 	58747
 replace UCID_1934_2003 = 	425512	 if UCID_1934_2003 == 	60862
 replace UCID_1934_2003 = 	89866	 if UCID_1934_2003 == 	64324
 replace UCID_1934_2003 = 	85545	 if UCID_1934_2003 == 	65209
 replace UCID_1934_2003 = 	75822	 if UCID_1934_2003 == 	67349
 replace UCID_1934_2003 = 	70877	 if UCID_1934_2003 == 	73492
 replace UCID_1934_2003 = 	128215	 if UCID_1934_2003 == 	110839
 replace UCID_1934_2003 = 	392699	 if UCID_1934_2003 == 	112228
 replace UCID_1934_2003 = 	153695	 if UCID_1934_2003 == 	113985
 replace UCID_1934_2003 = 	117665	 if UCID_1934_2003 == 	117333
 replace UCID_1934_2003 = 	127083	 if UCID_1934_2003 == 	133781
 replace UCID_1934_2003 = 	107721	 if UCID_1934_2003 == 	136705
 replace UCID_1934_2003 = 	266445	 if UCID_1934_2003 == 	137070
 replace UCID_1934_2003 = 	145293	 if UCID_1934_2003 == 	145700
 replace UCID_1934_2003 = 	79842	 if UCID_1934_2003 == 	170232
 replace UCID_1934_2003 = 	299214	 if UCID_1934_2003 == 	182141
 replace UCID_1934_2003 = 	312401	 if UCID_1934_2003 == 	187156
 replace UCID_1934_2003 = 	163341	 if UCID_1934_2003 == 	194363
 replace UCID_1934_2003 = 	287013	 if UCID_1934_2003 == 	202216
 replace UCID_1934_2003 = 	34887	 if UCID_1934_2003 == 	211941
 replace UCID_1934_2003 = 	432893	 if UCID_1934_2003 == 	213316
 replace UCID_1934_2003 = 	224001	 if UCID_1934_2003 == 	213819
 replace UCID_1934_2003 = 	227414	 if UCID_1934_2003 == 	216366
 replace UCID_1934_2003 = 	225489	 if UCID_1934_2003 == 	218431
 replace UCID_1934_2003 = 	225591	 if UCID_1934_2003 == 	222783
 replace UCID_1934_2003 = 	225668	 if UCID_1934_2003 == 	226251
 replace UCID_1934_2003 = 	26572	 if UCID_1934_2003 == 	227584
 replace UCID_1934_2003 = 	211212	 if UCID_1934_2003 == 	228563
 replace UCID_1934_2003 = 	243648	 if UCID_1934_2003 == 	244028
 replace UCID_1934_2003 = 	188552	 if UCID_1934_2003 == 	266019
 replace UCID_1934_2003 = 	267111	 if UCID_1934_2003 == 	266966
 replace UCID_1934_2003 = 	168653	 if UCID_1934_2003 == 	278050
 replace UCID_1934_2003 = 	277981	 if UCID_1934_2003 == 	278357
 replace UCID_1934_2003 = 	75562	 if UCID_1934_2003 == 	285975
 replace UCID_1934_2003 = 	304622	 if UCID_1934_2003 == 	294510
 replace UCID_1934_2003 = 	410835	 if UCID_1934_2003 == 	296649
 replace UCID_1934_2003 = 	74753	 if UCID_1934_2003 == 	297385
 replace UCID_1934_2003 = 	143914	 if UCID_1934_2003 == 	304684
 replace UCID_1934_2003 = 	254796	 if UCID_1934_2003 == 	306229
 replace UCID_1934_2003 = 	320189	 if UCID_1934_2003 == 	307944
 replace UCID_1934_2003 = 	291014	 if UCID_1934_2003 == 	308667
 replace UCID_1934_2003 = 	161040	 if UCID_1934_2003 == 	311045
 replace UCID_1934_2003 = 	320260	 if UCID_1934_2003 == 	313903
 replace UCID_1934_2003 = 	320067	 if UCID_1934_2003 == 	317726
 replace UCID_1934_2003 = 	233452	 if UCID_1934_2003 == 	324715
 replace UCID_1934_2003 = 	368881	 if UCID_1934_2003 == 	329030
 replace UCID_1934_2003 = 	337654	 if UCID_1934_2003 == 	331898
 replace UCID_1934_2003 = 	388707	 if UCID_1934_2003 == 	342449
 replace UCID_1934_2003 = 	326079	 if UCID_1934_2003 == 	348257
 replace UCID_1934_2003 = 	358454	 if UCID_1934_2003 == 	348318
 replace UCID_1934_2003 = 	347388	 if UCID_1934_2003 == 	349233
 replace UCID_1934_2003 = 	177088	 if UCID_1934_2003 == 	353279
 replace UCID_1934_2003 = 	367679	 if UCID_1934_2003 == 	362005
 replace UCID_1934_2003 = 	78967	 if UCID_1934_2003 == 	365369
 replace UCID_1934_2003 = 	370254	 if UCID_1934_2003 == 	371891
 replace UCID_1934_2003 = 	368726	 if UCID_1934_2003 == 	373337
 replace UCID_1934_2003 = 	402576	 if UCID_1934_2003 == 	391418
 replace UCID_1934_2003 = 	368902	 if UCID_1934_2003 == 	392588
 replace UCID_1934_2003 = 	384058	 if UCID_1934_2003 == 	403009
 replace UCID_1934_2003 = 	378568	 if UCID_1934_2003 == 	404471
 replace UCID_1934_2003 = 	417053	 if UCID_1934_2003 == 	417956
 replace UCID_1934_2003 = 	12078	 if UCID_1934_2003 == 	426228
 replace UCID_1934_2003 = 	90916	 if UCID_1934_2003 == 	427915
 replace UCID_1934_2003 = 	50518	 if UCID_1934_2003 == 	576891

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_70.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_69.dta"

*******************************************************************************
* 2000 & 2001

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_70.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 2000 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 2000
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 2000
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2000_remaining_2.dta", replace
	restore

*Second year 2001 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 2001
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 2001
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2001_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 2000 & 2001
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2000_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2001_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2000_2001_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2000_2001_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2000_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2001_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2000_2001_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2000_2001_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_2000_2001_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 2000 to 2001
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_70.dta", clear

replace UCID_1934_2003 = 323752 if UCID_1934_2003 == 323747

 replace UCID_1934_2003 = 	11395	 if UCID_1934_2003 == 	9382
 replace UCID_1934_2003 = 	193967	 if UCID_1934_2003 == 	9928
 replace UCID_1934_2003 = 	10965	 if UCID_1934_2003 == 	12122
 replace UCID_1934_2003 = 	12421	 if UCID_1934_2003 == 	12584
 replace UCID_1934_2003 = 	158138	 if UCID_1934_2003 == 	16067
 replace UCID_1934_2003 = 	403534	 if UCID_1934_2003 == 	16254
 replace UCID_1934_2003 = 	236689	 if UCID_1934_2003 == 	21104
 replace UCID_1934_2003 = 	231184	 if UCID_1934_2003 == 	27692
 replace UCID_1934_2003 = 	27657	 if UCID_1934_2003 == 	27695
 replace UCID_1934_2003 = 	30821	 if UCID_1934_2003 == 	31819
 replace UCID_1934_2003 = 	161	 if UCID_1934_2003 == 	41254
 replace UCID_1934_2003 = 	48138	 if UCID_1934_2003 == 	48643
 replace UCID_1934_2003 = 	66735	 if UCID_1934_2003 == 	50067
 replace UCID_1934_2003 = 	250717	 if UCID_1934_2003 == 	50605
 replace UCID_1934_2003 = 	106161	 if UCID_1934_2003 == 	52865
 replace UCID_1934_2003 = 	233393	 if UCID_1934_2003 == 	53182
 replace UCID_1934_2003 = 	96142	 if UCID_1934_2003 == 	62336
 replace UCID_1934_2003 = 	1670	 if UCID_1934_2003 == 	64124
 replace UCID_1934_2003 = 	70450	 if UCID_1934_2003 == 	68176
 replace UCID_1934_2003 = 	257543	 if UCID_1934_2003 == 	69076
 replace UCID_1934_2003 = 	349489	 if UCID_1934_2003 == 	71229
 replace UCID_1934_2003 = 	60925	 if UCID_1934_2003 == 	71277
 replace UCID_1934_2003 = 	72709	 if UCID_1934_2003 == 	74051
 replace UCID_1934_2003 = 	311885	 if UCID_1934_2003 == 	79538
 replace UCID_1934_2003 = 	60898	 if UCID_1934_2003 == 	82851
 replace UCID_1934_2003 = 	248943	 if UCID_1934_2003 == 	85528
 replace UCID_1934_2003 = 	306018	 if UCID_1934_2003 == 	86459
 replace UCID_1934_2003 = 	55728	 if UCID_1934_2003 == 	97274
 replace UCID_1934_2003 = 	54632	 if UCID_1934_2003 == 	103050
 replace UCID_1934_2003 = 	79272	 if UCID_1934_2003 == 	107736
 replace UCID_1934_2003 = 	11525	 if UCID_1934_2003 == 	112247
 replace UCID_1934_2003 = 	123819	 if UCID_1934_2003 == 	112742
 replace UCID_1934_2003 = 	118045	 if UCID_1934_2003 == 	113055
 replace UCID_1934_2003 = 	116784	 if UCID_1934_2003 == 	118409
 replace UCID_1934_2003 = 	234907	 if UCID_1934_2003 == 	118561
 replace UCID_1934_2003 = 	194022	 if UCID_1934_2003 == 	124908
 replace UCID_1934_2003 = 	67991	 if UCID_1934_2003 == 	145616
 replace UCID_1934_2003 = 	144893	 if UCID_1934_2003 == 	145831
 replace UCID_1934_2003 = 	292196	 if UCID_1934_2003 == 	155494
 replace UCID_1934_2003 = 	155409	 if UCID_1934_2003 == 	156794
 replace UCID_1934_2003 = 	329155	 if UCID_1934_2003 == 	165226
 replace UCID_1934_2003 = 	47016	 if UCID_1934_2003 == 	181914
 replace UCID_1934_2003 = 	164816	 if UCID_1934_2003 == 	182610
 replace UCID_1934_2003 = 	105109	 if UCID_1934_2003 == 	182660
 replace UCID_1934_2003 = 	293498	 if UCID_1934_2003 == 	186733
 replace UCID_1934_2003 = 	286482	 if UCID_1934_2003 == 	188880
 replace UCID_1934_2003 = 	110464	 if UCID_1934_2003 == 	199013
 replace UCID_1934_2003 = 	199814	 if UCID_1934_2003 == 	200067
 replace UCID_1934_2003 = 	321179	 if UCID_1934_2003 == 	204110
 replace UCID_1934_2003 = 	148485	 if UCID_1934_2003 == 	219031
 replace UCID_1934_2003 = 	226984	 if UCID_1934_2003 == 	221060
 replace UCID_1934_2003 = 	167995	 if UCID_1934_2003 == 	228703
 replace UCID_1934_2003 = 	245227	 if UCID_1934_2003 == 	243754
 replace UCID_1934_2003 = 	237006	 if UCID_1934_2003 == 	244289
 replace UCID_1934_2003 = 	155407	 if UCID_1934_2003 == 	244311
 replace UCID_1934_2003 = 	58536	 if UCID_1934_2003 == 	248422
 replace UCID_1934_2003 = 	265407	 if UCID_1934_2003 == 	253717
 replace UCID_1934_2003 = 	253312	 if UCID_1934_2003 == 	255659
 replace UCID_1934_2003 = 	255812	 if UCID_1934_2003 == 	256127
 replace UCID_1934_2003 = 	241379	 if UCID_1934_2003 == 	260050
 replace UCID_1934_2003 = 	254271	 if UCID_1934_2003 == 	264238
 replace UCID_1934_2003 = 	276398	 if UCID_1934_2003 == 	267210
 replace UCID_1934_2003 = 	51296	 if UCID_1934_2003 == 	276535
 replace UCID_1934_2003 = 	229921	 if UCID_1934_2003 == 	277756
 replace UCID_1934_2003 = 	309880	 if UCID_1934_2003 == 	291263
 replace UCID_1934_2003 = 	60173	 if UCID_1934_2003 == 	293025
 replace UCID_1934_2003 = 	415521	 if UCID_1934_2003 == 	302077
 replace UCID_1934_2003 = 	320199	 if UCID_1934_2003 == 	313920
 replace UCID_1934_2003 = 	320236	 if UCID_1934_2003 == 	319705
 replace UCID_1934_2003 = 	81033	 if UCID_1934_2003 == 	321325
 replace UCID_1934_2003 = 	416329	 if UCID_1934_2003 == 	323752
 replace UCID_1934_2003 = 	327083	 if UCID_1934_2003 == 	334016
 replace UCID_1934_2003 = 	348403	 if UCID_1934_2003 == 	339242
 replace UCID_1934_2003 = 	325574	 if UCID_1934_2003 == 	342222
 replace UCID_1934_2003 = 	314085	 if UCID_1934_2003 == 	342848
 replace UCID_1934_2003 = 	354635	 if UCID_1934_2003 == 	355467
 replace UCID_1934_2003 = 	344228	 if UCID_1934_2003 == 	355470
 replace UCID_1934_2003 = 	350850	 if UCID_1934_2003 == 	360774
 replace UCID_1934_2003 = 	14503	 if UCID_1934_2003 == 	361624
 replace UCID_1934_2003 = 	342892	 if UCID_1934_2003 == 	386760
 replace UCID_1934_2003 = 	21614	 if UCID_1934_2003 == 	392747
 replace UCID_1934_2003 = 	94779	 if UCID_1934_2003 == 	396466
 replace UCID_1934_2003 = 	249871	 if UCID_1934_2003 == 	412331
 replace UCID_1934_2003 = 	11016	 if UCID_1934_2003 == 	425866
 replace UCID_1934_2003 = 	20802	 if UCID_1934_2003 == 	426342

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_71.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_70.dta"

*******************************************************************************
* 2001 & 2002

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_71.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 2001 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 2001
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 2001
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2001_remaining_2.dta", replace
	restore

*Second year 2002 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 2002
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 2002
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2002_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 2001 & 2002
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2001_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2002_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2001_2002_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2001_2002_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2001_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2002_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2001_2002_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2001_2002_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_2001_2002_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 2001 to 2002
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_71.dta", clear

 replace UCID_1934_2003 = 	615	 if UCID_1934_2003 == 	76
 replace UCID_1934_2003 = 	3341	 if UCID_1934_2003 == 	2313
 replace UCID_1934_2003 = 	1251	 if UCID_1934_2003 == 	3714
 replace UCID_1934_2003 = 	48193	 if UCID_1934_2003 == 	6592
 replace UCID_1934_2003 = 	345745	 if UCID_1934_2003 == 	8899
 replace UCID_1934_2003 = 	991	 if UCID_1934_2003 == 	9252
 replace UCID_1934_2003 = 	12166	 if UCID_1934_2003 == 	11033
 replace UCID_1934_2003 = 	17039	 if UCID_1934_2003 == 	16756
 replace UCID_1934_2003 = 	200115	 if UCID_1934_2003 == 	22681
 replace UCID_1934_2003 = 	9603	 if UCID_1934_2003 == 	23439
 replace UCID_1934_2003 = 	325787	 if UCID_1934_2003 == 	23539
 replace UCID_1934_2003 = 	26822	 if UCID_1934_2003 == 	26972
 replace UCID_1934_2003 = 	254194	 if UCID_1934_2003 == 	31316
 replace UCID_1934_2003 = 	291469	 if UCID_1934_2003 == 	35637
 replace UCID_1934_2003 = 	38416	 if UCID_1934_2003 == 	37722
 replace UCID_1934_2003 = 	22104	 if UCID_1934_2003 == 	38672
 replace UCID_1934_2003 = 	64492	 if UCID_1934_2003 == 	52056
 replace UCID_1934_2003 = 	271412	 if UCID_1934_2003 == 	58088
 replace UCID_1934_2003 = 	104504	 if UCID_1934_2003 == 	60412
 replace UCID_1934_2003 = 	227552	 if UCID_1934_2003 == 	67449
 replace UCID_1934_2003 = 	142242	 if UCID_1934_2003 == 	70669
 replace UCID_1934_2003 = 	65360	 if UCID_1934_2003 == 	74803
 replace UCID_1934_2003 = 	403054	 if UCID_1934_2003 == 	76532
 replace UCID_1934_2003 = 	65322	 if UCID_1934_2003 == 	88736
 replace UCID_1934_2003 = 	58967	 if UCID_1934_2003 == 	93927
 replace UCID_1934_2003 = 	303824	 if UCID_1934_2003 == 	94911
 replace UCID_1934_2003 = 	53644	 if UCID_1934_2003 == 	98927
 replace UCID_1934_2003 = 	205360	 if UCID_1934_2003 == 	133122
 replace UCID_1934_2003 = 	130516	 if UCID_1934_2003 == 	133331
 replace UCID_1934_2003 = 	143999	 if UCID_1934_2003 == 	136128
 replace UCID_1934_2003 = 	144018	 if UCID_1934_2003 == 	137750
 replace UCID_1934_2003 = 	62959	 if UCID_1934_2003 == 	142812
 replace UCID_1934_2003 = 	55032	 if UCID_1934_2003 == 	147905
 replace UCID_1934_2003 = 	402027	 if UCID_1934_2003 == 	149072
 replace UCID_1934_2003 = 	148224	 if UCID_1934_2003 == 	153060
 replace UCID_1934_2003 = 	255151	 if UCID_1934_2003 == 	164454
 replace UCID_1934_2003 = 	188918	 if UCID_1934_2003 == 	167772
 replace UCID_1934_2003 = 	393672	 if UCID_1934_2003 == 	174854
 replace UCID_1934_2003 = 	294359	 if UCID_1934_2003 == 	180698
 replace UCID_1934_2003 = 	148627	 if UCID_1934_2003 == 	181839
 replace UCID_1934_2003 = 	164714	 if UCID_1934_2003 == 	185114
 replace UCID_1934_2003 = 	349032	 if UCID_1934_2003 == 	203318
 replace UCID_1934_2003 = 	227779	 if UCID_1934_2003 == 	218314
 replace UCID_1934_2003 = 	279259	 if UCID_1934_2003 == 	223591
 replace UCID_1934_2003 = 	212460	 if UCID_1934_2003 == 	226446
 replace UCID_1934_2003 = 	222849	 if UCID_1934_2003 == 	227288
 replace UCID_1934_2003 = 	230062	 if UCID_1934_2003 == 	228770
 replace UCID_1934_2003 = 	38539	 if UCID_1934_2003 == 	241523
 replace UCID_1934_2003 = 	254928	 if UCID_1934_2003 == 	254726
 replace UCID_1934_2003 = 	265575	 if UCID_1934_2003 == 	267021
 replace UCID_1934_2003 = 	271428	 if UCID_1934_2003 == 	270208
 replace UCID_1934_2003 = 	269219	 if UCID_1934_2003 == 	270474
 replace UCID_1934_2003 = 	198768	 if UCID_1934_2003 == 	271120
 replace UCID_1934_2003 = 	269838	 if UCID_1934_2003 == 	272744
 replace UCID_1934_2003 = 	106250	 if UCID_1934_2003 == 	274067
 replace UCID_1934_2003 = 	150343	 if UCID_1934_2003 == 	274549
 replace UCID_1934_2003 = 	274913	 if UCID_1934_2003 == 	275774
 replace UCID_1934_2003 = 	314718	 if UCID_1934_2003 == 	291260
 replace UCID_1934_2003 = 	317207	 if UCID_1934_2003 == 	309861
 replace UCID_1934_2003 = 	312818	 if UCID_1934_2003 == 	311545
 replace UCID_1934_2003 = 	37724	 if UCID_1934_2003 == 	313512
 replace UCID_1934_2003 = 	359232	 if UCID_1934_2003 == 	322861
 replace UCID_1934_2003 = 	323637	 if UCID_1934_2003 == 	339026
 replace UCID_1934_2003 = 	267897	 if UCID_1934_2003 == 	339254
 replace UCID_1934_2003 = 	321479	 if UCID_1934_2003 == 	356840
 replace UCID_1934_2003 = 	365592	 if UCID_1934_2003 == 	365833
 replace UCID_1934_2003 = 	342937	 if UCID_1934_2003 == 	368277
 replace UCID_1934_2003 = 	110778	 if UCID_1934_2003 == 	404756
 replace UCID_1934_2003 = 	250001	 if UCID_1934_2003 == 	415460
 replace UCID_1934_2003 = 	17625	 if UCID_1934_2003 == 	425655
 replace UCID_1934_2003 = 	234275	 if UCID_1934_2003 == 	432446

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_72.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_71.dta"

*******************************************************************************
* 2002 & 2003

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_72.dta", clear

* flag UCID_1934_2003 groups with gaps
preserve 
xtset UCID_1934_2003 periode
tsfill

	gen gap2 = 1 if missing(year)
	by UCID_1934_2003 (periode), sort: egen group_gap2 = max(gap2)
	by UCID_1934_2003 (periode), sort: gen group_gap2_full = 1 if group_gap2 == 1
	keep if group_gap2 == 1

	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)
	bysort UCID_1934_2003: replace year = year[_n-1] if missing(year)

	sort UCID_1934_2003 periode
	bysort UCID_1934_2003 year:  gen dup = cond(_N==1,0,_n)

	bysort UCID_1934_2003: replace CID = CID[_n+1] if missing(CID)
	bysort UCID_1934_2003 CID:  gen dup2 = cond(_N==1,0,_n)

	drop if dup == 0 & dup2 == 0
	sort UCID_1934_2003 periode
	drop if missing(cname_temp)
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps", replace
restore

*Subsample of single obs & start and end of series
sort UCID_1934_2003 year
gen tag = 1
egen obs_total = total(tag), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxpanel = max(dup_panel)
bysort UCID_1934_2003: egen minpanel = min(dup_panel)

preserve
keep if dup_panel == maxpanel | dup_panel == minpanel
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", replace
restore

* Subsample of potential matches
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_potential_match.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
drop tag obs_total group_gap2 gap2 dup dup2
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", replace	
	
*Split in crossections

/*
1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
*/

*First year 2002 : dropping start of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 2002
		drop if dup_panel == minpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 2002
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2002_remaining_2.dta", replace
	restore

*Second year 2003 : dropping end of series
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_remaining_2_obs.dta", clear
sort UCID_1934_2003 year

	preserve
		keep if year == 2003
		drop if dup_panel == maxpanel
		
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_obs_gaps"
	drop dup 
	
	bysort CID year:  gen dup = cond(_N==1,0,_n)
	drop if dup > 1
	keep if year == 2003
		
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2003_remaining_2.dta", replace
	restore


*Manual evaluation for each pair of years (no loop, to reduce amount of match)
*Matching 2002 & 2003
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2002_remaining_2.dta", clear
matchit UCID_1934_2003 cname_temp using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2003_remaining_2.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.80
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2002_2003_matched_remaining_2.dta", replace
	
* Add gdenr
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2002_2003_matched_remaining_2.dta", clear
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2002_remaining_2.dta"
order UCID_1934_2003 cname_temp year gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 _merge_1 fullname_1)
rename (UCID_1934_20031 cname_temp1)(UCID_1934_2003 cname_temp)
drop if _merge_1 == 2
merge m:1 UCID_1934_2003 using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2003_remaining_2.dta"
order UCID_1934_2003_1 cname_temp_1 year_1 fullname_1 gdenr_2018_1 gdename_1 UCID_1934_2003 cname_temp year fullname gdenr_2018 gdename
rename (UCID_1934_2003 cname_temp year gdenr_2018 gdename _merge fullname) (UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 _merge_2 fullname_2)
drop CID _merge_long_geo_comp max_id2 cname_orig gdenr gdenr_2012 gdename_orig cname_group owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract ctn E_CNTR N_CNTR nb_owners cname
keep if gdenr_2018_1 != gdenr_2018_2
drop if UCID_1934_2003_1 == UCID_1934_2003_2
sort UCID_1934_2003_1 similscore

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2002_2003_matched_remaining_2_toeval.dta", replace

*export to xlsx
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_2002_2003_matched_remaining_2_toeval.dta", clear
drop if missing(year_1)
sort similscore
gen correct = .
order UCID_1934_2003_1 cname_temp_1 year_1 gdenr_2018_1 gdename_1 correct UCID_1934_2003_2 cname_temp_2 year_2 gdenr_2018_2 gdename_2 fullname_1 fullname_2
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Yannick_Eval\firm_panel_2002_2003_matched_remaining_2_toeval.xlsx", replace firstrow(variables)

* Manual corrections based for 2002 to 2003
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_72.dta", clear

 replace UCID_1934_2003 = 	22726	 if UCID_1934_2003 == 	111051
 replace UCID_1934_2003 = 	340464	 if UCID_1934_2003 == 	338384
 replace UCID_1934_2003 = 	386979	 if UCID_1934_2003 == 	385989
 replace UCID_1934_2003 = 	50518	 if UCID_1934_2003 == 	473295

* Deal with new duplicates within one year
drop dup4
bysort UCID_1934_2003 year: gen dup4 = cond(_N==1,0,_n)
	*keep if dup4 > 0


* Deduplicate
foreach y in owners function signature page CID_str cname_group {
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_n-1] + "; " + `y' if dup4>0
	replace `y' = " "+`y' if dup4>0
	replace `y' = usubinstr(`y', " ; ", "",.) if dup4>0
	bysort UCID_1934_2003 year (max_id2): replace `y' = `y'[_N] if dup4>0
	}
	
	bysort UCID_1934_2003: replace gdenr = gdenr[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr = gdenr[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr) & dup4>0
	bysort UCID_1934_2003: replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr) & dup4>0

drop if dup4 > 1 
 
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_73.dta", replace
erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_72.dta"

********************************************************************************

********************************************************************************
***** Finale modifications

*** Gen New UCID_1934_2003 going from 1 to N

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_73.dta", clear

egen ID=group(UCID_1934_2003)
replace UCID_1934_2003 = ID
drop ID

***Missing gdenr, gdenr_2012, gdenr_2018

bysort UCID_1934_2003 (year): replace gdenr = gdenr[_n-1] if missing(gdenr)
bysort UCID_1934_2003 (year): replace gdenr_2012 = gdenr_2012[_n-1] if missing(gdenr_2012)
bysort UCID_1934_2003 (year): replace gdenr_2018 = gdenr_2018[_n-1] if missing(gdenr_2018)

forval i = 1/19 {
bysort UCID_1934_2003 (year): replace gdenr = gdenr[_n+1] if missing(gdenr)
bysort UCID_1934_2003 (year): replace gdenr_2012 = gdenr_2012[_n+1] if missing(gdenr_2012)
bysort UCID_1934_2003 (year): replace gdenr_2018 = gdenr_2018[_n+1] if missing(gdenr_2018)
}

***Corrections with Pully, Prilly and Fully

/*
preserve
keep if gdenr == 5589 | gdenr == 5590 | gdenr == 6133 | gdenr == 5855 | gdenr == 5530
bysort UCID_1934_2003 (year):  gen occurences = cond(_N==1,0,_n) //nb of appearance by ID, creation when = 1
bysort UCID_1934_2003 (year):  egen max_occurences = max(occurences) //years of appearance, death when nb of appearance = max appearance

* Creation of firms
by UCID_1934_2003 (year), sort: gen creation_2012 = 1 if occurences == 1 & year != 1934
by UCID_1934_2003 (year), sort: gen creation_2018 = 1 if occurences == 1 & year != 1934
by UCID_1934_2003 (year), sort: replace creation_2012 = 1 if occurences == 0 & year != 1934
by UCID_1934_2003 (year), sort: replace creation_2018 = 1 if occurences == 0 & year != 1934

* Deaths (or disappearing from our data)
by UCID_1934_2003 (year), sort: gen death_2012 = 1 if occurences == max_occurences & year != 2003
by UCID_1934_2003 (year), sort: gen death_2018 = 1 if occurences == max_occurences & year != 2003

* Relocation of firms
by UCID_1934_2003 (year), sort: gen relocation_in_2012 = 1 if gdenr_2012 != gdenr_2012[_n-1] & creation_2012 != 1
by UCID_1934_2003 (year), sort: gen relocation_in_2018 = 1 if gdenr_2018 != gdenr_2018[_n-1] & creation_2018 != 1
replace relocation_in_2012 = 0 if year == 1934
replace relocation_in_2018 = 0 if year == 1934

by UCID_1934_2003 (year), sort: gen relocation_out_2012 = 1 if gdenr_2012 != gdenr_2012[_n+1] & max_occurences != occurences
by UCID_1934_2003 (year), sort: gen relocation_out_2018 = 1 if gdenr_2018 != gdenr_2018[_n+1] & max_occurences != occurences

*export to eval
order CID UCID_1934_2003 year gdenr gdenr_2012 gdenr_2018 relocation*
*export

restore
*/

*Actual corrections

replace gdenr = 5590 if year == 1960 & UCID_1934_2003 == 33397
replace gdenr_2012 = 5590 if year == 1960 & UCID_1934_2003 == 33397
replace gdenr_2018 = 5590 if year == 1960 & UCID_1934_2003 == 33397

replace gdenr = 5590 if year == 1982 & UCID_1934_2003 == 35714
replace gdenr_2012 = 5590 if year == 1982 & UCID_1934_2003 == 35714
replace gdenr_2018 = 5590 if year == 1982 & UCID_1934_2003 == 35714

replace gdenr = 5590 if year == 1991 & UCID_1934_2003 == 35998
replace gdenr_2012 = 5590 if year == 1991 & UCID_1934_2003 == 35998
replace gdenr_2018 = 5590 if year == 1991 & UCID_1934_2003 == 35998

replace gdenr = 5590 if year == 1997 & UCID_1934_2003 == 36100
replace gdenr_2012 = 5590 if year == 1997 & UCID_1934_2003 == 36100
replace gdenr_2018 = 5590 if year == 1997 & UCID_1934_2003 == 36100

replace gdenr = 5590 if year == 1986 & UCID_1934_2003 == 36312
replace gdenr_2012 = 5590 if year == 1986 & UCID_1934_2003 == 36312
replace gdenr_2018 = 5590 if year == 1986 & UCID_1934_2003 == 36312

replace gdenr = 5589 if year == 1988 & UCID_1934_2003 == 36356
replace gdenr_2012 = 5589 if year == 1988 & UCID_1934_2003 == 36356
replace gdenr_2018 = 5589 if year == 1988 & UCID_1934_2003 == 36356

replace gdenr = 5590 if year == 1987 & UCID_1934_2003 == 36430
replace gdenr_2012 = 5590 if year == 1987 & UCID_1934_2003 == 36430
replace gdenr_2018 = 5590 if year == 1987 & UCID_1934_2003 == 36430

replace gdenr = 5590 if year == 1991 & UCID_1934_2003 == 36526
replace gdenr_2012 = 5590 if year == 1991 & UCID_1934_2003 == 36526
replace gdenr_2018 = 5590 if year == 1991 & UCID_1934_2003 == 36526

replace gdenr = 5590 if year == 1998 & UCID_1934_2003 == 36526
replace gdenr_2012 = 5590 if year == 1998 & UCID_1934_2003 == 36526
replace gdenr_2018 = 5590 if year == 1998 & UCID_1934_2003 == 36526

replace gdenr = 5589 if year == 1989 & UCID_1934_2003 == 36590
replace gdenr_2012 = 5589 if year == 1989 & UCID_1934_2003 == 36590
replace gdenr_2018 = 5589 if year == 1989 & UCID_1934_2003 == 36590

replace gdenr = 5590 if year == 1991 & UCID_1934_2003 == 36609
replace gdenr_2012 = 5590 if year == 1991 & UCID_1934_2003 == 36609
replace gdenr_2018 = 5590 if year == 1991 & UCID_1934_2003 == 36609

replace gdenr = 5590 if year == 1991 & UCID_1934_2003 == 36610
replace gdenr_2012 = 5590 if year == 1991 & UCID_1934_2003 == 36610
replace gdenr_2018 = 5590 if year == 1991 & UCID_1934_2003 == 36610

replace gdenr = 5590 if year == 1991 & UCID_1934_2003 == 36686
replace gdenr_2012 = 5590 if year == 1991 & UCID_1934_2003 == 36686
replace gdenr_2018 = 5590 if year == 1991 & UCID_1934_2003 == 36686

replace gdenr = 5590 if year == 1995 & UCID_1934_2003 == 36883
replace gdenr_2012 = 5590 if year == 1995 & UCID_1934_2003 == 36883
replace gdenr_2018 = 5590 if year == 1995 & UCID_1934_2003 == 36883

replace gdenr = 5590 if year == 1997 & UCID_1934_2003 == 36883
replace gdenr_2012 = 5590 if year == 1997 & UCID_1934_2003 == 36883
replace gdenr_2018 = 5590 if year == 1997 & UCID_1934_2003 == 36883

replace gdenr = 5590 if year == 1991 & UCID_1934_2003 == 36933
replace gdenr_2012 = 5590 if year == 1991 & UCID_1934_2003 == 36933
replace gdenr_2018 = 5590 if year == 1991 & UCID_1934_2003 == 36933

replace gdenr = 5590 if year == 1990 & UCID_1934_2003 == 37065
replace gdenr_2012 = 5590 if year == 1990 & UCID_1934_2003 == 37065
replace gdenr_2018 = 5590 if year == 1990 & UCID_1934_2003 == 37065

replace gdenr = 5590 if year == 1991 & UCID_1934_2003 == 37131
replace gdenr_2012 = 5590 if year == 1991 & UCID_1934_2003 == 37131
replace gdenr_2018 = 5590 if year == 1991 & UCID_1934_2003 == 37131

replace gdenr = 5590 if year == 1991 & UCID_1934_2003 == 37461
replace gdenr_2012 = 5590 if year == 1991 & UCID_1934_2003 == 37461
replace gdenr_2018 = 5590 if year == 1991 & UCID_1934_2003 == 37461

replace gdenr = 5590 if year == 1996 & UCID_1934_2003 == 37464
replace gdenr_2012 = 5590 if year == 1996 & UCID_1934_2003 == 37464
replace gdenr_2018 = 5590 if year == 1996 & UCID_1934_2003 == 37464

replace gdenr = 5590 if year == 1983 & UCID_1934_2003 == 37465
replace gdenr_2012 = 5590 if year == 1983 & UCID_1934_2003 == 37465
replace gdenr_2018 = 5590 if year == 1983 & UCID_1934_2003 == 37465

replace gdenr = 5590 if year == 1964 & UCID_1934_2003 == 37492
replace gdenr_2012 = 5590 if year == 1964 & UCID_1934_2003 == 37492
replace gdenr_2018 = 5590 if year == 1964 & UCID_1934_2003 == 37492

replace gdenr = 5590 if year == 1965 & UCID_1934_2003 == 37492
replace gdenr_2012 = 5590 if year == 1965 & UCID_1934_2003 == 37492
replace gdenr_2018 = 5590 if year == 1965 & UCID_1934_2003 == 37492

replace gdenr = 5590 if year == 1969 & UCID_1934_2003 == 37501
replace gdenr_2012 = 5590 if year == 1969 & UCID_1934_2003 == 37501
replace gdenr_2018 = 5590 if year == 1969 & UCID_1934_2003 == 37501

replace gdenr = 5590 if year == 1987 & UCID_1934_2003 == 37503
replace gdenr_2012 = 5590 if year == 1987 & UCID_1934_2003 == 37503
replace gdenr_2018 = 5590 if year == 1987 & UCID_1934_2003 == 37503

replace gdenr = 5590 if year == 1990 & UCID_1934_2003 == 37504
replace gdenr_2012 = 5590 if year == 1990 & UCID_1934_2003 == 37504
replace gdenr_2018 = 5590 if year == 1990 & UCID_1934_2003 == 37504

replace gdenr = 5590 if year == 1995 & UCID_1934_2003 == 37511
replace gdenr_2012 = 5590 if year == 1995 & UCID_1934_2003 == 37511
replace gdenr_2018 = 5590 if year == 1995 & UCID_1934_2003 == 37511

replace gdenr = 5590 if year == 1984 & UCID_1934_2003 == 37513
replace gdenr_2012 = 5590 if year == 1984 & UCID_1934_2003 == 37513
replace gdenr_2018 = 5590 if year == 1984 & UCID_1934_2003 == 37513

replace gdenr = 5590 if year == 1997 & UCID_1934_2003 == 37514
replace gdenr_2012 = 5590 if year == 1997 & UCID_1934_2003 == 37514
replace gdenr_2018 = 5590 if year == 1997 & UCID_1934_2003 == 37514

replace gdenr = 5590 if year == 1991 & UCID_1934_2003 == 37516
replace gdenr_2012 = 5590 if year == 1991 & UCID_1934_2003 == 37516
replace gdenr_2018 = 5590 if year == 1991 & UCID_1934_2003 == 37516

replace gdenr = 5590 if year == 1999 & UCID_1934_2003 == 77010
replace gdenr_2012 = 5590 if year == 1999 & UCID_1934_2003 == 77010
replace gdenr_2018 = 5590 if year == 1999 & UCID_1934_2003 == 77010

replace gdenr = 5590 if year == 1999 & UCID_1934_2003 == 183467
replace gdenr_2012 = 5590 if year == 1999 & UCID_1934_2003 == 183467
replace gdenr_2018 = 5590 if year == 1999 & UCID_1934_2003 == 183467

replace gdenr = 5590 if year == 1986 & UCID_1934_2003 == 278501
replace gdenr_2012 = 5590 if year == 1986 & UCID_1934_2003 == 278501
replace gdenr_2018 = 5590 if year == 1986 & UCID_1934_2003 == 278501

replace gdenr = 5590 if year == 1997 & UCID_1934_2003 == 303617
replace gdenr_2012 = 5590 if year == 1997 & UCID_1934_2003 == 303617
replace gdenr_2018 = 5590 if year == 1997 & UCID_1934_2003 == 303617

replace gdenr = 5589 if year == 1995 & UCID_1934_2003 == 304416
replace gdenr_2012 = 5589 if year == 1995 & UCID_1934_2003 == 304416
replace gdenr_2018 = 5589 if year == 1995 & UCID_1934_2003 == 304416

replace gdenr = 5589 if year == 1991 & UCID_1934_2003 == 304813
replace gdenr_2012 = 5589 if year == 1991 & UCID_1934_2003 == 304813
replace gdenr_2018 = 5589 if year == 1991 & UCID_1934_2003 == 304813

replace gdenr = 5590 if year == 1999 & UCID_1934_2003 == 305430
replace gdenr_2012 = 5590 if year == 1999 & UCID_1934_2003 == 305430
replace gdenr_2018 = 5590 if year == 1999 & UCID_1934_2003 == 305430

replace gdenr = 5590 if UCID_1934_2003 == 305472 & year == 1962
replace gdenr_2012 = 5590 if UCID_1934_2003 == 305472 & year == 1962
replace gdenr_2018 = 5590 if UCID_1934_2003 == 305472 & year == 1962

replace gdenr = 5590 if UCID_1934_2003 == 305472 & year == 1964
replace gdenr_2012 = 5590 if UCID_1934_2003 == 305472 & year == 1964
replace gdenr_2018 = 5590 if UCID_1934_2003 == 305472 & year == 1964

replace gdenr = 5590 if UCID_1934_2003 == 305472 & year == 1965
replace gdenr_2012 = 5590 if UCID_1934_2003 == 305472 & year == 1965
replace gdenr_2018 = 5590 if UCID_1934_2003 == 305472 & year == 1965

replace gdenr = 5590 if UCID_1934_2003 == 305472 & year == 1990
replace gdenr_2012 = 5590 if UCID_1934_2003 == 305472 & year == 1990
replace gdenr_2018 = 5590 if UCID_1934_2003 == 305472 & year == 1990

replace gdenr = 5590 if year == 1991 & UCID_1934_2003 == 306290
replace gdenr_2012 = 5590 if year == 1991 & UCID_1934_2003 == 306290
replace gdenr_2018 = 5590 if year == 1991 & UCID_1934_2003 == 306290

replace gdenr = 5590 if UCID_1934_2003 == 307179 & year == 1992
replace gdenr_2012 = 5590 if UCID_1934_2003 == 307179 & year == 1992
replace gdenr_2018 = 5590 if UCID_1934_2003 == 307179 & year == 1992

replace gdenr = 5590 if UCID_1934_2003 == 307179 & year == 1998
replace gdenr_2012 = 5590 if UCID_1934_2003 == 307179 & year == 1998
replace gdenr_2018 = 5590 if UCID_1934_2003 == 307179 & year == 1998

replace gdenr = 5590 if UCID_1934_2003 == 307895 & year == 1996
replace gdenr_2012 = 5590 if UCID_1934_2003 == 307895 & year == 1996
replace gdenr_2018 = 5590 if UCID_1934_2003 == 307895 & year == 1996

replace gdenr = 5589 if UCID_1934_2003 == 309626 & year == 1989
replace gdenr_2012 = 5589 if UCID_1934_2003 == 309626 & year == 1989
replace gdenr_2018 = 5589 if UCID_1934_2003 == 309626 & year == 1989

replace gdenr = 5590 if UCID_1934_2003 == 309966 & year == 1984
replace gdenr_2012 = 5590 if UCID_1934_2003 == 309966 & year == 1984
replace gdenr_2018 = 5590 if UCID_1934_2003 == 309966 & year == 1984

replace gdenr = 5590 if UCID_1934_2003 == 310078 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 310078 & year == 1991 
replace gdenr_2018 = 5590 if UCID_1934_2003 == 310078 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 310826 & year == 1998
replace gdenr_2012 = 5590 if UCID_1934_2003 == 310826 & year == 1998
replace gdenr_2018 = 5590 if UCID_1934_2003 == 310826 & year == 1998

replace gdenr = 5590 if UCID_1934_2003 == 313309 & year == 1987
replace gdenr_2012 = 5590 if UCID_1934_2003 == 313309 & year == 1987
replace gdenr_2018 = 5590 if UCID_1934_2003 == 313309 & year == 1987

replace gdenr = 5590 if UCID_1934_2003 == 313483 & year == 1995
replace gdenr_2012 = 5590 if UCID_1934_2003 == 313483 & year == 1995
replace gdenr_2018 = 5590 if UCID_1934_2003 == 313483 & year == 1995

replace gdenr = 5590 if UCID_1934_2003 == 313722 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 313722 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 313722 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 313732 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 313732 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 313732 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 313753 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 313753 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 313753 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 313794 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 313794 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 313794 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 314050 & year == 1997
replace gdenr_2012 = 5590 if UCID_1934_2003 == 314050 & year == 1997
replace gdenr_2018 = 5590 if UCID_1934_2003 == 314050 & year == 1997

replace gdenr = 5590 if UCID_1934_2003 == 314117 & year == 1965
replace gdenr_2012 = 5590 if UCID_1934_2003 == 314117 & year == 1965
replace gdenr_2018 = 5590 if UCID_1934_2003 == 314117 & year == 1965

replace gdenr = 5590 if UCID_1934_2003 == 314396 & year == 1995
replace gdenr_2012 = 5590 if UCID_1934_2003 == 314396 & year == 1995
replace gdenr_2018 = 5590 if UCID_1934_2003 == 314396 & year == 1995

replace gdenr = 5590 if UCID_1934_2003 == 315081 & year == 1966
replace gdenr_2012 = 5590 if UCID_1934_2003 == 315081 & year == 1966
replace gdenr_2018 = 5590 if UCID_1934_2003 == 315081 & year == 1966

replace gdenr = 5589 if UCID_1934_2003 == 315946 & year == 1995
replace gdenr_2012 = 5589 if UCID_1934_2003 == 315946 & year == 1995
replace gdenr_2018 = 5589 if UCID_1934_2003 == 315946 & year == 1995

replace gdenr = 5589 if UCID_1934_2003 == 315950 & year == 1966
replace gdenr_2012 = 5589 if UCID_1934_2003 == 315950 & year == 1966
replace gdenr_2018 = 5589 if UCID_1934_2003 == 315950 & year == 1966

replace gdenr = 5589 if UCID_1934_2003 == 315961 & year == 1988
replace gdenr_2012 = 5589 if UCID_1934_2003 == 315961 & year == 1988
replace gdenr_2018 = 5589 if UCID_1934_2003 == 315961 & year == 1988

replace gdenr = 5589 if UCID_1934_2003 == 315965 & year == 1964
replace gdenr_2012 = 5589 if UCID_1934_2003 == 315965 & year == 1964
replace gdenr_2018 = 5589 if UCID_1934_2003 == 315965 & year == 1964

replace gdenr = 5589 if UCID_1934_2003 == 315966 & year == 1988
replace gdenr_2012 = 5589 if UCID_1934_2003 == 315966 & year == 1988
replace gdenr_2018 = 5589 if UCID_1934_2003 == 315966 & year == 1988

replace gdenr = 5589 if UCID_1934_2003 == 315969 & year == 2003
replace gdenr_2012 = 5589 if UCID_1934_2003 == 315969 & year == 2003
replace gdenr_2018 = 5589 if UCID_1934_2003 == 315969 & year == 2003

replace gdenr = 5589 if UCID_1934_2003 == 315977 & year == 1987
replace gdenr_2012 = 5589 if UCID_1934_2003 == 315977 & year == 1987
replace gdenr_2018 = 5589 if UCID_1934_2003 == 315977 & year == 1987

replace gdenr = 5589 if UCID_1934_2003 == 315985 & year == 2000
replace gdenr_2012 = 5589 if UCID_1934_2003 == 315985 & year == 2000
replace gdenr_2018 = 5589 if UCID_1934_2003 == 315985 & year == 2000

replace gdenr = 5589 if UCID_1934_2003 == 315996 & year == 1990
replace gdenr_2012 = 5589 if UCID_1934_2003 == 315996 & year == 1990
replace gdenr_2018 = 5589 if UCID_1934_2003 == 315996 & year == 1990

replace gdenr = 5589 if UCID_1934_2003 == 315996 & year == 2002
replace gdenr_2012 = 5589 if UCID_1934_2003 == 315996 & year == 2002
replace gdenr_2018 = 5589 if UCID_1934_2003 == 315996 & year == 2002

replace gdenr = 5589 if UCID_1934_2003 == 316006 & year == 1980
replace gdenr_2012 = 5589 if UCID_1934_2003 == 316006 & year == 1980
replace gdenr_2018 = 5589 if UCID_1934_2003 == 316006 & year == 1980

replace gdenr = 5589 if UCID_1934_2003 == 316010 & year == 1965
replace gdenr_2012 = 5589 if UCID_1934_2003 == 316010 & year == 1965
replace gdenr_2018 = 5589 if UCID_1934_2003 == 316010 & year == 1965

replace gdenr = 5589 if UCID_1934_2003 == 316029 & year == 1964
replace gdenr_2012 = 5589 if UCID_1934_2003 == 316029 & year == 1964
replace gdenr_2018 = 5589 if UCID_1934_2003 == 316029 & year == 1964

replace gdenr = 5589 if UCID_1934_2003 == 316032 & year == 1962
replace gdenr_2012 = 5589 if UCID_1934_2003 == 316032 & year == 1962 
replace gdenr_2018 = 5589 if UCID_1934_2003 == 316032 & year == 1962

replace gdenr = 5589 if UCID_1934_2003 == 316040 & year == 1995
replace gdenr_2012 = 5589 if UCID_1934_2003 == 316040 & year == 1995
replace gdenr_2018 = 5589 if UCID_1934_2003 == 316040 & year == 1995

replace gdenr = 5589 if UCID_1934_2003 == 316078 & year == 1990
replace gdenr_2012 = 5589 if UCID_1934_2003 == 316078 & year == 1990
replace gdenr_2018 = 5589 if UCID_1934_2003 == 316078 & year == 1990

replace gdenr = 5589 if UCID_1934_2003 == 316150 & year == 1963
replace gdenr_2012 = 5589 if UCID_1934_2003 == 316150 & year == 1963
replace gdenr_2018 = 5589 if UCID_1934_2003 == 316150 & year == 1963

replace gdenr = 5589 if UCID_1934_2003 == 316150 & year == 1965
replace gdenr_2012 = 5589 if UCID_1934_2003 == 316150 & year == 1965
replace gdenr_2018 = 5589 if UCID_1934_2003 == 316150 & year == 1965

replace gdenr = 5590 if UCID_1934_2003 == 316166 & year == 1992
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316166 & year == 1992
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316166 & year == 1992

replace gdenr = 5590 if UCID_1934_2003 == 316169 & year == 1980
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316169 & year == 1980
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316169 & year == 1980

replace gdenr = 5590 if UCID_1934_2003 == 316187 & year == 1990
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316187 & year == 1990
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316187 & year == 1990

replace gdenr = 5590 if UCID_1934_2003 == 316189 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316189 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316189 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 316190 & year == 1983
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316190 & year == 1983
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316190 & year == 1983

replace gdenr = 5590 if UCID_1934_2003 == 316195 & year == 1997
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316195 & year == 1997
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316195 & year == 1997

replace gdenr = 5590 if UCID_1934_2003 == 316200 & year == 1987
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316200 & year == 1987
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316200 & year == 1987

replace gdenr = 5590 if UCID_1934_2003 == 316201 & year == 1992
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316201 & year == 1992
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316201 & year == 1992

replace gdenr = 5590 if UCID_1934_2003 == 316206 & year == 1984
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316206 & year == 1984
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316206 & year == 1984

replace gdenr = 5590 if UCID_1934_2003 == 316216 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316216 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316216 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 316217 & year == 1990
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316217 & year == 1990
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316217 & year == 1990

replace gdenr = 5590 if UCID_1934_2003 == 316238 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316238 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316238 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 316239 & year == 1993
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316239 & year == 1993
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316239 & year == 1993

replace gdenr = 5590 if UCID_1934_2003 == 316260 & year == 1998
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316260 & year == 1998
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316260 & year == 1998

replace gdenr = 5590 if UCID_1934_2003 == 316264 & year == 1987
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316264 & year == 1987
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316264 & year == 1987

replace gdenr = 5590 if UCID_1934_2003 == 316286 & year == 1995
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316286 & year == 1995
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316286 & year == 1995

replace gdenr = 5590 if UCID_1934_2003 == 316300 & year == 1981
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316300 & year == 1981
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316300 & year == 1981

replace gdenr = 5590 if UCID_1934_2003 == 316301 & year == 1989
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316301 & year == 1989
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316301 & year == 1989

replace gdenr = 5590 if UCID_1934_2003 == 316306 & year == 1995
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316306 & year == 1995
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316306 & year == 1995

replace gdenr = 5590 if UCID_1934_2003 == 316318 & year == 1963
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316318 & year == 1963
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316318 & year == 1963

replace gdenr = 5590 if UCID_1934_2003 == 316318 & year == 1964
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316318 & year == 1964
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316318 & year == 1964

replace gdenr = 5590 if UCID_1934_2003 == 316326 & year == 1995
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316326 & year == 1995
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316326 & year == 1995

replace gdenr = 5590 if UCID_1934_2003 == 316335 & year == 1987
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316335 & year == 1987
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316335 & year == 1987

replace gdenr = 5590 if UCID_1934_2003 == 316336 & year == 1992
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316336 & year == 1992
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316336 & year == 1992

replace gdenr = 5590 if UCID_1934_2003 == 316344 & year == 1999
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316344 & year == 1999
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316344 & year == 1999

replace gdenr = 5590 if UCID_1934_2003 == 316358 & year == 2001
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316358 & year == 2001
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316358 & year == 2001

replace gdenr = 5590 if UCID_1934_2003 == 316377 & year == 1984
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316377 & year == 1984
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316377 & year == 1984

replace gdenr = 5590 if UCID_1934_2003 == 316379 & year == 1993
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316379 & year == 1993
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316379 & year == 1993

replace gdenr = 5590 if UCID_1934_2003 == 316403 & year == 1982
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316403 & year == 1982
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316403 & year == 1982

replace gdenr = 5590 if UCID_1934_2003 == 316418 & year == 1988
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316418 & year == 1988
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316418 & year == 1988

replace gdenr = 5590 if UCID_1934_2003 == 316418 & year == 1989
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316418 & year == 1989
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316418 & year == 1989

replace gdenr = 5590 if UCID_1934_2003 == 316426 & year == 1963
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316426 & year == 1963
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316426 & year == 1963

replace gdenr = 5590 if UCID_1934_2003 == 316430 & year == 1983
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316430 & year == 1983
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316430 & year == 1983

replace gdenr = 5590 if UCID_1934_2003 == 316453 & year == 1989
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316453 & year == 1989
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316453 & year == 1989

replace gdenr = 5590 if UCID_1934_2003 == 316483 & year == 1998
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316483 & year == 1998
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316483 & year == 1998

replace gdenr = 5590 if UCID_1934_2003 == 316497 & year == 1995
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316497 & year == 1995
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316497 & year == 1995

replace gdenr = 5590 if UCID_1934_2003 == 316497 & year == 1996
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316497 & year == 1996
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316497 & year == 1996

replace gdenr = 5590 if UCID_1934_2003 == 316508 & year == 1985
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316508 & year == 1985
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316508 & year == 1985

replace gdenr = 5590 if UCID_1934_2003 == 316513 & year == 1964
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316513 & year == 1964
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316513 & year == 1964

replace gdenr = 5590 if UCID_1934_2003 == 316513 & year == 1975
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316513 & year == 1975
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316513 & year == 1975

replace gdenr = 5590 if UCID_1934_2003 == 316513 & year == 1997
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316513 & year == 1997
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316513 & year == 1997

replace gdenr = 5590 if UCID_1934_2003 == 316516 & year == 1990
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316516 & year == 1990
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316516 & year == 1990

replace gdenr = 5590 if UCID_1934_2003 == 316517 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316517 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316517 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 316530 & year == 1984
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316530 & year == 1984
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316530 & year == 1984

replace gdenr = 5590 if UCID_1934_2003 == 316530 & year == 1987
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316530 & year == 1987
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316530 & year == 1987

replace gdenr = 5590 if UCID_1934_2003 == 316531 & year == 1964
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316531 & year == 1964
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316531 & year == 1964

replace gdenr = 5590 if UCID_1934_2003 == 316543 & year == 1999
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316543 & year == 1999
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316543 & year == 1999

replace gdenr = 5590 if UCID_1934_2003 == 316543 & year == 2000
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316543 & year == 2000
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316543 & year == 2000

replace gdenr = 5590 if UCID_1934_2003 == 316545 & year == 1999
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316545 & year == 1999
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316545 & year == 1999

replace gdenr = 5590 if UCID_1934_2003 == 316550 & year == 1965
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316550 & year == 1965
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316550 & year == 1965

replace gdenr = 5590 if UCID_1934_2003 == 316561 & year == 1999
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316561 & year == 1999
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316561 & year == 1999

replace gdenr = 5590 if UCID_1934_2003 == 316573 & year == 1984
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316573 & year == 1984
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316573 & year == 1984

replace gdenr = 5590 if UCID_1934_2003 == 316600 & year == 1984
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316600 & year == 1984
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316600 & year == 1984

replace gdenr = 5590 if UCID_1934_2003 == 316616 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316616 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316616 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 316638 & year == 1965
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316638 & year == 1965
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316638 & year == 1965

replace gdenr = 5590 if UCID_1934_2003 == 316638 & year == 1966
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316638 & year == 1966
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316638 & year == 1966

replace gdenr = 5590 if UCID_1934_2003 == 316638 & year == 1999
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316638 & year == 1999
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316638 & year == 1999

replace gdenr = 5590 if UCID_1934_2003 == 316646 & year == 1997
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316646 & year == 1997 
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316646 & year == 1997

replace gdenr = 5590 if UCID_1934_2003 == 316650 & year == 1965
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316650 & year == 1965
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316650 & year == 1965

replace gdenr = 5590 if UCID_1934_2003 == 316652 & year == 1962
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316652 & year == 1962
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316652 & year == 1962

replace gdenr = 5590 if UCID_1934_2003 == 316655 & year == 1997
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316655 & year == 1997
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316655 & year == 1997

replace gdenr = 5590 if UCID_1934_2003 == 316663 & year == 1966
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316663 & year == 1966
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316663 & year == 1966

replace gdenr = 5590 if UCID_1934_2003 == 316664 & year == 1992
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316664 & year == 1992
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316664 & year == 1992

replace gdenr = 5590 if UCID_1934_2003 == 316673 & year == 1965
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316673 & year == 1965
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316673 & year == 1965

replace gdenr = 5590 if UCID_1934_2003 == 316689 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316689 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316689 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 316692 & year == 1994
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316692 & year == 1994
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316692 & year == 1994

replace gdenr = 5590 if UCID_1934_2003 == 316711 & year == 2000
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316711 & year == 2000
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316711 & year == 2000

replace gdenr = 5590 if UCID_1934_2003 == 316722 & year == 1962
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316722 & year == 1962
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316722 & year == 1962

replace gdenr = 5590 if UCID_1934_2003 == 316730 & year == 1997
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316730 & year == 1997
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316730 & year == 1997

replace gdenr = 5590 if UCID_1934_2003 == 316744 & year == 1963
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316744 & year == 1963
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316744 & year == 1963

replace gdenr = 5590 if UCID_1934_2003 == 316770 & year == 1990
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316770 & year == 1990
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316770 & year == 1990

replace gdenr = 5589 if UCID_1934_2003 == 316774 & year == 1965
replace gdenr_2012 = 5589 if UCID_1934_2003 == 316774 & year == 1965
replace gdenr_2018 = 5589 if UCID_1934_2003 == 316774 & year == 1965

replace gdenr = 5590 if UCID_1934_2003 == 316775 & year == 1990
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316775 & year == 1990
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316775 & year == 1990

replace gdenr = 5590 if UCID_1934_2003 == 316822 & year == 1983
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316822 & year == 1983
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316822 & year == 1983

replace gdenr = 5590 if UCID_1934_2003 == 316875 & year == 1963
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316875 & year == 1963
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316875 & year == 1963

replace gdenr = 5590 if UCID_1934_2003 == 316877 & year == 1965
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316877 & year == 1965
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316877 & year == 1965

replace gdenr = 5590 if UCID_1934_2003 == 316927 & year == 1962
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316927 & year == 1962
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316927 & year == 1962

replace gdenr = 5590 if UCID_1934_2003 == 316928 & year == 1975
replace gdenr_2012 = 5590 if UCID_1934_2003 == 316928 & year == 1975
replace gdenr_2018 = 5590 if UCID_1934_2003 == 316928 & year == 1975

replace gdenr = 5590 if UCID_1934_2003 == 322436 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 322436 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 322436 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 325014 & year == 1995
replace gdenr_2012 = 5590 if UCID_1934_2003 == 325014 & year == 1995
replace gdenr_2018 = 5590 if UCID_1934_2003 == 325014 & year == 1995

replace gdenr = 5590 if UCID_1934_2003 == 325200 & year == 1960
replace gdenr_2012 = 5590 if UCID_1934_2003 == 325200 & year == 1960
replace gdenr_2018 = 5590 if UCID_1934_2003 == 325200 & year == 1960

replace gdenr = 5590 if UCID_1934_2003 == 326718 & year == 1988
replace gdenr_2012 = 5590 if UCID_1934_2003 == 326718 & year == 1988
replace gdenr_2018 = 5590 if UCID_1934_2003 == 326718 & year == 1988

replace gdenr = 5590 if UCID_1934_2003 == 330980 & year == 1966
replace gdenr_2012 = 5590 if UCID_1934_2003 == 330980 & year == 1966
replace gdenr_2018 = 5590 if UCID_1934_2003 == 330980 & year == 1966

replace gdenr = 5590 if UCID_1934_2003 == 330980 & year == 1982
replace gdenr_2012 = 5590 if UCID_1934_2003 == 330980 & year == 1982
replace gdenr_2018 = 5590 if UCID_1934_2003 == 330980 & year == 1982

replace gdenr = 5590 if UCID_1934_2003 == 330983 & year == 1963
replace gdenr_2012 = 5590 if UCID_1934_2003 == 330983 & year == 1963
replace gdenr_2018 = 5590 if UCID_1934_2003 == 330983 & year == 1963

replace gdenr = 5590 if UCID_1934_2003 == 330989 & year == 1963
replace gdenr_2012 = 5590 if UCID_1934_2003 == 330989 & year == 1963
replace gdenr_2018 = 5590 if UCID_1934_2003 == 330989 & year == 1963

replace gdenr = 5590 if UCID_1934_2003 == 331013 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 331013 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 331013 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 331013 & year == 1997
replace gdenr_2012 = 5590 if UCID_1934_2003 == 331013 & year == 1997
replace gdenr_2018 = 5590 if UCID_1934_2003 == 331013 & year == 1997

replace gdenr = 5590 if UCID_1934_2003 == 331016 & year == 1988
replace gdenr_2012 = 5590 if UCID_1934_2003 == 331016 & year == 1988
replace gdenr_2018 = 5590 if UCID_1934_2003 == 331016 & year == 1988

replace gdenr = 5590 if UCID_1934_2003 == 331039 & year == 1979
replace gdenr_2012 = 5590 if UCID_1934_2003 == 331039 & year == 1979
replace gdenr_2018 = 5590 if UCID_1934_2003 == 331039 & year == 1979

replace gdenr = 5590 if UCID_1934_2003 == 331039 & year == 1985
replace gdenr_2012 = 5590 if UCID_1934_2003 == 331039 & year == 1985
replace gdenr_2018 = 5590 if UCID_1934_2003 == 331039 & year == 1985

replace gdenr = 5589 if UCID_1934_2003 == 331040 & year == 1960
replace gdenr_2012 = 5589 if UCID_1934_2003 == 331040 & year == 1960
replace gdenr_2018 = 5589 if UCID_1934_2003 == 331040 & year == 1960

replace gdenr = 5590 if UCID_1934_2003 == 331074 & year == 1994
replace gdenr_2012 = 5590 if UCID_1934_2003 == 331074 & year == 1994
replace gdenr_2018 = 5590 if UCID_1934_2003 == 331074 & year == 1994

replace gdenr = 6133 if UCID_1934_2003 == 331088 & year == 1962
replace gdenr_2012 = 6133 if UCID_1934_2003 == 331088 & year == 1962
replace gdenr_2018 = 6133 if UCID_1934_2003 == 331088 & year == 1962

replace gdenr = 5590 if UCID_1934_2003 == 331093 & year == 1986
replace gdenr_2012 = 5590 if UCID_1934_2003 == 331093 & year == 1986
replace gdenr_2018 = 5590 if UCID_1934_2003 == 331093 & year == 1986

replace gdenr = 5590 if UCID_1934_2003 == 331098 & year == 1982
replace gdenr_2012 = 5590 if UCID_1934_2003 == 331098 & year == 1982
replace gdenr_2018 = 5590 if UCID_1934_2003 == 331098 & year == 1982

replace gdenr = 5590 if UCID_1934_2003 == 355823 & year == 1998
replace gdenr_2012 = 5590 if UCID_1934_2003 == 355823 & year == 1998
replace gdenr_2018 = 5590 if UCID_1934_2003 == 355823 & year == 1998

replace gdenr = 5590 if UCID_1934_2003 == 390713 & year == 1998
replace gdenr_2012 = 5590 if UCID_1934_2003 == 390713 & year == 1998
replace gdenr_2018 = 5590 if UCID_1934_2003 == 390713 & year == 1998

replace gdenr = 5590 if UCID_1934_2003 == 390713 & year == 2003
replace gdenr_2012 = 5590 if UCID_1934_2003 == 390713 & year == 2003
replace gdenr_2018 = 5590 if UCID_1934_2003 == 390713 & year == 2003

replace gdenr = 5590 if UCID_1934_2003 == 392331 & year == 1984
replace gdenr_2012 = 5590 if UCID_1934_2003 == 392331 & year == 1984
replace gdenr_2018 = 5590 if UCID_1934_2003 == 392331 & year == 1984

replace gdenr = 5590 if UCID_1934_2003 == 392331 & year == 1986
replace gdenr_2012 = 5590 if UCID_1934_2003 == 392331 & year == 1986
replace gdenr_2018 = 5590 if UCID_1934_2003 == 392331 & year == 1986

replace gdenr = 5590 if UCID_1934_2003 == 392331 & year == 1987
replace gdenr_2012 = 5590 if UCID_1934_2003 == 392331 & year == 1987
replace gdenr_2018 = 5590 if UCID_1934_2003 == 392331 & year == 1987

replace gdenr = 5590 if UCID_1934_2003 == 392834 & year == 1963
replace gdenr_2012 = 5590 if UCID_1934_2003 == 392834 & year == 1963
replace gdenr_2018 = 5590 if UCID_1934_2003 == 392834 & year == 1963

replace gdenr = 5590 if UCID_1934_2003 == 393452 & year == 1997
replace gdenr_2012 = 5590 if UCID_1934_2003 == 393452 & year == 1997
replace gdenr_2018 = 5590 if UCID_1934_2003 == 393452 & year == 1997

replace gdenr = 5590 if UCID_1934_2003 == 393726 
replace gdenr_2012 = 5590 if UCID_1934_2003 == 393726 
replace gdenr_2018 = 5590 if UCID_1934_2003 == 393726

replace gdenr = 5590 if UCID_1934_2003 == 393893 & year == 1962
replace gdenr_2012 = 5590 if UCID_1934_2003 == 393893 & year == 1962
replace gdenr_2018 = 5590 if UCID_1934_2003 == 393893 & year == 1962

replace gdenr = 6133 if UCID_1934_2003 == 393894 & year == 1966
replace gdenr_2012 = 6133 if UCID_1934_2003 == 393894 & year == 1966
replace gdenr_2018 = 6133 if UCID_1934_2003 == 393894 & year == 1966

replace gdenr = 5590 if UCID_1934_2003 == 394216 & year == 1964
replace gdenr_2012 = 5590 if UCID_1934_2003 == 394216 & year == 1964
replace gdenr_2018 = 5590 if UCID_1934_2003 == 394216 & year == 1964

replace gdenr = 5590 if UCID_1934_2003 == 394216 & year == 1982
replace gdenr_2012 = 5590 if UCID_1934_2003 == 394216 & year == 1982
replace gdenr_2018 = 5590 if UCID_1934_2003 == 394216 & year == 1982

replace gdenr = 5590 if UCID_1934_2003 == 394216 & year == 1996
replace gdenr_2012 = 5590 if UCID_1934_2003 == 394216 & year == 1996
replace gdenr_2018 = 5590 if UCID_1934_2003 == 394216 & year == 1996

replace gdenr = 5590 if UCID_1934_2003 == 394620 & year == 1995
replace gdenr_2012 = 5590 if UCID_1934_2003 == 394620 & year == 1995
replace gdenr_2018 = 5590 if UCID_1934_2003 == 394620 & year == 1995

replace gdenr = 5590 if UCID_1934_2003 == 396141 & year == 1995
replace gdenr_2012 = 5590 if UCID_1934_2003 == 396141 & year == 1995
replace gdenr_2018 = 5590 if UCID_1934_2003 == 396141 & year == 1995

replace gdenr = 5590 if UCID_1934_2003 == 396141 & year == 1997
replace gdenr_2012 = 5590 if UCID_1934_2003 == 396141 & year == 1997
replace gdenr_2018 = 5590 if UCID_1934_2003 == 396141 & year == 1997

replace gdenr = 5590 if UCID_1934_2003 == 396146 & year == 1991
replace gdenr_2012 = 5590 if UCID_1934_2003 == 396146 & year == 1991
replace gdenr_2018 = 5590 if UCID_1934_2003 == 396146 & year == 1991

replace gdenr = 5590 if UCID_1934_2003 == 396186 & year == 1984
replace gdenr_2012 = 5590 if UCID_1934_2003 == 396186 & year == 1984
replace gdenr_2018 = 5590 if UCID_1934_2003 == 396186 & year == 1984

replace gdenr = 5590 if UCID_1934_2003 == 396690 & year == 1998
replace gdenr_2012 = 5590 if UCID_1934_2003 == 396690 & year == 1998
replace gdenr_2018 = 5590 if UCID_1934_2003 == 396690 & year == 1998

replace gdenr = 5590 if UCID_1934_2003 == 400225 & year == 1984
replace gdenr_2012 = 5590 if UCID_1934_2003 == 400225 & year == 1984
replace gdenr_2018 = 5590 if UCID_1934_2003 == 400225 & year == 1984

replace gdenr = 5590 if UCID_1934_2003 == 400229 & year == 1983
replace gdenr_2012 = 5590 if UCID_1934_2003 == 400229 & year == 1983
replace gdenr_2018 = 5590 if UCID_1934_2003 == 400229 & year == 1983

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_finale.dta", replace
