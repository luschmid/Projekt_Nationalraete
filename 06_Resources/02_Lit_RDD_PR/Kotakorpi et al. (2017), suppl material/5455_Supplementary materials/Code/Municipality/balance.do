* BALANCE of background variables --- RD graphs

quietly rdob_mod2 $OUTCOME x
global mainband = r(h_opt) /* use IK bandwidth from main specification */
*global mainband = 30 /* use a fixed bw */
*global mainband = 10^10 /* use case-by-case IK bandwidth */
global forcebw = 1 /* forcebw=1 forces all bandwidths equal to $mainband, under
					  forcebw=0 $mainband is used as an upper bound for bandwidth on top of case-by-case bw */
global XTITLE=""
global YTITLE=""

eststo clear
					  
* Group 1 (binary variables)

foreach depvar in incumbent female KESK KOK SDP OTHER  {
quietly replace y=`depvar'
quietly rdob_mod2 y x if $OUTCOME!=.
global band0=r(h_opt)
if $forcebw {
	global band = $mainband
	} 
else {
	global band = min($band0,$mainband)
	}
quietly rd y x if $OUTCOME!=., bw($band) n(100) mbw(100)
quietly eststo
estadd scalar CaseBandwidth = round($band0,.001) 
estadd scalar Bandwidth = round($band,.001) 
global TITLE="`: var label `depvar''"
global GPHSTRING="$GRAPHPATH/`depvar'.gph"
MainRDgraph
}

esttab using $OUTPATH/RD_Balance_1.smcl, replace type label nodepvar nonumber se  ///
 stats(N Bandwidth CaseBandwidth)  ///
 mtitles("Incumbent" "Female" "Centre" "NCP" "SDP" "OTHER") ///
 title("Balance of predetermined variables - Group 1" )

eststo clear

* Group 2 ("continuous" variables)

foreach depvar in last1_3tulo past_avg voteshare age yob year nseats nvotes{
quietly replace y=`depvar'
quietly rdob_mod2 y x if $OUTCOME!=.
global band0=r(h_opt)
if $forcebw {
	global band = $mainband
	} 
else {
	global band = min($band0,$mainband)
	}
quietly rd y x if $OUTCOME!=., bw($band) n(100) mbw(100)
quietly eststo
estadd scalar CaseBandwidth = round($band0,.001) 
estadd scalar Bandwidth = round($band,.001) 
global TITLE="`: var label `depvar''"
global GPHSTRING="$GRAPHPATH/`depvar'.gph"
MainRDgraph
}

esttab using $OUTPATH/RD_Balance_2.smcl, replace type label nodepvar nonumber se  ///
 stats(N Bandwidth CaseBandwidth )  ///
 mtitles("last1_3tulo" "past_avg" "voteshare" "age" "yob"  "year" "nseats" "nvotes") ///
 title("Balance of predetermined variables - Group 2" )


graph combine  $GRAPHPATH/last1_3tulo.gph $GRAPHPATH/past_avg.gph, title("Average earnings before the election") rows(1) ycommon
graph export "$OUTPATH/Figure_muni_balance_priorearnings.pdf", fontface(Helvetica) replace
graph export "$OUTPATH/Figure_muni_balance_priorearnings.eps", fontface(Helvetica) replace

graph combine  $GRAPHPATH/incumbent.gph  $GRAPHPATH/voteshare.gph $GRAPHPATH/nseats.gph ///
 $GRAPHPATH/female.gph $GRAPHPATH/age.gph  ///
 $GRAPHPATH/KESK.gph $GRAPHPATH/KOK.gph $GRAPHPATH/SDP.gph  ///
, b2title("") rows(3) xcommon
graph export "$OUTPATH/Figure_muni_balance_background.pdf", fontface(Helvetica)  replace
graph export "$OUTPATH/Figure_muni_balance_background.eps", fontface(Helvetica)  replace
