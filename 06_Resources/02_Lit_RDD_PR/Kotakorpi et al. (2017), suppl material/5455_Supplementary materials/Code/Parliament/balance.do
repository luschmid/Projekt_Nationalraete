* Balance of background variables --- RD graphs of other falsification tests 

cap drop y
cap gen y = $OUTCOME

* Bandwidth options
quietly rdob_mod2 y x
global mainband=r(h_opt) /* IK-optimal bandwidth for main specification */
quietly if ($forceband!=.) global mainband=$forceband  /* forcing the main bandwidth */

global samebw = 1 /* samebw=1 forces all bandwidths to equal $mainband,
					 samebw=0 use case-by-case IK-optimal bw */ 
if $samebw==1 global casenote = "All cases estimated with bandwidth $mainband" 
if $samebw==0 global casenote = "Estimated with case-by-case IK optimal bandwidth"

*Options for MainRDgraph
global XTITLE=""
global YTITLE=""
global XMAX=50
global MARKERSIZE="tiny"

eststo clear
					  
* Group 1 (binary variables)

foreach depvar in incumbent neverelected female KESK KOK SDP OTHER south {
quietly replace y=`depvar'
quietly rdob_mod2 y x if $OUTCOME!=.
global ikband=r(h_opt)
if $samebw {
	global band = $mainband
	} 
else {
	global band = $ikband
	}
quietly rd y x if $OUTCOME!=., bw($band) n(100) mbw(100)
quietly eststo
estadd scalar CaseBandwidth = round($ikband,.001) 
estadd scalar Bandwidth = round($band,.001) 
global TITLE="`: var label `depvar''"
global GPHSTRING="$GRAPHPATH/`depvar'.gph"
MainRDgraph
}

esttab using $OUTPATH/RD_Balance_1.smcl, replace type label nodepvar nonumber se  ///
 stats(N Bandwidth CaseBandwidth)  ///
 mtitles("Incumbent" "Never elected" "Female" "Centre" "NCP" "SDP" "OTHER" "South") ///
 title("Balance of predetermined variables - Group 1" ) ///
 addnote("$casenote")

eststo clear

* corr incumbent neverelected female KESK KOK SDP OTHER south if $OUTCOME!=. 

* Group 2 ("continuous" variables)

foreach depvar in last1_3tulo past_avg voteshare age yob year nseats nvotes{
quietly replace y=`depvar'
quietly rdob_mod2 y x if $OUTCOME!=.
global ikband=r(h_opt)
if $samebw {
	global band = $mainband
	} 
else {
	global band = $ikband
	}
quietly rd y x if $OUTCOME!=., bw($band) n(100) mbw(100)
quietly eststo
estadd scalar CaseBandwidth = round($ikband,.001) 
estadd scalar Bandwidth = round($band,.001) 
global TITLE="`: var label `depvar''"
global GPHSTRING="$GRAPHPATH/`depvar'.gph"
MainRDgraph
}

esttab using $OUTPATH/RD_Balance_2.smcl, replace type label nodepvar nonumber se  ///
 stats(N Bandwidth CaseBandwidth )  ///
 mtitles("last1_3tulo" "past_avg" "voteshare" "age" "yob" "year" "nseats" "nvotes") ///
 title("Balance of predetermined variables - Group 2") ///
  addnote("$casenote")


graph combine  $GRAPHPATH/last1_3tulo.gph $GRAPHPATH/past_avg.gph, title("Average earnings before the election") rows(1) ycommon
graph export "$OUTPATH/Figure_parliament_balance_priorearnings.pdf", fontface(Helvetica) replace
graph export "$OUTPATH/Figure_parliament_balance_priorearnings.eps",  fontface(Helvetica) replace

graph combine  $GRAPHPATH/incumbent.gph  $GRAPHPATH/neverelected.gph $GRAPHPATH/voteshare.gph ///
 $GRAPHPATH/female.gph $GRAPHPATH/age.gph $GRAPHPATH/south.gph ///
 $GRAPHPATH/KESK.gph $GRAPHPATH/KOK.gph $GRAPHPATH/SDP.gph  ///
 ,b2title("") rows(3) xcommon
graph export "$OUTPATH/Figure_parliament_balance_background.pdf", fontface(Helvetica) replace
graph export "$OUTPATH/Figure_parliament_balance_background.eps", fontface(Helvetica) replace

