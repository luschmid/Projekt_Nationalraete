*INCUMBENCY EFFECT
eststo clear

quietly replace y=nextElected
quietly replace y=0 if y==. & year<2008
quietly rdob_mod2 y x
global band=r(h_opt)
quietly rd y x, bw($band) n(100) mbw(100)
eststo
estadd scalar Bandwidth = round($band,.001) 
global TITLE="Incumbency effect (Municipal elections)"
global XTITLE="Electoral closeness"
global YTITLE="Elected in subsequent election"
global B2TITLE=""
global GPHSTRING="$GRAPHPATH/IncumbencyEffect.gph"
MainRDgraph

esttab  using $OUTPATH/Incumbency.smcl, replace type label nonumber nodepvar se stats(N Bandwidth) mtitles("") ///
 title("Incumbency effect in municipal elections" ) 

 
* Bandwidth robustness
* Note: bwplot drops observations with missing $OUTCOME
clear
use $DATAPATH/tempdata
global oldOutcome = "$OUTCOME"
global OUTCOME = "nextElected"
bwplot nextElected x "Incumbency effect" "Bandwidth" " "
global OUTCOME=  "$oldOutcome"
global drop oldOutcome

