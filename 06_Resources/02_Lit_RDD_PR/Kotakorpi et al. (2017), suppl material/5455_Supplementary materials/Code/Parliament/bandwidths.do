* Investigating the impact of bandwidth on results 

* quietly rdob_mod2 $OUTCOME x
* global mainband = r(h_opt) /* IK bandwidth from main specification */
bwplot $OUTCOME x "" "Bandwidth" "Euros/year"
* Impact of bandwidth on main estimate. Vertical line marks the Imbens-Kalyanaraman optimal bandiwidth. 
graph export "$OUTPATH/BandwidthRobustness - parliament.pdf", as(pdf) fontface(Helvetica) replace


				  
local depvarlist incumbent neverelected female KESK KOK SDP OTHER south ///
	  age voteshare nseats nvotes last1_3tulo past_avg

foreach depvar in `depvarlist' {
bwplot `depvar' x "`: var label `depvar''" " " " "
}

*
graph combine  $GRAPHPATH/bw_last1_3tulo.gph $GRAPHPATH/bw_past_avg.gph, ///
	title("Average earnings before the election") b2title("Bandwidth") rows(1) ycommon
graph export "$OUTPATH/Figure_parliament_bw_priorearnings.pdf", fontface(Helvetica) replace
graph export "$OUTPATH/Figure_parliament_bw_priorearnings.eps", fontface(Helvetica) replace

graph combine  $GRAPHPATH/bw_incumbent.gph $GRAPHPATH/bw_neverelected.gph $GRAPHPATH/bw_voteshare.gph ///
 $GRAPHPATH/bw_female.gph $GRAPHPATH/bw_age.gph $GRAPHPATH/bw_south.gph ///
 $GRAPHPATH/bw_KESK.gph $GRAPHPATH/bw_KOK.gph $GRAPHPATH/bw_SDP.gph  ///
 , b2title("Bandwidth") rows(3) xcommon
graph export "$OUTPATH/Figure_parliament_bw_background.pdf", fontface(Helvetica)  replace
graph export "$OUTPATH/Figure_parliament_bw_background.eps", fontface(Helvetica)  replace
