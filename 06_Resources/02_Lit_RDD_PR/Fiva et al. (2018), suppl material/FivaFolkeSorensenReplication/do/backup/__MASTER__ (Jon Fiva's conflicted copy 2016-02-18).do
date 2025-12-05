clear
cd C:\Users\A0910846\Dropbox\PartiesReplication\do

/* GRAPH PROG */
run "RD graph Programs.do"

/* DATA PREPARATION */
run PanelParties.do

/* SEAT ALLOCATIONS*/
run SeatAllocations.do /* this is MandatfordelningOlikaMetoderMedFelcheck.do and Index.do and StepPanel*/
run DistansKod.do /* creating variables for Folke (2014) method --- SeatAllocations.dta is the input */

/* MERGING POLITICAL AND FISCAL POLICY DATA */
/* the following is SpendingPanel_V2 and MergeData_2 */
run MergedData.do 	 /* Creating dependent variables  - based on Fiva, Halse and Natvik (2012) (http://www.jon.fiva.no/data.htm)  */ 

/* SIMULATIONS */
** run "coalition simulation norway ver comp".do  /* MODIFIED SAINT LAGUE METHOD */
** run "coalition simulation norway ver 2 DH".do  /* D'HONDT */

/* FIGURE 1 - ILLUSTRATION OF PR SYSTEM */

use ..\dta\IllustrerendeFigur.dta, clear
twoway (line seats votes) (line prop votes, lpattern(dash)), xscale(range(0 1)) xmlabel(0.1667 0.5000 0.8333) legend(off) xline(0.1567, lwidth(medium) lpattern(dash) lcolor(blue)) xline(0.1767, lwidth(medium) lpattern(dash) lcolor(blue)) ytitle(Share of representatives for party X) xtitle(Share of votes for party X) text(0.8 0.1 "0.1567", box) text(0.8 0.23 "0.1767", box)
graph export "..\figures\Figure1.eps", replace

/* FIGURE 2 - PARTY POSITIONS */
run Partypositions.do

/* FIGURE 3 -  MEDIA ANALYSIS */
run Retriever.do

/* TABLE 1: Descriptive statistics on fiscal policy outcomes  */
run DescriptivesPolicy.do

/* TABLE 2: BASELINE RESULTS */
run Analysis.do

/* FIGURE 4 - RD PLOTS FISCAL POLICY*/
run RDplots.do /* ****** THIS ONE--> run Coalitions_Figures_Binscatter        /* Figure 4 */ */

/* TABLE B.1: Descriptive statistics: Seat shares in the local council */
run DescriptivesParties.do

/* TABLE B.2 IS BASED ON IND. LEVEL SURVEY DATA AND IS NOT INCLUDED IN REPLICATION FILE */

/* TABLE B.3 FIRST STAGE */
run FirstStage

/* TABLE B.4: BASELINE W. DEMOGRAPHIC CONTROLS */
run AnalysisControls

/* TABLE B.5: Estimated effects of left-right index and right wing-majority on demographic control variables */
run PreDetermined

/* TABLE B.6 and B.7 */
run ExternalValidity /*   Table B.6 and B.7 */

/* FIGURE B.1 - RD PLOTS ON CONTROLS */
run RDplotsControls.do

/* FIGURE B.2. - VARYING BANDWIDTH FOR TWO ALTERNATIVE PTAX VARAIBLES */
run RDrobust.do /* --> Coalitions_RdRobust        /* FIGURE 5 */ */

/* FIGURE B.3 - BALANCE ON VOTESHARES */ 
run RDplotsVoteshares.do

/* FIGURE B.4 - DENSITY TESTS */ 
run DensityTest.do

/* CURRENT ONLINE APPENDIX C */
run RDplotsIndividualParties.do
