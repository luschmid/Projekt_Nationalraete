clear
local path="C:\Users\a0910846\Dropbox\FivaFolkeSorensonReplication" /* INSERT DIRECTORY HERE */
cap mkdir "`path'\tables"
cap mkdir "`path'\log"
cap mkdir "`path'\figures"
cap mkdir "`path'\figures\gph"
cd `path'\do

/* INSTALLLING STATA PACKAGES - if needed */
ssc install sutex

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
twoway (line seats votes) (line prop votes, lpattern(dash) lcolor(gray*0.75)), xscale(range(0 1)) xmlabel(0.167 0.500 0.833) legend(off) xline(0.157, lwidth(medium) lpattern(dash) lcolor(gray)) xline(0.177, lwidth(medium) lpattern(dash) lcolor(gray)) ytitle(Share of representatives for party X) xtitle(Share of votes for party X) graphregion(fcolor(white)) scheme(s2mono)
graph export "..\figures\Figure1.tif", replace
graph export "..\figures\Figure1.eps", replace

/* FIGURE 2 - PARTY POSITIONS */
run Partypositions.do

/* FIGURE 3 - RD PLOTS FISCAL POLICY */
run RDplots.do 

/* FIGURE 4 -  MEDIA ANALYSIS */
run Retriever_v2.do     

/* TABLE 1: Descriptive statistics on fiscal policy outcomes  */
run DescriptivesPolicy.do

/* TABLE 2: BASELINE RESULTS */
run Analysis.do

/* APPENDIX A: Descriptive statistics: Seat shares in the local council */
run DescriptivesParties.do

/* Online Appendix B */
run OnlineAppendixB.do

/* TABLE C.1 IS BASED ON IND. LEVEL SURVEY DATA AND IS NOT INCLUDED IN REPLICATION FILE */

/* TABLE C.2 FIRST STAGE */
run FirstStage

/* TABLE C.3: BASELINE W. DEMOGRAPHIC CONTROLS */
run AnalysisControls

/* TABLE C.4: Estimated effects of left-right index and right wing-majority on demographic control variables */
run PreDetermined

/* TABLE C.5 and C.6 */
run ExternalValidity 

/* FIGURE C.1 - RD PLOTS ON CONTROLS     */
run RDplotsControls.do

/* FIGURE C.2. - VARYING BANDWIDTH FOR TWO ALTERNATIVE PTAX VARAIBLES */
run RDrobust.do 

/* FIGURE C.3 - DENSITY TESTS */ 
run DensityTest.do

/* FIGURE C.4 - BALANCE ON VOTESHARES */ 
run RDplotsVoteshares.do

/* FIGURE C.5 - C.15 */
run RDplotsIndividualParties.do
