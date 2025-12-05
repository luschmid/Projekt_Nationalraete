*** Replication materials for:
*** Cox, Gary W., Jon H. Fiva, and Daniel M. Smith
*** "Measuring the Competitiveness of Elections"
*** Political Analysis

*set local directory as replication archive"
*cd "/Users/~/Cox-Fiva-Smith-PA-Replication"

***** packages
* ssc install binscatter

***** data files required
* dta/Turnout_District_Level  /* replication file from Cox-Fiva-Smith JoP 2016 */
* dta/AlternativeMargin       /* replication file from Cox-Fiva-Smith JoP 2016 */
* dta/AlternativeMargin2      /* replication file from Cox-Fiva-Smith JoP 2016 */
* dta/illustration            /* For figure 1*/
* dta/ESP77-05parties.dta     /* From Grofman-Selb */
* dta/CH71-03parties.dta 	  /* From Grofman-Selb */ 
* dta/495_Selects_CumulativeFile_Data_1971-2015_v1.0.dta  /* Swiss National Election Studies from FORS: https://forsbase.unil.ch/project/study-public-overview/14738/0/ */
* dta/620_CCS_Data_Wave1_v4.0.dta  /* The Comparative Candidates Survey (CCS) from FORS: https://forsbase.unil.ch/project/study-public-overview/15760/0/ */
* dta/1965.dta   /* Norwegian election survey 1965, available for download here: https://nsd.no/nsddata/serier/norske_valgundersokelser.html */
* dta/1969.dta   /* Norwegian election survey 1969, available for download here: https://nsd.no/nsddata/serier/norske_valgundersokelser.html */
* dta/1965_magnitude.dta /* district magnitude as of 1965 */

**** Additional do files: Distance.do  /* PROGRAM FILE FOR COMPUTING DISTANCE TO THRESHOLD */
**** Additional do files: Grofman-Selb-Original.do  /* ORIGINAL FILE FROM GROFMAN AND SELB */

**** Data preparation -- Spain and Switzerland ******
run do/Spain_GS.do   		/* Use dta/ESP77-05parties.dta to create GS measure */
run do/Spain_CFS.do  		/* Use dta/ESP77-05parties.dta to create traditional measure */
run do/Switzerland_GS.do 	/* Use use dta/CH71-03parties.dta to create GS measure */
run do/Switzerland_CFS.do 	/* Use use dta/CH71-03parties.dta to create traditional measure */
run do/Panel.do				/* Merge all files to create a panel of Spain and Switzerland */

**** Data preparation -- survey data *****
run do/ImportSurveySwi.do 
run do/ImportSurveyNor.do 
run do/ImportSurveyEur.do 

*** Figure 1 -- Effort, votes, and seats 
run do/Figure1.do 

***  Figure2 -- Voter survey with respondents from Norway and Switzerland
run do/Figure2.do

*** Figure 3 -- Candidate survey with respondents from European countries
run do/Figure3.do

*** Figure 4 -- Comparing across measurement decisions  
run do/Figure4.do  

*** Figure 5 -- Alternative measures of competitiveness and their relationship with district magnitude
run do/Figure5.do 

*** Figure 6 -- Alternative measures of competitiveness and their relationship with voter turnout
run do/Figure6.do 

*** Appendix Figures
run do/AppendixFigures.do

*** Appendix Tables
run do/AppendixTables.do
