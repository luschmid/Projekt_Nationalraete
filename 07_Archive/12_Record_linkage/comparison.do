* 1. Set directories 

*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

insheet using "$path\02_Processed_data\12_Record_linkage\03_Data_analysis\bisnode_results_long_Generation_7_optimal.csv", ///
	delimiter(";") clear

drop if category=="Verwaltungsrat_Praesident"		
merge 1:1 id_0 year using "$path\02_Processed_data\12_Record_linkage\03_Data_analysis\data_SL.dta"