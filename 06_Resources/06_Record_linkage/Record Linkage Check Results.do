*******************************
*                             *
* Intercoder reliability test *
*                             *
*******************************

* This do-file computes measures for the accuracy of the RecordLinkage Procedure performed in R
* coder reliability

* 1. Read-in dataset

*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte\1 Data\9 Record Linkage"
global path "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\9 Record Linkage"

import delimited "$path\Results\df_final.csv", delim(";") clear


* 2. Compare coding

egen id_num=group(id)

bysort idfinal: egen id_num_sd=sd(id_num)
bysort id_num: egen idfinal_sd=sd(idfinal)

drop different_coding
gen different_coding=0
replace different_coding=1 if (id_num_sd>0 | idfinal_sd>0) & !missing(idfinal_sd) & !missing(id_num_sd)

tab different_coding
