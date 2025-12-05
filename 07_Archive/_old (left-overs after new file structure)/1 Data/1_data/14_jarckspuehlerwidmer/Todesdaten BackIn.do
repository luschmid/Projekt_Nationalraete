

global path "D:\SchmidLu\Dropbox\Projekt Nationalräte\_old\1 Data\1_data\14_jarckspuehlerwidmer"
global path "C:\Users\Lukas\Dropbox\Projekt Nationalräte\_old\1 Data\1_data\14_jarckspuehlerwidmer"
*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte\1 Data\1_data\14_jarckspuehlerwidmer\Outfiles\"



* In the first round we have used three sources: 1. Historisches Lexikon der Schweiz, 
* 2. Todesanzeigenportal, 3. Wikipedia

* 1. Combine datasets

cd "$path\Infiles"
quietly fs *.csv 		// ssc install fs
display `r(files)'


foreach file in `r(files)'{
di "`file'"
import delimited using `file',  clear  delimiters(";")
gen coder="`file'"
sum bmonth
tostring geburtsdatumneu todesdatumneu dday, replace
capture append using "$path\Infiles\help.dta", force
save "$path\Infiles\help.dta", replace
}

erase "$path\Infiles\help.dta"

split coder,  parse("_")  
gen coder_id=substr(coder2,1,1)+substr(coder3,1,1)
gen batch=substr(coder4,1,1)

keep id geburtsdatumneu - quelletodesdatumneu unsicher coder coder_id batch
// goal: 27,604 obs = 26,644 unique obs + 960 doubletten

replace todesdatumneu= "lebend" if todesdatumneu=="lebendig" | todesdatumneu=="lebendug" | todesdatumneu=="Lebend" | todesdatumneu=="Noch lebend" | todesdatumneu=="noch lebend"
replace dday= "lebend" if dday=="Lebend" | dday=="Noch lebend" | dday=="noch lebend" | dday=="Lebend " 

gen lebend=.
replace lebend=1 if todesdatumneu=="lebend" | dday=="lebend"
destring dday, ignore("lebend") replace



* 2. Quality check of round 1


preserve
duplicates tag id, gen(dups)
tab dups
keep if dups!=0

replace coder_id="ys2" if coder=="geburtstodesdaten_yves_spuehler_5.csv" | coder=="geburtstodesdaten_yves_spuehler_6.csv" // Zivi package for Yves Spühler
drop coder

replace coder_id="zv" if coder_id=="fm" | coder_id=="rn" | coder_id=="ys2"  | coder_id=="jc" | coder_id=="jp"

drop geburtsdatumneu todesdatumneu

reshape wide   bday bmonth byear quellegeburtsdatumneu lebend  dday dmonth dyear quelletodesdatumneu unsicher batch, i(id) j(coder_id) string

foreach var in bday bmonth byear quellegeburtsdatumneu lebend dday dmonth dyear quelletodesdatumneu{
list id `var'jj `var'lw `var'ys `var'zv   if `var'jj!=. | `var'lw!=. | `var'ys!=. | `var'zv!=. 

}
restore


foreach var in bday bmonth byear lebend dday dmonth dyear{
bysort id: egen `var'sd=sd(`var')
list id `var' if !missing(`var'sd) & `var'sd>0
drop `var'sd
}

// Result: Overall good quality. If non-missing entry, figures are identical across coders. We modified seven entries. 

* 3. Keep most complete observations from duplicates 

egen nononmiss = rownonmiss(bday bmonth byear lebend dday dmonth dyear)
bysort id: egen maxnononmiss = max(nononmiss)
keep if nononmiss == maxnononmiss
duplicates drop id, force
drop nononmiss maxnononmiss

* 4. In-files for second round

preserve

cd "$path\Infiles2"
quietly fs 	*.xlsx	// ssc install fs
display `r(files)'


foreach file in `r(files)'{
display "`file'"
import excel  using `file',  clear first case(lower)
describe bday bmonth byear quellegeburtsdatumneu dmonth dyear quelletodesdatumneu
tostring geburtsdatumneu todesdatumneu dday, replace
replace geburtsdatumneu="" if geburtsdatumneu=="."
replace todesdatumneu="" if todesdatumneu=="."
replace dday="" if dday=="."
keep id geburtsdatumneu - quelletodesdatumneu
gen coder="`file'"
capture append using "$path\Infiles2\help.dta"
save "$path\Infiles2\help.dta", replace
}


split coder,  parse("_")  
gen coder_id=substr(coder2,1,1)+substr(coder3,1,1)
keep id geburtsdatumneu-quelletodesdatumneu  coder_id

gen lebend=.
replace lebend=1 if todesdatumneu=="lebend" | dday=="lebend"
destring dday, ignore("lebend") replace

save "$path\Infiles2\help.dta", replace
restore


append using "$path\Infiles2\help.dta"
erase "$path\Infiles2\help.dta"


* 5. Quality check of consistency between round 1 and 2

foreach var in bday bmonth byear lebend dday dmonth dyear{
bysort id: egen `var'sd=sd(`var')
list id `var' coder_id batch if !missing(`var'sd) & `var'sd>0
drop `var'sd
}

drop if id=="ZG-1999-0012" & coder_id=="jp"
drop if id=="ZH-1939-0066" & coder_id=="jc"
drop if id=="ZH-1943-0085" & coder_id=="jc"
drop if id=="ZH-1971-0062" & coder_id=="jc" & batch==""
drop if id=="ZH-1947-0074" & coder_id=="jc" 
drop if id=="ZH-1951-0163" & coder_id=="rn" 
drop if id=="ZH-1975-0493" & coder_id=="jc" 
drop if id=="ZG-2003-0015" & coder_id=="jp"
drop if id=="ZG-2011-0033" & coder_id=="jp"
drop if id=="ZH-1931-0040" & coder_id=="ys"
drop if id=="ZH-1931-0079" & coder_id=="jc"
drop if id=="ZH-1931-0039" & coder_id=="rn"
drop if id=="ZH-1935-9005" & coder_id=="jc"
drop if id=="ZH-1943-0221" & coder_id=="jc"
drop if id=="ZH-1947-0069" & coder_id=="jc"
drop if id=="ZH-1967-0245" & coder_id=="jc"
drop if id=="ZH-1987-0429" & coder_id=="jj"
drop if id=="ZH-1931-0019" & coder_id=="jc"
drop if id=="ZH-1931-0057" & coder_id=="jc"
drop if id=="ZH-1931-0067" & coder_id=="jj"
drop if id=="ZH-1943-0152" & coder_id=="jc"
drop if id=="GE-1935-0025" & coder_id=="jp" & batch==""
drop if id=="GE-1971-0027" & coder_id=="jp"
drop if id=="GL-2011-0002" & coder_id=="jp"
drop if id=="LU-1975-0018" & coder_id=="rn"
drop if id=="ZH-1955-0200" & coder_id=="jc"
drop if id=="ZH-1959-9002" & coder_id=="jc"
drop if id=="ZH-1971-0007" & coder_id=="jc"

foreach var in bday bmonth byear lebend dday dmonth dyear{
bysort id: egen `var'sd=sd(`var')
list id `var' coder_id batch if !missing(`var'sd) & `var'sd>0
drop `var'sd
}

egen nononmiss = rownonmiss(bday bmonth byear lebend dday dmonth dyear)
bysort id: egen maxnononmiss = max(nononmiss)
keep if nononmiss == maxnononmiss
duplicates drop id, force
drop maxnononmiss
drop if id == ""


/*
* 7. Generate outfiles for third round

* (a) Recode cases with death information

gen IndicatorDeathDate=0

foreach var of varlist todesdatumneu dday{
replace IndicatorDeathDate=1 if `var'!="" & `var'!="."
}

foreach var of varlist dmonth dyear{
replace IndicatorDeathDate=1 if `var'!=.
}


bysort id: egen IndicatorDeathDateMax=max(IndicatorDeathDate)

preserve
keep if IndicatorDeathDateMax==1
keep id
duplicates drop id, force
save "$path\IDsDeathInformation.dta", replace
restore

* (b) Outfile 


import excel  using "$path\JarckSpuehlerWidmer_out_wide.xlsx",  clear first case(lower)
merge 1:1 id using "$path\IDsDeathInformation.dta", gen(merge_deathdate)
erase "$path\IDsDeathInformation.dta"
drop if merge_deathdate==3 // drop all observations for which we have death date information from rounds 1 and 2

gen bday=.
gen bmonth=.
gen byear=.
gen dday=.
gen dmonth=.
gen dyear=.
gen unsicher=.

order id-geburtsdatumneu bday bmonth byear quellegeburtsdatumneu todesdatumneu dday dmonth dyear quelletodesdatumneu unsicher

export excel  using "$path\Zivis_out_wide.xlsx",  replace   firstrow(varlabels)

*/

* 8. Save death information of round 1 and 2

keep if nononmiss>0
drop nononmiss

save "$path\DeathInformation_Round_1_and_2.dta", replace

import excel  using "$path\JarckSpuehlerWidmer_out_wide.xlsx",  clear first case(lower)
drop geburtsdatumneu quellegeburtsdatumneu todesdatumneu quelletodesdatumneu
merge 1:1 id using "$path\DeathInformation_Round_1_and_2.dta", gen(merge_deathdate)
erase "$path\DeathInformation_Round_1_and_2.dta"
keep if merge_deathdate==3

gen name=name_1
replace name=name+", "+name_2 if name_2!=""
replace name=name+name_3 if name_3!=""
replace firstname=firstname+", "+h if h!=""
replace firstname=firstname+", "+i if i!=""

replace job=job+", "+q if q!=""
replace job=job+", "+r if r!=""
replace job=job+", "+s if s!=""
replace job=job+", "+t if t!=""
replace job=job+", "+u if u!=""
replace job=job+", "+v if v!=""
replace job=job+", "+w if w!=""
replace job=job+", "+x if x!=""
replace job=job+", "+y if y!=""

replace partyname=partyname+", "+aa if aa!=""
replace partyname=partyname+", "+ab if ab!=""
replace partyname=partyname+", "+ac if ac!=""
replace partyname=partyname+", "+ad if ad!=""

replace origin=origin+", "+af if af!=""
replace origin=origin+", "+ag if ag!=""
replace origin=origin+", "+ah if ah!=""
replace origin=origin+", "+ai if ai!=""
replace origin=origin+", "+aj if aj!=""
replace origin=origin+", "+ak if ak!=""
replace origin=origin+", "+al if al!=""

replace municipality=municipality+", "+an if an!=""
replace municipality=municipality+", "+ao if ao!=""
replace municipality=municipality+", "+ap if ap!=""
replace municipality=municipality+", "+aq if aq!=""

order id canton sex name firstname job partyname municipality origin birthyear bday bmonth byear quellegeburtsdatumneu lebend dday dmonth dyear quelletodesdatumneu  unsicher coder_id batch lebend

keep id canton sex name firstname job partyname municipality origin birthyear bday bmonth byear quellegeburtsdatumneu lebend dday dmonth dyear quelletodesdatumneu  unsicher coder_id batch

export excel  using "$path\GeburtsTodesdaten_Check_Round_1_and_2.xlsx",  replace   firstrow(varlabels)

// Note: This file was checked by Simon Lüchinger and saved as GeburtsTodesdaten_Check_Round_1_and_2_20181203.xslx.  This file is then merged to the main dataset in 01_Data Construction.do




* 9. Read in round 3

capture erase "$path\Infiles3\help.dta"

cd "$path\Infiles3"
quietly fs *.xlsx 		
display `r(files)'

quietly {
foreach file in `r(files)'{
di "`file'"
import excel using `file',  clear  first
//noisily tab dyear
destring bday bmonth byear dmonth dyear, replace
d, s
local obs = r(N)
local vars = r(k)
gen coder="`file'"
sum bmonth
local found = r(N)
tostring GeburtsdatumNeu TodesdatumNeu dday L, replace 
capture append using "$path\Infiles3\help.dta", force
save "$path\Infiles3\help.dta", replace
noisily di "`file'" _col(50) `obs' _col(70) `vars' _col(90) `found'
}
}

split coder,  parse("_")  
gen coder_id=substr(coder2,1,1)+substr(coder3,1,1)
gen batch=substr(coder4,1,1)

keep ID GeburtsdatumNeu - QuelleTodesdatumNeu unsicher coder coder_id batch

replace dday=strtrim(dday)
replace TodesdatumNeu=strtrim(TodesdatumNeu)
replace dday= "lebend 2017" if dday=="2017 lebend" 
replace dday= "lebend" if dday=="LEBEND" | dday=="Lebend" | dday=="Lebened"  

replace TodesdatumNeu=dday if substr(dday,1,6)=="lebend" 

replace dday="." if substr(dday,1,6)=="lebend" 

gen lebend=.
replace lebend=1 if substr(TodesdatumNeu,1,6)=="lebend"  & TodesdatumNeu!="lebend 1992" 

destring dday, ignore(".") replace


* 10. Quality check of consistency round 3

foreach var in bday bmonth byear lebend dday dmonth dyear{
bysort ID: egen `var'sd=sd(`var')
list ID `var' coder_id batch if !missing(`var'sd) & `var'sd>0
drop `var'sd
}

bysort ID: gen indi=_N
br if indi>1

* (a) First half up to sixth batch

drop if ID=="AG-1995-0002" & coder_id=="jc"
drop if ID=="AG-1999-0215" & coder_id=="jp"
drop if ID=="AG-2011-0141" & coder_id=="jp"
drop if ID=="BEJU-1943-0244" & coder_id=="jp"
drop if ID=="BEJU-1951-0084" & coder_id=="ys"
drop if ID=="BEJU-1971-0052" & coder_id=="jj"
drop if ID=="BEJU-1975-0496" & coder_id=="ys"
drop if ID=="BEJU-1983-0333" & coder_id=="jj"
drop if ID=="BEJU-1987-0056" & coder_id=="jj"
drop if ID=="BEJU-1987-0290" & coder_id=="jj"
drop if ID=="BEJU-1991-0164" & coder_id=="lw"
drop if ID=="BEJU-1999-0289" & coder_id=="lw"
drop if ID=="BEJU-2007-0082" & coder_id=="lw"
drop if ID=="BEJU-2007-0082" & coder_id=="lw"
drop if ID=="BEJU-2011-0108" & coder_id=="lw"
drop if ID=="BEJU-2011-0377" & coder_id=="lw"
drop if ID=="BEJU-2011-0431" & coder_id=="lw"
drop if ID=="BL-2015-0066" & coder_id=="ki"
drop if ID=="BS-1931-0035" & coder_id=="jc"
drop if ID=="BS-1971-0032" & coder_id=="ki"
drop if ID=="BS-2011-0101" & coder_id=="jc"
drop if ID=="BS-2015-0047" & coder_id=="jc"
drop if ID=="FR-1999-0037" & coder_id=="jc"
drop if ID=="GE-2015-0113" & coder_id=="ki"
drop if ID=="GR-1947-0010" & coder_id=="ki"
drop if ID=="GR-2003-0015" & coder_id=="ki"
drop if ID=="GR-2015-0011" & coder_id=="ki"
drop if ID=="LU-2007-0117" & coder_id=="jc"
drop if ID=="LU-2011-0089" & coder_id=="jc"
drop if ID=="NW-1995-0004" & coder_id=="jc"
drop if ID=="SG-1935-0043" & coder_id=="ki"
drop if ID=="SG-1987-0056" & coder_id=="ki"
drop if ID=="SG-2011-0049" & coder_id=="ki"
drop if ID=="SG-2015-0198" & coder_id=="jj"
drop if ID=="SH-1995-0006" & coder_id=="ki"
drop if ID=="SH-2007-9016" & coder_id=="ki"
drop if ID=="TI-2003-0023" & coder_id=="ki" & batch=="5"
drop if ID=="TI-2015-0047" & coder_id=="ki" & batch=="5"
drop if ID=="VD-1983-0034" & coder_id=="ki" & batch=="5"
drop if ID=="VD-2011-0269" & coder_id=="jj"
drop if ID=="VS-1999-9050" & coder_id=="jj"
drop if ID=="VS-2007-0110" & coder_id=="jj"
drop if ID=="VS-2015-0065" & coder_id=="jj"
drop if ID=="VS-2015-0086" & coder_id=="jj"
drop if ID=="ZG-2011-0035" & coder_id=="tm"
drop if ID=="ZH-1975-0480" & coder_id=="lw"
drop if ID=="ZH-1979-0026" & coder_id=="js"
drop if ID=="ZH-1979-0565" & coder_id=="js"
drop if ID=="ZH-1987-0137" & coder_id=="js"
drop if ID=="ZH-1991-0491" & coder_id=="js"
drop if ID=="ZH-1991-9238" & coder_id=="js"
drop if ID=="ZH-1995-0236" & coder_id=="js"
drop if ID=="ZH-1999-0341" & coder_id=="ki" & batch=="4"
drop if ID=="ZH-1999-0609" & coder_id=="ki" 
drop if ID=="ZH-2003-0817" & coder_id=="ki" & batch=="8"
drop if ID=="ZH-2003-0817" & coder_id=="js" 
drop if ID=="ZH-2007-0025" & coder_id=="ki"

replace dyear=1946 if ID=="BEJU-1943-0244" & coder_id=="jc"
replace lebend=1 if ID=="BEJU-1987-0056" & coder_id=="ys"
replace TodesdatumNeu="lebend" if ID=="BEJU-1987-0056" & coder_id=="ys"
replace lebend=1 if ID=="BEJU-1987-0290" & coder_id=="ys"
replace TodesdatumNeu="lebend" if ID=="BEJU-1987-0290" & coder_id=="ys"
replace byear=1950 if ID=="BEJU-1991-0164" & coder_id=="jp"  
replace QuelleGeburtsdatumNeu=6 if ID=="BEJU-1991-0164" & coder_id=="jp"  
replace byear=1985 if ID=="BEJU-2011-0108" & coder_id=="jp"  
replace QuelleGeburtsdatumNeu=6 if ID=="BEJU-2011-0108" & coder_id=="jp"  
replace byear=1955 if ID=="BEJU-2011-0431" & coder_id=="jp"  
replace QuelleGeburtsdatumNeu=6 if ID=="BEJU-2011-0431" & coder_id=="jp"  

egen nononmiss = rownonmiss(bday bmonth byear lebend dday dmonth dyear)
duplicates drop ID, force
drop if ID == ""

keep if nononmiss>0
drop nononmiss

save "$path\DeathInformation_Round_3.dta", replace

import excel  using "$path\JarckSpuehlerWidmer_out_wide.xlsx",  clear first case(lower)
gen ID=id
merge 1:1 ID using "$path\DeathInformation_Round_3.dta", gen(merge_deathdate)
erase "$path\DeathInformation_Round_3.dta"

keep if merge_deathdate==3

replace name=name+", "+e if e!=""
replace name=name+f if f!=""
replace firstname=firstname+", "+h if h!=""
replace firstname=firstname+", "+i if i!=""

replace job=job+", "+m if m!=""
replace job=job+", "+n if n!=""
replace job=job+", "+o if o!=""
replace job=job+", "+p if p!=""
replace job=job+", "+q if q!=""
replace job=job+", "+r if r!=""
replace job=job+", "+s if s!=""
replace job=job+", "+t if t!=""
replace job=job+", "+u if u!=""

replace partyname=partyname+", "+w if w!=""
replace partyname=partyname+", "+x if x!=""
replace partyname=partyname+", "+y if y!=""
replace partyname=partyname+", "+z if z!=""

replace origin=origin+", "+ab if ab!=""
replace origin=origin+", "+ac if ac!=""
replace origin=origin+", "+ad if ad!=""
replace origin=origin+", "+ae if ae!=""
replace origin=origin+", "+af if af!=""
replace origin=origin+", "+ag if ag!=""
replace origin=origin+", "+ah if ah!=""

replace municipality=municipality+", "+aj if aj!=""
replace municipality=municipality+", "+ak if ak!=""
replace municipality=municipality+", "+al if al!=""
replace municipality=municipality+", "+am if am!=""

order id canton sex name firstname job partyname municipality origin birthyear bday bmonth byear GeburtsdatumNeu bday bmonth byear QuelleGeburtsdatumNeu TodesdatumNeu dday dmonth dyear QuelleTodesdatumNeu unsicher coder_id batch lebend
keep id canton sex name firstname job partyname municipality origin birthyear bday bmonth byear GeburtsdatumNeu bday bmonth byear QuelleGeburtsdatumNeu TodesdatumNeu dday dmonth dyear QuelleTodesdatumNeu unsicher coder_id batch lebend

export excel  using "$path\GeburtsTodesdaten_Check_Round_3.xlsx",  replace   firstrow(varlabels)

// Note: This file was checked by Simon Lüchinger and saved as .xslx.  This file is then merged to the main dataset in 01_Data Construction.do


// out for Simon adding remaining files from round 3 (jj, js, ki)

preserve
import excel  using "$path\GeburtsTodesdaten_Check_Round_3_20181227.xlsx",  clear first case(lower) cellrange(A2:Z3245)
save "$path\help_luesi.dta", replace
restore

merge 1:1 id using "$path\help_luesi.dta", gen(merge_luesi)
erase "$path\help_luesi.dta"

keep if merge_luesi==1

order id	canton	sex	name	firstname	job	partyname	municipality	origin	birthyear	GeburtsdatumNeu	bday	bmonth	byear	QuelleGeburtsdatumNeu	TodesdatumNeu	lebend	dday	dmonth	dyear	QuelleTodesdatumNeu	unsicher	coder_id	batch	checked	mistakes
// Hinweis: Diese Daten werden in der Datei GeburtsTodesdaten_Check_Round_3.xlsx unten angefügt (31.1.2019)

* 11. Out for round 4, C and D

use "$path\Infiles3\help.dta", clear
capture erase "$path\Infiles3\help.dta"

split coder,  parse("_")  
gen coder_id=substr(coder2,1,1)+substr(coder3,1,1)
gen batch=substr(coder4,1,1)

destring dday, force replace

egen nononmiss = rownonmiss(bday bmonth byear dday dmonth dyear)
bysort ID: egen maxnononmiss = max(nononmiss)
duplicates drop ID, force
keep if nononmiss == maxnononmiss
drop maxnononmiss
drop if ID == ""

keep if nononmiss==0

tab coder_id
/*

   coder_id |      Freq.     Percent        Cum.
------------+-----------------------------------
         jc |      3,177       15.29       15.29
         jj |      4,950       23.82       39.12
         jp |      1,925        9.27       48.38
         ki |      5,798       27.91       76.29
         lw |      1,749        8.42       84.70
         tm |      2,284       10.99       95.70
         ys |        894        4.30      100.00
------------+-----------------------------------
      Total |     20,777      100.00

*/

gen round3b="C" if coder_id=="ki" 
replace round3b="D" if coder_id=="lw"

set seed 1234
generate random = runiform() 
replace round3b="C" if random<0.307 
replace round3b="D" if random>=0.307 & round3b==""

bysort round3b: tab coder_id

sort coder_id random // sort such that jc and jj are the first (avoid ki and lw)
gen id_help=_n  if _n<=216


sort id_help coder_id random round3b  
bysort round3b: egen id_help2=rank(random) if id_help==. 
bysort round3b: replace id_help=id_help2+216 if id_help==.

gen name=name_1
replace name=name+", "+name_2 if name_2!=""
replace firstname=firstname+", "+I if I!=""



replace job=job+", "+Y if Y!=""
replace job=job+", "+Z if Z!=""
replace job=job+", "+AA if AA!=""
replace job=job+", "+AB if AB!=""
replace job=job+", "+AC if AC!=""
replace job=job+", "+AD if AD!=""
replace job=job+", "+AE if AE!=""


replace partyname=partyname+", "+AI if AI!=""
replace partyname=partyname+", "+AJ if AJ!=""

replace origin=origin+", "+AN if AN!=""
replace origin=origin+", "+AO if AO!=""
replace origin=origin+", "+AP if AP!=""


replace municipality=municipality+", "+AV if AV!=""
replace municipality=municipality+", "+AW if AW!=""
replace municipality=municipality+", "+AX if AX!=""

rename Canton canton
rename Sex sex

order ID canton sex name firstname birthyear GeburtsdatumNeu bday bmonth byear QuelleGeburtsdatumNeu TodesdatumNeu dday dmonth dyear QuelleTodesdatumNeu job partyname municipality origin  id_help round3b
keep ID canton sex name firstname birthyear GeburtsdatumNeu bday bmonth byear QuelleGeburtsdatumNeu TodesdatumNeu dday dmonth dyear QuelleTodesdatumNeu job partyname municipality origin  id_help round3b

replace QuelleGeburtsdatumNeu=.
replace QuelleTodesdatumNeu=.

destring TodesdatumNeu GeburtsdatumNeu , replace


label var job "job"
label var partyname "partyname"
label var municipality "municipality"
label var origin "origin"
label var birthyear "birthyear"

// Note: We need 10 batches b/c di 9*976=8784 is slightly higher than 8798
/*
    round3b |      Freq.     Percent        Cum.
------------+-----------------------------------
          C |      8,838       50.11       50.11
          D |      8,798       49.89      100.00
------------+-----------------------------------
      Total |     17,636      100.00
*/

forvalues i=1(1)10{
preserve
keep if (inrange(id_help,(`i'-1)*24+1,(`i')*24)) |  ( round3b=="C" & inrange(id_help,(`i'-1)*976+217,(`i')*976+216))
drop id_help round3b
sort ID
export excel  using "$path\Outfiles3\GeburtsTodesdaten3_C_`i'.xlsx",  replace   firstrow(varlabels)
restore

preserve
keep if (inrange(id_help,(`i'-1)*24+1,(`i')*24)) |  ( round3b=="D" & inrange(id_help,(`i'-1)*976+217,(`i')*976+216))
drop id_help round3b
sort ID
export excel  using "$path\Outfiles3\GeburtsTodesdaten3_D_`i'.xlsx",  replace   firstrow(varlabels)
restore
}



* 12. Read in round 4 (without files by new zivis gonzalez and dellagiacoma)


capture erase "$path\Infiles4\help.dta"

cd "$path\Infiles4"
quietly fs *.xlsx 		
display `r(files)'

quietly {
foreach file in `r(files)'{
di "`file'"
import excel using `file',  clear  first
//noisily tab dyear
d, s
local obs = r(N)
local vars = r(k)
gen coder="`file'"
//noisily di "`file'" 
//noisily sum bday
sum bmonth
local found = r(N)
tostring GeburtsdatumNeu TodesdatumNeu QuelleGeburtsdatumNeu QuelleTodesdatumNeu bday bmonth byear dday dmonth dyear, replace 
capture append using "$path\Infiles4\help.dta", force
save "$path\Infiles4\help.dta", replace
noisily di "`file'" _col(50) `obs' _col(70) `vars' _col(90) `found'
}
}

drop if dday!="." & dday!=""
drop if dmonth!="." & dmonth!=""
drop if dyear!="." & dyear!=""


replace name=name_1 if name==""
replace name=name+", "+name_2 if name_2!=""


replace firstname=firstname+", "+I if I !=""

replace job=job+", "+ Y if Y!=""
replace job=job+", "+ Z if Z!=""

replace partyname=partyname+", "+AI if AI!=""


replace origin=origin+", "+AN if AN!=""
replace origin=origin+", "+AO if AO!=""
replace origin=origin+", "+AP if AP!=""

replace municipality=municipality+", "+AV if AV!=""
replace municipality=municipality+", "+AW if AW!=""

rename Canton canton
rename Sex sex

replace bday=""
replace bmonth=""
replace byear=""
replace QuelleGeburtsdatumNeu=""
replace GeburtsdatumNeu=""
replace dday=""
replace dmonth=""
replace dyear=""
replace QuelleTodesdatumNeu=""
replace TodesdatumNeu=""
replace Unsicher=.

drop if ID=="" & firstname==""



order ID canton sex name firstname birthyear GeburtsdatumNeu bday bmonth byear QuelleGeburtsdatumNeu TodesdatumNeu dday dmonth dyear QuelleTodesdatumNeu Unsicher partyname municipality job origin    

keep ID canton sex name firstname birthyear GeburtsdatumNeu bday bmonth byear QuelleGeburtsdatumNeu TodesdatumNeu dday dmonth dyear QuelleTodesdatumNeu Unsicher partyname municipality job origin   

set seed 1234
generate random = runiform() 
sort random
gen id_help=_n

forvalues i=1(1)16{
preserve
keep if inrange(id_help,(`i'-1)*1000+1,(`i')*1000)
drop random id_help 
sort ID
export excel  using "$path\Outfiles5\GeburtsTodesdaten5_`i'.xlsx",  replace   firstrow(varlabels)
restore
}
