version 16
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schelkem\Dropbox\Projekt Nationalräte\"

*global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"

**************************************************
*** Import corrected first and last name files ***
**************************************************

*** First names *** 

gen count = .
gen firstname = ""
gen newfirstname = ""
gen newlastname = ""

save "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\firstname_unique_R1.dta", replace

foreach i in @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z { 
	import excel using "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\firstname_unique`i'_out.xlsx", clear firstrow
	gen serie = "`i'"
	cap gen newlastname = ""
	gen nonewlastname = _rc  // error code for capture (110: newlastname exists already)
	save "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\firstname_unique`i'.dta", replace
}
use "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\firstname_unique@.dta", clear
foreach i in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z { 
	append using "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\firstname_unique`i'.dta"
	save "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\firstname_unique_R1.dta", replace
}

** normalize coding
tab serie nonewlastname
/*
A - H, N, Y: no column "newlastname"
A - C, E, Y: other denomination (Name or lastname)
D, F, G, H, N: no corrections of lastnames in firstname-files
*/
replace newlastname = lastname if newlastname == "" & lastname != ""  & (serie == "A" | serie == "B" | serie == "C")
drop lastname
drop E // one comment in columne E (serie C): "Jean-Marc or Jean-Luc" not useful comment
replace newlastname = Name if newlastname == "" & Name != ""  & serie == "E"
drop Name
replace newlastname = D if newlastname == "" & D != ""  & serie == "Y"
drop D


** Clean special cases (e.g., "gen.", "dit", "dette")
gen firstname_rev = firstname
replace firstname_rev = newfirstname if newfirstname != ""

gen flag1 = 0
gen flag_help = ustrpos(firstname_rev,"gen.")   // find specific string position
replace flag1 = flag_help if flag_help != 0 & flag1 == 0
drop flag_help
gen flag_help = ustrpos(firstname_rev,"gen Hans")
replace flag1 = flag_help if flag_help != 0 & flag1 == 0
drop flag_help
gen flag_help = ustrpos(firstname_rev,"gen Mar")
replace flag1 = flag_help if flag_help != 0 & flag1 == 0
drop flag_help
gen flag_help = ustrpos(firstname_rev,"gen Ruedi")
replace flag1 = flag_help if flag_help != 0 & flag1 == 0
drop flag_help
gen flag_help = ustrpos(firstname_rev,"dit ")
replace flag1 = flag_help if flag_help != 0 & flag1 == 0
drop flag_help
gen flag_help = ustrpos(firstname_rev,"dette ")
replace flag1 = flag_help if flag_help != 0 & flag1 == 0
drop flag_help
gen flag_help = ustrpos(firstname_rev,"detto ")
replace flag1 = flag_help if flag_help != 0 & flag1 == 0
drop flag_help
gen flag_help = ustrpos(firstname_rev," di ")
replace flag1 = flag_help if flag_help != 0 & flag1 == 0
drop flag_help

tab flag1

* drop everything after some attribute (e.g. "gen.")
gen firstname_rev1 = firstname_rev
replace firstname_rev1 = usubstr(firstname_rev, 1, ustrpos(firstname_rev, "gen.") - 1)  if ustrpos(firstname_rev, "gen.")
replace firstname_rev1 = usubstr(firstname_rev, 1, ustrpos(firstname_rev, "Gen.") - 1)  if ustrpos(firstname_rev, "Gen.")
replace firstname_rev1 = usubstr(firstname_rev, 1, ustrpos(firstname_rev, "gen Hans") - 1) if ustrpos(firstname_rev, "gen Hans")
replace firstname_rev1 = usubstr(firstname_rev, 1, ustrpos(firstname_rev, "gen Mar") - 1) if ustrpos(firstname_rev, "gen Mar")
replace firstname_rev1 = usubstr(firstname_rev, 1, ustrpos(firstname_rev, "gen Ruedi") - 1) if ustrpos(firstname_rev, "gen Ruedi")
replace firstname_rev1 = usubstr(firstname_rev, 1, ustrpos(firstname_rev, "dit ") - 1) if ustrpos(firstname_rev, "dit ")
replace firstname_rev1 = usubstr(firstname_rev, 1, ustrpos(firstname_rev, "dette ") - 1) if ustrpos(firstname_rev, "dette ")
replace firstname_rev1 = usubstr(firstname_rev, 1, ustrpos(firstname_rev, "detto ") - 1) if ustrpos(firstname_rev, "detto ")
replace firstname_rev1 = usubstr(firstname_rev, 1, ustrpos(firstname_rev, " di ") - 1) if ustrpos(firstname_rev, " di ")
replace firstname_rev1 = usubinstr(firstname_rev1, "(","",.)
replace firstname_rev1 = ustrtrim(firstname_rev1)

* correct mistakes
replace firstname_rev1 = "Bernadette H.G." if firstname_rev1 == "Berna" & firstname_rev == "Bernadette H.G."
replace firstname_rev1 = "Bernadette Hedwig" if firstname_rev1 == "Berna" & firstname_rev == "Bernadette Hedwig"
replace firstname_rev1 = "Bernadette Hélène" if firstname_rev1 == "Berna" & firstname_rev == "Bernadette Hélène"
replace firstname_rev1 = "Bernadette Josephine" if firstname_rev1 == "Berna" & firstname_rev == "Bernadette Josephine"
replace firstname_rev1 = "Bernadette L." if firstname_rev1 == "Berna" & firstname_rev == "Bernadette L."
replace firstname_rev1 = "Bernadette Maria" if firstname_rev1 == "Berna" & firstname_rev == "Bernadette Maria"
replace firstname_rev1 = "Bernadette T." if firstname_rev1 == "Berna" & firstname_rev == "Bernadette T."
replace firstname_rev1 = "Bernadette Theresia" if firstname_rev1 == "Berna" & firstname_rev == "Bernadette Theresia"
replace firstname_rev1 = "Lisa" if firstname_rev1 == "" & firstname_rev == "gen. Lisa"
replace firstname_rev1 = "René" if firstname_rev1 == "" & firstname_rev == "gen. René"
replace firstname_rev1 = "Odette Alice" if firstname_rev == "Odette Alice"
replace firstname_rev1 = "Odette D." if firstname_rev == "Odette D."
replace firstname_rev1 = "Odette F." if firstname_rev == "Odette F."
replace firstname_rev1 = "Odette Hélène" if firstname_rev == "Odette Hélène"
replace firstname_rev1 = "Odette Marcelle" if firstname_rev == "Odette Marcelle"
replace firstname_rev1 = "Odette de" if firstname_rev == "Odette de"
replace firstname_rev1 = "" if firstname_rev == "Yor detto Milano Yor"
replace firstname_rev1 = "Elisabeth" if firstname_rev == "Rotz Elisabeth gen. Lis"
replace newlastname = "Rotz" if firstname_rev == "Rotz Elisabeth gen. Lis"
replace firstname_rev1 = "Josef Anton" if firstname_rev == "Rotz Josef Anton gen. Jo"
replace newlastname = "Rotz" if firstname_rev == "Rotz Josef Anton gen. Jo"
replace firstname_rev1 = "Domenico" if firstname_rev == "DZucchi Domenico detto Sergio"
replace newlastname = "Zucchi" if firstname_rev == "DZucchi Domenico detto Sergio"
replace firstname_rev1 = "Ulrike" if firstname_rev == "Zimmermann Ulrike gen. UIl"
replace newlastname = "Zimmermann" if firstname_rev == "Zimmermann Ulrike gen. UIl"
replace firstname_rev1 = "Ulrike" if firstname_rev == "Zimmermann Ulrike gen. Uli"
replace newlastname = "Zimmermann" if firstname_rev == "Zimmermann Ulrike gen. Uli"
replace firstname_rev1 = "Efisio Cao di San" if firstname_rev == "Efisio Cao di San"
replace firstname_rev1 = "Andrea" if firstname_rev == "Galamini di Recanati Andrea"
replace newlastname = "Galamini di Recanati" if firstname_rev == "Galamini di Recanati Andrea"
replace firstname_rev1 = "Gertrud" if firstname_rev == "Ge. rud gen. Maria Petra"
replace firstname_rev1 = "Georg" if firstname_rev == "Geigy Georg gen. Georges"
replace newlastname = "Geigy" if firstname_rev == "Geigy Georg gen. Georges"
replace firstname_rev1 = "Alois" if firstname_rev == "Gruber Alois gen. Louis"
replace newlastname = "Gruber" if firstname_rev == "Gruber Alois gen. Louis"
replace firstname_rev1 = "Georges" if firstname_rev == "Guerin Georges dit Georgy"
replace newlastname = "Guerin" if firstname_rev == "Guerin Georges dit Georgy"
replace firstname_rev1 = "Johann P." if firstname_rev == "Götte Johann P. gen. Hans"
replace newlastname = "Götte" if firstname_rev == "Götte Johann P. gen. Hans"
replace firstname_rev1 = "HeiNick Heinrich detto Harry" if firstname_rev == "HeiNick Heinrich detto Harry"
replace firstname_rev1 = "Hildegard" if firstname_rev == "Ladner Hildegard gen. Angela"
replace newlastname = "Ladner" if firstname_rev == "Ladner Hildegard gen. Angela"
replace firstname_rev1 = "Mariei" if firstname_rev == "MLüdi Marie gen. Trudi"
replace newlastname = "Lüdi" if firstname_rev == "MLüdi Marie gen. Trudi"
replace firstname_rev1 = "Johann Ulrich" if firstname_rev == "Marti Johann Ulrich gen. Hans-Ueli"
replace newlastname = "Marti " if firstname_rev == "Marti Johann Ulrich gen. Hans-Ueli"
replace firstname_rev1 = "Luigi" if firstname_rev == "Matinoni Luigi detto Gigi"
replace newlastname = "Matinoni" if firstname_rev == "Matinoni Luigi detto Gigi"
replace firstname_rev1 = "Gertrud" if firstname_rev == "Meier Gertrud gen. Trudy"
replace newlastname = "Meier" if firstname_rev == "Meier Gertrud gen. Trudy"
replace firstname_rev1 = "Giuseppe" if firstname_rev == "Mondini Giuseppe detto Saverio"
replace newlastname = "Mondini" if firstname_rev == "Mondini Giuseppe detto Saverio"
replace firstname_rev1 = "Josef" if firstname_rev == "Moser Josef gen. Sepp"
replace newlastname = "Moser" if firstname_rev == "Moser Josef gen. Sepp"
replace firstname_rev1 = "Patricia" if firstname_rev == "Patricia gen. Zoebeli-Arnold"
replace newlastname = "Zoebeli-Arnold" if firstname_rev == "Patricia gen. Zoebeli-Arnold"
replace firstname_rev1 = "Pia icta gen. Schwester Bened" if firstname_rev == "Pia icta gen. Schwester Bened"
replace firstname_rev1 = "Richard dit Bressel-Braunschweiler René W." if firstname_rev == "Richard dit Bressel-Braunschweiler René W."  // this is correct! I checked on image!


* flag sen./jun./etc.
gen flag2 = 0
gen flag_help = ustrpos(firstname_rev, "sen.")
replace flag2 = 1 if flag_help != 0 & flag2 == 0
drop flag_help
gen flag_help = ustrpos(firstname_rev, "jun.")
replace flag2 = 1 if flag_help != 0 & flag2 == 0
drop flag_help
gen flag_help = ustrpos(firstname_rev, "Dr.")
replace flag2 = 1 if flag_help != 0 & flag2 == 0
drop flag_help
gen flag_help = ustrpos(firstname_rev, "Dott.")
replace flag2 = 1 if flag_help != 0 & flag2 == 0
drop flag_help
gen flag_help = ustrpos(firstname_rev, "Prof.")
replace flag2 = 1 if flag_help != 0 & flag2 == 0
drop flag_help

* delete sen./jun./etc.
replace firstname_rev1 = usubinstr(firstname_rev, "sen.", "",.) if flag2 == 1
replace firstname_rev1 = usubinstr(firstname_rev, "jun.", "",.) if flag2 == 1
replace firstname_rev1 = usubinstr(firstname_rev, "Dr.", "",.) if flag2 == 1
replace firstname_rev1 = usubinstr(firstname_rev, "Dott.", "",.) if flag2 == 1
replace firstname_rev1 = usubinstr(firstname_rev, "Prof.", "",.) if flag2 == 1
replace firstname_rev1 = usubinstr(firstname_rev, "Prof ", "",.) if flag2 == 1
replace firstname_rev1 = usubinstr(firstname_rev1, "  ", " ",.)
replace firstname_rev1 = ustrtrim(firstname_rev1)

* implement corrections in newfirstname
replace newfirstname = firstname_rev1 if flag1 != 0 | flag2 == 1
drop firstname_rev firstname_rev1
replace newfirstname = ustrtrim(newfirstname)
replace newlastname = ustrtrim(newlastname)

** indicators for newfirstname issues
gen indfirst = .
replace indfirst = 2 if newlastname != "" & newfirstname == ""  // lastname in firstname column
replace indfirst = 3 if newfirstname != ""  // corrections of only the first name (--> standard case)
replace indfirst = 4 if newfirstname != "" & newlastname != ""  // there is first and last name info in the firstname variable
replace indfirst = 1 if newfirstname == "0" | newfirstname == "1" | newfirstname == "?" | newfirstname == "Ortsfragment" | newfirstname == "ortsfragment" | newfirstname == "Namensfragment" | newfirstname == "namensfragment" | newfirstname == "Fragment" | newfirstname == "fragment" | newfirstname == "Firmeninformation" | newfirstname == "firmeninformation" | newfirstname == "Firmenbezeichnung" | newfirstname == "firmenbezeichnung" | newfirstname == "Firstname fragment" | newfirstname == "firstname fragment" | newfirstname == "Lastname fragment" | newfirstname == "lastname fragment"  | newlastname == "Lastname fragment" | newlastname == "lastname fragment" | newlastname == "lastnamefragment" // flagged cases by coders
label var indfirst "Indicator firstname issue: 1(fragment), 2(lastname in first), 3(corr firstname only - std), 4(corr first&last)"

replace newfirstname = "" if indfirst == 1 // these are flags
drop if firstname == "" & count == . // there are empty cells from excel-import

tab indfirst

drop flag1 flag2 nonewlastname

gen fid = _n
order fid
save "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\firstname_unique_R1.dta", replace

foreach i in @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z { 
	erase "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\firstname_unique`i'.dta"
}

** Read-out second-round checks
use "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\firstname_unique_R1.dta", clear
keep if indfirst == 2 | indfirst == 3

preserve
gen n = 1
collapse (sum) n, by(serie)
order serie n
export excel using "$path\02_Processed_data\10_Directors_1934_2003\03_Round2_Outfiles\overview_first.xlsx", first(var) replace 
restore

gen errorcode = .
drop count indfirst
order fid firstname newfirstname newlastname errorcode serie

foreach i in @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z {
	preserve
	keep if serie == "`i'"
	export excel using "$path\02_Processed_data\10_Directors_1934_2003\03_Round2_Outfiles\firstname_unique`i'_c2.xlsx", first(var) replace 
	restore
}


*** Lastnames ***

clear 
gen count = .
gen firstname = ""
gen newfirstname = ""
gen newlastname = ""

save "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\lastname_unique_R1.dta", replace

foreach i in @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z { 
	import excel using "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\lastname_unique`i'_out.xlsx", clear firstrow
	gen serie = "`i'"
	cap gen newfirstname = ""
	gen nonewfirstname = _rc  // error code for capture (110: newfirstname exists already)
	save "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\lastname_unique`i'.dta", replace
}
use "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\lastname_unique@.dta", clear
foreach i in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z {  
	append using "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\lastname_unique`i'.dta", force
	save "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\lastname_unique_R1.dta", replace
}

** normalize coding
tab serie nonewfirstname
/*
A - E, H, I, K - N, Q, U - X: no "newfirstname"
Y: column "D" with comments, but no "newfirstname" indication
W: "firstname" contains first names from lastname variable
*/
replace newfirstname = firstname if firstname != "" & serie == "W"
drop firstname

* indicators for newfirstname issues
gen indlast = .
replace indlast = 2 if newlastname == "" & newfirstname != ""  // firstname in lastname column
replace indlast = 3 if newlastname != ""  // corrections of only the last name (--> standard case)
replace indlast = 4 if newlastname != "" & newfirstname != ""  // there is first and last name info in the lastname variable
replace indlast = 1 if newlastname == "0" | newlastname == "1" | newlastname == "?" | newlastname == "Ortsfragment" | newlastname == "ortsfragment" | newlastname == "newlastname" | newlastname == "firmeninformation" | newlastname == "Firmeninformation" | newlastname == "Firmeninfortmation" | newlastname == "Firmenbezeichnung" | newlastname == "firmenbezeichnung" | newfirstname == "Firstname fragment" | newfirstname == "firstname fragment"  // flagged cases by coders
label var indlast "Indicator lastname issue: 1(fragment), 2(firstname in last), 3(corr lastname only - std), 4(corr first&last)"

replace newlastname = "" if indlast == 1  // these are flags
drop if lastname == "" & count == . // there are empty cells from excel-import
replace newfirstname = ustrtrim(newfirstname)
replace newlastname = ustrtrim(newlastname)

tab indlast
drop nonewfirstname

gen lid = _n

save "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\lastname_unique_R1.dta", replace

foreach i in @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z { 
	erase "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\lastname_unique`i'.dta"
}

** Read-out second-round checks
use "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\lastname_unique_R1.dta", clear
keep if indlast == 2 | indlast == 3

preserve
gen n = 1
collapse (sum) n, by(serie)
order serie n
export excel using "$path\02_Processed_data\10_Directors_1934_2003\03_Round2_Outfiles\overview_last.xlsx", first(var) replace 
restore

gen errorcode = .
drop count indlast
order lid lastname newlastname newfirstname errorcode serie

foreach i in @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z {
	preserve
	keep if serie == "`i'"
	export excel using "$path\02_Processed_data\10_Directors_1934_2003\03_Round2_Outfiles\lastname_unique`i'_c2.xlsx", first(var) replace 
	restore
}

