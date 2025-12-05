
******************************************************
* (0) Set directories
******************************************************

	clear
	set more off

	capture cd "C:\Dropbox\Incumbency Advantage"
	if _rc==0{
		global hauptpfad "C:\Dropbox\Projekt Nationalräte\1 Data"
		global datapath "C:\Dropbox\Projekt Nationalräte\1 Data\1_data"
		global path_original_data "C:\Dropbox\Projekt Nationalräte\1 Data\0_original_data"
		global path_fuzzy_data "C:\Dropbox\Projekt Nationalräte\1 Data\1_data\fuzzy merge datasets"
		global path_fuzzy "C:\Dropbox\Projekt Nationalräte\1 Data\3_fuzzy merge do files"

	}
	
    capture cd "D:\SchmidLu\Dropbox\Projekt Nationalräte"

	if _rc==0{ 
	global hauptpfad "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data"
		global datapath "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data"
		global path_original_data "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\0_original_data"
		global path_fuzzy_data "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data\fuzzy merge datasets"
		global path_fuzzy "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\3_fuzzy merge do files"
	}
	
	 *net install rdrobust, from(https://sites.google.com/site/rdpackages/rdrobust/stata) replace


******************************************************
* (A) Merge Todesdaten to NR dataset and generate variables
******************************************************
		
		use "$datapath\nationalraete_1931_2015.dta", clear
		
	
		bysort ID: egen ever_elected=max(elected)
		
		collapse death_info ever_elected (first) firstname name canton year, by(ID)
		
		
		tab death_info ever_elected
		ttest death_info, by(ever_elected)
		
		br if ever_elected==1 & death_info==0
		
		tab year if ever_elected==1 & death_info==0
		
		gen dead_2017=0
		replace dead_2017=1 if dday!=.
		
		

******************************************************
* (B) Generate running variable
******************************************************
	
use "$datapath\nationalraete_1931_2015.dta", clear

bysort canton year: egen totseat=sum(elected)
bysort canton year list: egen votes_rank=rank(votes), field

sort canton year list votes_rank
keep  canton year list votes_rank

sum totseat
local imax = r(max)
forvalues i = 1(1)`imax' {
generate Bruchzahl`i' = votes/`i'
}

reshape long Bruchzahl, i(party) j(number)
destring number, replace

* (iv) Assign seats to lists

egen Bruchzahl_rank=rank(Bruchzahl), field

gen Seat=0 
replace Seat=1 if Bruchzahl_rank<=totseat

gsort -Bruchzahl

bysort number: egen SumVotes=sum(votes)

* Running varibale: how many more votes party j would have needed/lost to capture a specific seat?

gen partymargin_ji = .

sum party
local jmax = r(max)

forvalues j = 1(1)`jmax' {
	gen votes_n`j' = Bruchzahl if party!=`j'   // take vector of bruchzahlen of other than your own party (votes-j)
	egen order_n`j' = rank(-votes_n`j')        // rank() assigns average if two have the same rank
	replace order_n`j' = int(order_n`j')       // int(rank()) gives the same absolute rank to all with same rank
	forvalues i = 1(1)`imax' {
		local k = `imax'-`i'+1
		gen help_votes_nj_rank_`j'`i' = votes_n`j' if order_n`j' == `k' 
		replace help_votes_nj_rank_`j'`i' = votes_n`j' if order_n`j' == `k'-1 & help_votes_nj_rank_`j'`i'  == .  // if two have the same rank, take the previous
		replace help_votes_nj_rank_`j'`i' = votes_n`j' if order_n`j' == `k'-2 & help_votes_nj_rank_`j'`i'  == .  // if three have the same rank, take the previous
		replace help_votes_nj_rank_`j'`i' = votes_n`j' if order_n`j' == `k'-3 & help_votes_nj_rank_`j'`i'  == .  // if four have the same rank, take the previous
		egen votes_nj_rank_`j'`i' = min(help_votes_nj_rank_`j'`i')
		drop help_votes_nj_rank_`j'`i'
		replace partymargin_ji = votes - `i' * votes_nj_rank_`j'`i' if party == `j' & number == `i'
	}
}

replace SumVotes=SumVotes-partymargin_ji

keep party number votes partymargin_ji Bruchzahl_rank SumVotes
gsort party number


timer off 1

timer list 1


		
		
		
		
	