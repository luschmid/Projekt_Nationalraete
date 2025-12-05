clear all

*global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schelkem\Dropbox\Projekt Nationalräte\"
global path "C:\Users\schmutzy\Dropbox\Projekt Nationalräte\" /* Office */
global dump "C:\Users\schmutzy\Desktop\Projekt Nationalräte\" /* Office Dumpster */
*global path "C:\Users\yanni\Dropbox\Projekt Nationalräte\" /* Home */
*global dump "C:\Users\yanni\Desktop\Projekt Nationalräte\" /* Home Dumpster */

*
*Loading final panel data at the end of "sugarcube_gaps_281124.do"
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_finale_3.dta", clear
* Finale panel / raw (rematch with initial dataset to have clean names, etc.)

********************************************************************************
* Gen relocation, creation and deaths of firms

bysort UCID_1934_2003 (year):  gen occurences = cond(_N==1,0,_n) //nb of appearance by ID, creation when = 1
bysort UCID_1934_2003 (year):  egen max_occurences = max(occurences) //years of appearance, death when nb of appearance = max appearance

* Creation of firms
by UCID_1934_2003 (year), sort: gen creation_2012 = 1 if occurences == 1 & year != 1934
by UCID_1934_2003 (year), sort: gen creation_2018 = 1 if occurences == 1 & year != 1934
by UCID_1934_2003 (year), sort: replace creation_2012 = 1 if occurences == 0 & year != 1934
by UCID_1934_2003 (year), sort: replace creation_2018 = 1 if occurences == 0 & year != 1934

* Deaths (or disappearing from our data)
by UCID_1934_2003 (year), sort: gen death_2012 = 1 if occurences == max_occurences & year != 2003
by UCID_1934_2003 (year), sort: gen death_2018 = 1 if occurences == max_occurences & year != 2003

* Relocation of firms
by UCID_1934_2003 (year), sort: gen relocation_in_2012 = 1 if gdenr_2012 != gdenr_2012[_n-1] & creation_2012 != 1
by UCID_1934_2003 (year), sort: gen relocation_in_2018 = 1 if gdenr_2018 != gdenr_2018[_n-1] & creation_2018 != 1
replace relocation_in_2012 = 0 if year == 1934
replace relocation_in_2018 = 0 if year == 1934

by UCID_1934_2003 (year), sort: gen relocation_out_2012 = 1 if gdenr_2012 != gdenr_2012[_n+1] & max_occurences != occurences
by UCID_1934_2003 (year), sort: gen relocation_out_2018 = 1 if gdenr_2018 != gdenr_2018[_n+1] & max_occurences != occurences


*Dummy to count firms with a certain nb of owners
		replace owners = usubinstr(owners, ";", ",",.)
		replace owners = usubinstr(owners, "  ", " ",.)
		replace owners = usubinstr(owners, " ", "",.)
		drop nb_owners
		gen nb_owners = length(owners) - length(subinstr(owners, ",", "", .)) + 1
		gen single_owner = 1 if nb_owners == 1
		gen multiple_owner = 1 if nb_owners > 1
		gen two_owners = 1 if nb_owners == 2
		gen three_owners = 1 if nb_owners == 3
		gen four_owners = 1 if nb_owners == 4
		gen five_owners = 1 if nb_owners == 5
		gen two_four_owners = 1 if nb_owners > 1 & nb_owners < 5
		gen five_more_owners = 1 if nb_owners > 4
		gen six_more_owners = 1 if nb_owners > 5
		
*Age of firms
* Identify the first year of each ID (birth year)
bysort UCID_1934_2003 (year): gen birth_year = year if _n == 1
bysort UCID_1934_2003 (year): replace birth_year = birth_year[_n-1] if missing(birth_year)
gen firm_age = (year - birth_year) + 1

********************************************************************************
* Split in cross-sections to win time later 
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp{
	preserve
		keep if year == `a'
		save "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`a'_movements.dta", replace
	restore
	drop if year == `a'
}

********************************************************************************
* gdenr
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'_movements.dta", clear

	*Stock of firms
	gen tag = 1
	egen ndistinct = total(tag), by(gdenr)
	drop tag
	
	*Owners
		egen total_owners = total(nb_owners), by(gdenr)
		gen avg_owners = total_owners/ndistinct
		egen min_owners = min(nb_owners), by(gdenr)
		egen max_owners = max(nb_owners), by(gdenr)
		gen log_avg_owners = ln(avg_owners)
		
		summarize nb_owners, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_o = r(p`p') 
		}
		
		gen top50_dummy_o = 1 if nb_owners>p50_thresh_o
		gen top75_dummy_o = 1 if nb_owners>p75_thresh_o
		gen top90_dummy_o = 1 if nb_owners>p90_thresh_o
		gen top99_dummy_o = 1 if nb_owners>p99_thresh_o
		
		egen total_nb_top50_o = sum(top50_dummy), by(gdenr)
		egen total_nb_top75_o = sum(top75_dummy), by(gdenr)
		egen total_nb_top90_o = sum(top90_dummy), by(gdenr)
		egen total_nb_top99_o = sum(top99_dummy), by(gdenr)
		
		gen total_nb_top50_to_75_o = total_nb_top50_o - total_nb_top75_o
		gen total_nb_top75_to_90_o = total_nb_top75_o - total_nb_top90_o
		gen total_nb_top90_to_99_o = total_nb_top90_o - total_nb_top99_o
		gen total_nb_top99_to_100_o = total_nb_top99_o
		
		collapse (mean) ndistinct total_owners avg_owners min_owners max_owners log_avg_owners p1_thresh_o p5_thresh_o p10_thresh_o p25_thresh_o p50_thresh_o p75_thresh_o p90_thresh_o p99_thresh_o top50_dummy_o top75_dummy_o top90_dummy_o top99_dummy_o total_nb_top50_o total_nb_top75_o total_nb_top90_o total_nb_top99_o total_nb_top50_to_75_o total_nb_top75_to_90_o total_nb_top90_to_99_o total_nb_top99_to_100_o firm_age (sum) single_owner multiple_owner two_owners three_owners four_owners five_owners two_four_owners five_more_owners six_more_owners occurences relocation_in_2012 relocation_in_2018 relocation_out_2012 relocation_out_2018 creation_2012 creation_2018 max_occurences death_2012 death_2018, by(gdenr)
		gen year = `y'
	save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_final.dta", replace
}

* Append all years (gdenr)
use "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_count_mun_final.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_final.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_2003_count_owners_mun_final.dta", replace

********************************************************************************
* gdenr_2012
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'_movements.dta", clear

	*Stock of firms
	gen tag = 1
	egen ndistinct = total(tag), by(gdenr_2012)
	drop tag
	
	*Owners
		egen total_owners = total(nb_owners), by(gdenr_2012)
		gen avg_owners = total_owners/ndistinct
		egen min_owners = min(nb_owners), by(gdenr_2012)
		egen max_owners = max(nb_owners), by(gdenr_2012)
		gen log_avg_owners = ln(avg_owners)
		
		summarize nb_owners, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_o = r(p`p') 
		}
		
		gen top50_dummy_o = 1 if nb_owners>p50_thresh_o
		gen top75_dummy_o = 1 if nb_owners>p75_thresh_o
		gen top90_dummy_o = 1 if nb_owners>p90_thresh_o
		gen top99_dummy_o = 1 if nb_owners>p99_thresh_o
		
		egen total_nb_top50_o = sum(top50_dummy), by(gdenr_2012)
		egen total_nb_top75_o = sum(top75_dummy), by(gdenr_2012)
		egen total_nb_top90_o = sum(top90_dummy), by(gdenr_2012)
		egen total_nb_top99_o = sum(top99_dummy), by(gdenr_2012)
		
		gen total_nb_top50_to_75_o = total_nb_top50_o - total_nb_top75_o
		gen total_nb_top75_to_90_o = total_nb_top75_o - total_nb_top90_o
		gen total_nb_top90_to_99_o = total_nb_top90_o - total_nb_top99_o
		gen total_nb_top99_to_100_o = total_nb_top99_o
		
		collapse (mean) ndistinct total_owners avg_owners min_owners max_owners log_avg_owners p1_thresh_o p5_thresh_o p10_thresh_o p25_thresh_o p50_thresh_o p75_thresh_o p90_thresh_o p99_thresh_o top50_dummy_o top75_dummy_o top90_dummy_o top99_dummy_o total_nb_top50_o total_nb_top75_o total_nb_top90_o total_nb_top99_o total_nb_top50_to_75_o total_nb_top75_to_90_o total_nb_top90_to_99_o total_nb_top99_to_100_o firm_age (sum) single_owner multiple_owner two_owners three_owners four_owners five_owners two_four_owners five_more_owners six_more_owners occurences relocation_in_2012 relocation_in_2018 relocation_out_2012 relocation_out_2018 creation_2012 creation_2018 max_occurences death_2012 death_2018, by(gdenr_2012)
		gen year = `y'
	save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_final_12.dta", replace
}

* Append all years (gdenr_2012)
use "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_count_mun_final_12.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_final_12.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_2003_count_owners_mun_12.dta", replace

********************************************************************************
* gdenr_2018
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'_movements.dta", clear

	*Stock of firms
	gen tag = 1
	egen ndistinct = total(tag), by(gdenr_2018)
	drop tag
	
	*Owners
		egen total_owners = total(nb_owners), by(gdenr_2018)
		gen avg_owners = total_owners/ndistinct
		egen min_owners = min(nb_owners), by(gdenr_2018)
		egen max_owners = max(nb_owners), by(gdenr_2018)
		gen log_avg_owners = ln(avg_owners)
		
		summarize nb_owners, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_o = r(p`p') 
		}
		
		gen top50_dummy_o = 1 if nb_owners>p50_thresh_o
		gen top75_dummy_o = 1 if nb_owners>p75_thresh_o
		gen top90_dummy_o = 1 if nb_owners>p90_thresh_o
		gen top99_dummy_o = 1 if nb_owners>p99_thresh_o
		
		egen total_nb_top50_o = sum(top50_dummy), by(gdenr_2018)
		egen total_nb_top75_o = sum(top75_dummy), by(gdenr_2018)
		egen total_nb_top90_o = sum(top90_dummy), by(gdenr_2018)
		egen total_nb_top99_o = sum(top99_dummy), by(gdenr_2018)
		
		gen total_nb_top50_to_75_o = total_nb_top50_o - total_nb_top75_o
		gen total_nb_top75_to_90_o = total_nb_top75_o - total_nb_top90_o
		gen total_nb_top90_to_99_o = total_nb_top90_o - total_nb_top99_o
		gen total_nb_top99_to_100_o = total_nb_top99_o
		
		collapse (mean) ndistinct total_owners avg_owners min_owners max_owners log_avg_owners p1_thresh_o p5_thresh_o p10_thresh_o p25_thresh_o p50_thresh_o p75_thresh_o p90_thresh_o p99_thresh_o top50_dummy_o top75_dummy_o top90_dummy_o top99_dummy_o total_nb_top50_o total_nb_top75_o total_nb_top90_o total_nb_top99_o total_nb_top50_to_75_o total_nb_top75_to_90_o total_nb_top90_to_99_o total_nb_top99_to_100_o firm_age (sum) single_owner multiple_owner two_owners three_owners four_owners five_owners two_four_owners five_more_owners six_more_owners occurences relocation_in_2012 relocation_in_2018 relocation_out_2012 relocation_out_2018 creation_2012 creation_2018 max_occurences death_2012 death_2018, by(gdenr_2018)
		gen year = `y'
	save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_final_18.dta", replace
}

* Append all years (gdenr_2018)
use "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_count_mun_final_18.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_final_18.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_2003_count_owners_mun_18.dta", replace

********************************************************************************
* Gen data on nominal capital

*** Corrections of capital (gdenr)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {

use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'_movements.dta", clear

rename capital capital_orig
gen capital = capital_orig
replace capital = ustrlower(capital)
replace capital = ustrregexra(capital,"[a-z]","")
replace capital = usubinstr(capital, ")", "",.)
replace capital = usubinstr(capital, "í", "",.)
replace capital = usubinstr(capital, "'", "",.)
replace capital = usubinstr(capital, "-", "",.)
replace capital = usubinstr(capital, ",", "",.)
replace capital = usubinstr(capital, "!", "",.)
replace capital = usubinstr(capital, "`", "",.)
replace capital = usubinstr(capital, "�", "",.)
replace capital = usubinstr(capital, "/", "",.)
replace capital = usubinstr(capital, "?", "",.)
replace capital = usubinstr(capital, ":", "",.)
replace capital = usubinstr(capital, "è", "",.)
replace capital = usubinstr(capital, "ó", "",.)
replace capital = usubinstr(capital, "ö", "",.)
replace capital = usubinstr(capital, "ü", "",.)
replace capital = usubinstr(capital, "ä", "",.)
replace capital = usubinstr(capital, "   ", " ",.)
replace capital = usubinstr(capital, "  ", " ",.)
replace capital = usubinstr(capital, " ", "",.)
replace capital = ustrtrim(capital)
destring capital, replace force 
replace capital = 0.5 if capital_orig == "0.S"
replace capital = 0.5 if capital_orig == "0.5)."
replace capital = 0.05 if capital_orig == "0:05"
replace capital = 0.3 if capital_orig == "0:3"
replace capital = 0.3 if capital_orig == "0.3) 5"
replace capital = 0.1 if capital_orig == "0-1"
replace capital = 0.1 if capital_orig == "0.'"
replace capital = 0.051 if capital_orig == "0.05l"
replace capital = 32.66 if capital_orig == "32.6602).5"
replace capital = 0.1 if capital_orig == "/0,1"
replace capital = 145.0 if capital_orig == "145:0"
replace capital = 0.1 if capital_orig == "0.1) 5"
replace capital = 0.98 if capital_orig == "O,98"
replace capital = 0.06 if capital_orig == "0.C6"
replace capital = . if capital_orig == "Liesta140,1"
replace capital = . if capital_orig == ",4,0"
replace capital = . if capital_orig == "41,"
replace capital = 0.1 if capital_orig == "0.i"
replace capital = 0.05 if capital_orig == "0..05"
replace capital = 0.05 if capital_orig == "0.0.5"
replace capital = 0.05 if capital_orig == "0.05) Pr.-"
replace capital = 0.1 if capital_orig == "0.1) Vr."
replace capital = 0.1 if capital_orig == "0.1."
replace capital = 0.15 if capital_orig == "0.15."
replace capital = 0.1 if capital_orig == "0.25)' t..."

drop if missing(capital) | capital == .
		
		gen tag = 1
		egen ndistinct = total(tag), by(gdenr)
		drop tag
		egen total_capital = sum(capital), by(gdenr)
		gen avg_capital = total_capital/ndistinct
		egen min_capital = min(capital), by(gdenr)
		egen max_capital = max(capital), by(gdenr)
		gen log_avg_capital = ln(avg_capital)
		
		summarize capital, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_c = r(p`p') 
		}
		
		gen top50_dummy_c = 1 if capital>p50_thresh_c
		gen top75_dummy_c = 1 if capital>p75_thresh_c
		gen top90_dummy_c = 1 if capital>p90_thresh_c
		gen top99_dummy_c = 1 if capital>p99_thresh_c
		*gen bottom50_dummy_c = 1 if capital<=p50_thresh_c
		
		egen total_nb_top50_c = sum(top50_dummy_c), by(gdenr)
		egen total_nb_top75_c = sum(top75_dummy_c), by(gdenr)
		egen total_nb_top90_c = sum(top90_dummy_c), by(gdenr)
		egen total_nb_top99_c = sum(top99_dummy_c), by(gdenr)
		
		gen total_nb_top50_to_75_c = total_nb_top50_c - total_nb_top75_c
		gen total_nb_top75_to_90_c = total_nb_top75_c - total_nb_top90_c
		gen total_nb_top90_to_99_c = total_nb_top90_c - total_nb_top99_c
		gen total_nb_top99_to_100_c = total_nb_top99_c
		
		collapse (mean) ndistinct total_capital avg_capital min_capital max_capital log_avg_capital p1_thresh_c p5_thresh_c p10_thresh_c p25_thresh_c p50_thresh_c p75_thresh_c p90_thresh_c p99_thresh_c top50_dummy_c top75_dummy_c top90_dummy_c top99_dummy_c total_nb_top50_c total_nb_top75_c total_nb_top90_c total_nb_top99_c total_nb_top50_to_75_c total_nb_top75_to_90_c total_nb_top90_to_99_c total_nb_top99_to_100_c, by(gdenr)
		gen year = `y'

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun.dta", replace
}

* Append all years (gdenr)
use "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_capital_wstat_mun.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_2003_capital_wstat_mun.dta", replace

*** Corrections of capital (gdenr_2012)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {

use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'_movements.dta", clear

rename capital capital_orig
gen capital = capital_orig
replace capital = ustrlower(capital)
replace capital = ustrregexra(capital,"[a-z]","")
replace capital = usubinstr(capital, ")", "",.)
replace capital = usubinstr(capital, "í", "",.)
replace capital = usubinstr(capital, "'", "",.)
replace capital = usubinstr(capital, "-", "",.)
replace capital = usubinstr(capital, ",", "",.)
replace capital = usubinstr(capital, "!", "",.)
replace capital = usubinstr(capital, "`", "",.)
replace capital = usubinstr(capital, "�", "",.)
replace capital = usubinstr(capital, "/", "",.)
replace capital = usubinstr(capital, "?", "",.)
replace capital = usubinstr(capital, ":", "",.)
replace capital = usubinstr(capital, "è", "",.)
replace capital = usubinstr(capital, "ó", "",.)
replace capital = usubinstr(capital, "ö", "",.)
replace capital = usubinstr(capital, "ü", "",.)
replace capital = usubinstr(capital, "ä", "",.)
replace capital = usubinstr(capital, "   ", " ",.)
replace capital = usubinstr(capital, "  ", " ",.)
replace capital = usubinstr(capital, " ", "",.)
replace capital = ustrtrim(capital)
destring capital, replace force 
replace capital = 0.5 if capital_orig == "0.S"
replace capital = 0.5 if capital_orig == "0.5)."
replace capital = 0.05 if capital_orig == "0:05"
replace capital = 0.3 if capital_orig == "0:3"
replace capital = 0.3 if capital_orig == "0.3) 5"
replace capital = 0.1 if capital_orig == "0-1"
replace capital = 0.1 if capital_orig == "0.'"
replace capital = 0.051 if capital_orig == "0.05l"
replace capital = 32.66 if capital_orig == "32.6602).5"
replace capital = 0.1 if capital_orig == "/0,1"
replace capital = 145.0 if capital_orig == "145:0"
replace capital = 0.1 if capital_orig == "0.1) 5"
replace capital = 0.98 if capital_orig == "O,98"
replace capital = 0.06 if capital_orig == "0.C6"
replace capital = . if capital_orig == "Liesta140,1"
replace capital = . if capital_orig == ",4,0"
replace capital = . if capital_orig == "41,"
replace capital = 0.1 if capital_orig == "0.i"
replace capital = 0.05 if capital_orig == "0..05"
replace capital = 0.05 if capital_orig == "0.0.5"
replace capital = 0.05 if capital_orig == "0.05) Pr.-"
replace capital = 0.1 if capital_orig == "0.1) Vr."
replace capital = 0.1 if capital_orig == "0.1."
replace capital = 0.15 if capital_orig == "0.15."
replace capital = 0.1 if capital_orig == "0.25)' t..."

drop if missing(capital) | capital == .
		
		gen tag = 1
		egen ndistinct = total(tag), by(gdenr_2012)
		drop tag
		egen total_capital = sum(capital), by(gdenr_2012)
		gen avg_capital = total_capital/ndistinct
		egen min_capital = min(capital), by(gdenr_2012)
		egen max_capital = max(capital), by(gdenr_2012)
		gen log_avg_capital = ln(avg_capital)
		
		summarize capital, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_c = r(p`p') 
		}
		
		gen top50_dummy_c = 1 if capital>p50_thresh_c
		gen top75_dummy_c = 1 if capital>p75_thresh_c
		gen top90_dummy_c = 1 if capital>p90_thresh_c
		gen top99_dummy_c = 1 if capital>p99_thresh_c
		*gen bottom50_dummy_c = 1 if capital<=p50_thresh_c
		
		egen total_nb_top50_c = sum(top50_dummy_c), by(gdenr_2012)
		egen total_nb_top75_c = sum(top75_dummy_c), by(gdenr_2012)
		egen total_nb_top90_c = sum(top90_dummy_c), by(gdenr_2012)
		egen total_nb_top99_c = sum(top99_dummy_c), by(gdenr_2012)
		
		gen total_nb_top50_to_75_c = total_nb_top50_c - total_nb_top75_c
		gen total_nb_top75_to_90_c = total_nb_top75_c - total_nb_top90_c
		gen total_nb_top90_to_99_c = total_nb_top90_c - total_nb_top99_c
		gen total_nb_top99_to_100_c = total_nb_top99_c
		
		collapse (mean) ndistinct total_capital avg_capital min_capital max_capital log_avg_capital p1_thresh_c p5_thresh_c p10_thresh_c p25_thresh_c p50_thresh_c p75_thresh_c p90_thresh_c p99_thresh_c top50_dummy_c top75_dummy_c top90_dummy_c top99_dummy_c total_nb_top50_c total_nb_top75_c total_nb_top90_c total_nb_top99_c total_nb_top50_to_75_c total_nb_top75_to_90_c total_nb_top90_to_99_c total_nb_top99_to_100_c, by(gdenr_2012)
		gen year = `y'

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_12.dta", replace
}

* Append all years (gdenr_2012)
use "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_capital_wstat_mun_12.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_12.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_2003_capital_wstat_mun_12.dta", replace

*** Corrections of capital (gdenr_2018)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {

use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'_movements.dta", clear

rename capital capital_orig
gen capital = capital_orig
replace capital = ustrlower(capital)
replace capital = ustrregexra(capital,"[a-z]","")
replace capital = usubinstr(capital, ")", "",.)
replace capital = usubinstr(capital, "í", "",.)
replace capital = usubinstr(capital, "'", "",.)
replace capital = usubinstr(capital, "-", "",.)
replace capital = usubinstr(capital, ",", "",.)
replace capital = usubinstr(capital, "!", "",.)
replace capital = usubinstr(capital, "`", "",.)
replace capital = usubinstr(capital, "�", "",.)
replace capital = usubinstr(capital, "/", "",.)
replace capital = usubinstr(capital, "?", "",.)
replace capital = usubinstr(capital, ":", "",.)
replace capital = usubinstr(capital, "è", "",.)
replace capital = usubinstr(capital, "ó", "",.)
replace capital = usubinstr(capital, "ö", "",.)
replace capital = usubinstr(capital, "ü", "",.)
replace capital = usubinstr(capital, "ä", "",.)
replace capital = usubinstr(capital, "   ", " ",.)
replace capital = usubinstr(capital, "  ", " ",.)
replace capital = usubinstr(capital, " ", "",.)
replace capital = ustrtrim(capital)
destring capital, replace force 
replace capital = 0.5 if capital_orig == "0.S"
replace capital = 0.5 if capital_orig == "0.5)."
replace capital = 0.05 if capital_orig == "0:05"
replace capital = 0.3 if capital_orig == "0:3"
replace capital = 0.3 if capital_orig == "0.3) 5"
replace capital = 0.1 if capital_orig == "0-1"
replace capital = 0.1 if capital_orig == "0.'"
replace capital = 0.051 if capital_orig == "0.05l"
replace capital = 32.66 if capital_orig == "32.6602).5"
replace capital = 0.1 if capital_orig == "/0,1"
replace capital = 145.0 if capital_orig == "145:0"
replace capital = 0.1 if capital_orig == "0.1) 5"
replace capital = 0.98 if capital_orig == "O,98"
replace capital = 0.06 if capital_orig == "0.C6"
replace capital = . if capital_orig == "Liesta140,1"
replace capital = . if capital_orig == ",4,0"
replace capital = . if capital_orig == "41,"
replace capital = 0.1 if capital_orig == "0.i"
replace capital = 0.05 if capital_orig == "0..05"
replace capital = 0.05 if capital_orig == "0.0.5"
replace capital = 0.05 if capital_orig == "0.05) Pr.-"
replace capital = 0.1 if capital_orig == "0.1) Vr."
replace capital = 0.1 if capital_orig == "0.1."
replace capital = 0.15 if capital_orig == "0.15."
replace capital = 0.1 if capital_orig == "0.25)' t..."

drop if missing(capital) | capital == .
		
		gen tag = 1
		egen ndistinct = total(tag), by(gdenr_2018)
		drop tag
		egen total_capital = sum(capital), by(gdenr_2018)
		gen avg_capital = total_capital/ndistinct
		egen min_capital = min(capital), by(gdenr_2018)
		egen max_capital = max(capital), by(gdenr_2018)
		gen log_avg_capital = ln(avg_capital)
		
		summarize capital, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_c = r(p`p') 
		}
		
		gen top50_dummy_c = 1 if capital>p50_thresh_c
		gen top75_dummy_c = 1 if capital>p75_thresh_c
		gen top90_dummy_c = 1 if capital>p90_thresh_c
		gen top99_dummy_c = 1 if capital>p99_thresh_c
		*gen bottom50_dummy_c = 1 if capital<=p50_thresh_c
		
		egen total_nb_top50_c = sum(top50_dummy_c), by(gdenr_2018)
		egen total_nb_top75_c = sum(top75_dummy_c), by(gdenr_2018)
		egen total_nb_top90_c = sum(top90_dummy_c), by(gdenr_2018)
		egen total_nb_top99_c = sum(top99_dummy_c), by(gdenr_2018)
		
		gen total_nb_top50_to_75_c = total_nb_top50_c - total_nb_top75_c
		gen total_nb_top75_to_90_c = total_nb_top75_c - total_nb_top90_c
		gen total_nb_top90_to_99_c = total_nb_top90_c - total_nb_top99_c
		gen total_nb_top99_to_100_c = total_nb_top99_c
		
		collapse (mean) ndistinct total_capital avg_capital min_capital max_capital log_avg_capital p1_thresh_c p5_thresh_c p10_thresh_c p25_thresh_c p50_thresh_c p75_thresh_c p90_thresh_c p99_thresh_c top50_dummy_c top75_dummy_c top90_dummy_c top99_dummy_c total_nb_top50_c total_nb_top75_c total_nb_top90_c total_nb_top99_c total_nb_top50_to_75_c total_nb_top75_to_90_c total_nb_top90_to_99_c total_nb_top99_to_100_c, by(gdenr_2018)
		gen year = `y'

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_18.dta", replace
}

* Append all years (gdenr_2018)
use "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_capital_wstat_mun_18.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_18.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_2003_capital_wstat_mun_18.dta", replace



********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************



* Gen same data w/o single obs in panel (obs outside series)

*Priority now !!!

********************************************************************************
* gdenr
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'_movements.dta", clear
	drop if max_occurences == 0
	
	*Stock of firms
	gen tag = 1
	egen ndistinct = total(tag), by(gdenr)
	drop tag
	
	*Owners
		egen total_owners = total(nb_owners), by(gdenr)
		gen avg_owners = total_owners/ndistinct
		egen min_owners = min(nb_owners), by(gdenr)
		egen max_owners = max(nb_owners), by(gdenr)
		gen log_avg_owners = ln(avg_owners)
		
		summarize nb_owners, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_o = r(p`p') 
		}
		
		gen top50_dummy_o = 1 if nb_owners>p50_thresh_o
		gen top75_dummy_o = 1 if nb_owners>p75_thresh_o
		gen top90_dummy_o = 1 if nb_owners>p90_thresh_o
		gen top99_dummy_o = 1 if nb_owners>p99_thresh_o
		
		egen total_nb_top50_o = sum(top50_dummy), by(gdenr)
		egen total_nb_top75_o = sum(top75_dummy), by(gdenr)
		egen total_nb_top90_o = sum(top90_dummy), by(gdenr)
		egen total_nb_top99_o = sum(top99_dummy), by(gdenr)
		
		gen total_nb_top50_to_75_o = total_nb_top50_o - total_nb_top75_o
		gen total_nb_top75_to_90_o = total_nb_top75_o - total_nb_top90_o
		gen total_nb_top90_to_99_o = total_nb_top90_o - total_nb_top99_o
		gen total_nb_top99_to_100_o = total_nb_top99_o
		
		collapse (mean) ndistinct total_owners avg_owners min_owners max_owners log_avg_owners p1_thresh_o p5_thresh_o p10_thresh_o p25_thresh_o p50_thresh_o p75_thresh_o p90_thresh_o p99_thresh_o top50_dummy_o top75_dummy_o top90_dummy_o top99_dummy_o total_nb_top50_o total_nb_top75_o total_nb_top90_o total_nb_top99_o total_nb_top50_to_75_o total_nb_top75_to_90_o total_nb_top90_to_99_o total_nb_top99_to_100_o firm_age (sum) single_owner multiple_owner two_owners three_owners four_owners five_owners two_four_owners five_more_owners six_more_owners occurences relocation_in_2012 relocation_in_2018 relocation_out_2012 relocation_out_2018 creation_2012 creation_2018 max_occurences death_2012 death_2018, by(gdenr)
		gen year = `y'
	save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_final_nosingle.dta", replace
}

* Append all years (gdenr)
use "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_count_mun_final_nosingle.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_final_nosingle.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_2003_count_owners_mun_final_nosingle.dta", replace

********************************************************************************
* gdenr_2012
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'_movements.dta", clear
	drop if max_occurences == 0
	
	*Stock of firms
	gen tag = 1
	egen ndistinct = total(tag), by(gdenr_2012)
	drop tag
	
	*Owners
		egen total_owners = total(nb_owners), by(gdenr_2012)
		gen avg_owners = total_owners/ndistinct
		egen min_owners = min(nb_owners), by(gdenr_2012)
		egen max_owners = max(nb_owners), by(gdenr_2012)
		gen log_avg_owners = ln(avg_owners)
		
		summarize nb_owners, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_o = r(p`p') 
		}
		
		gen top50_dummy_o = 1 if nb_owners>p50_thresh_o
		gen top75_dummy_o = 1 if nb_owners>p75_thresh_o
		gen top90_dummy_o = 1 if nb_owners>p90_thresh_o
		gen top99_dummy_o = 1 if nb_owners>p99_thresh_o
		
		egen total_nb_top50_o = sum(top50_dummy), by(gdenr_2012)
		egen total_nb_top75_o = sum(top75_dummy), by(gdenr_2012)
		egen total_nb_top90_o = sum(top90_dummy), by(gdenr_2012)
		egen total_nb_top99_o = sum(top99_dummy), by(gdenr_2012)
		
		gen total_nb_top50_to_75_o = total_nb_top50_o - total_nb_top75_o
		gen total_nb_top75_to_90_o = total_nb_top75_o - total_nb_top90_o
		gen total_nb_top90_to_99_o = total_nb_top90_o - total_nb_top99_o
		gen total_nb_top99_to_100_o = total_nb_top99_o
		
		collapse (mean) ndistinct total_owners avg_owners min_owners max_owners log_avg_owners p1_thresh_o p5_thresh_o p10_thresh_o p25_thresh_o p50_thresh_o p75_thresh_o p90_thresh_o p99_thresh_o top50_dummy_o top75_dummy_o top90_dummy_o top99_dummy_o total_nb_top50_o total_nb_top75_o total_nb_top90_o total_nb_top99_o total_nb_top50_to_75_o total_nb_top75_to_90_o total_nb_top90_to_99_o total_nb_top99_to_100_o firm_age (sum) single_owner multiple_owner two_owners three_owners four_owners five_owners two_four_owners five_more_owners six_more_owners occurences relocation_in_2012 relocation_in_2018 relocation_out_2012 relocation_out_2018 creation_2012 creation_2018 max_occurences death_2012 death_2018, by(gdenr_2012)
		gen year = `y'
	save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_final_12_nosingle.dta", replace
}

* Append all years (gdenr_2012)
use "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_count_mun_final_12_nosingle.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_final_12_nosingle.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_2003_count_owners_mun_12_nosingle.dta", replace

********************************************************************************
* gdenr_2018
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'_movements.dta", clear
	drop if max_occurences == 0
	
	*Stock of firms
	gen tag = 1
	egen ndistinct = total(tag), by(gdenr_2018)
	drop tag
	
	*Owners
		egen total_owners = total(nb_owners), by(gdenr_2018)
		gen avg_owners = total_owners/ndistinct
		egen min_owners = min(nb_owners), by(gdenr_2018)
		egen max_owners = max(nb_owners), by(gdenr_2018)
		gen log_avg_owners = ln(avg_owners)
		
		summarize nb_owners, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_o = r(p`p') 
		}
		
		gen top50_dummy_o = 1 if nb_owners>p50_thresh_o
		gen top75_dummy_o = 1 if nb_owners>p75_thresh_o
		gen top90_dummy_o = 1 if nb_owners>p90_thresh_o
		gen top99_dummy_o = 1 if nb_owners>p99_thresh_o
		
		egen total_nb_top50_o = sum(top50_dummy), by(gdenr_2018)
		egen total_nb_top75_o = sum(top75_dummy), by(gdenr_2018)
		egen total_nb_top90_o = sum(top90_dummy), by(gdenr_2018)
		egen total_nb_top99_o = sum(top99_dummy), by(gdenr_2018)
		
		gen total_nb_top50_to_75_o = total_nb_top50_o - total_nb_top75_o
		gen total_nb_top75_to_90_o = total_nb_top75_o - total_nb_top90_o
		gen total_nb_top90_to_99_o = total_nb_top90_o - total_nb_top99_o
		gen total_nb_top99_to_100_o = total_nb_top99_o
		
		collapse (mean) ndistinct total_owners avg_owners min_owners max_owners log_avg_owners p1_thresh_o p5_thresh_o p10_thresh_o p25_thresh_o p50_thresh_o p75_thresh_o p90_thresh_o p99_thresh_o top50_dummy_o top75_dummy_o top90_dummy_o top99_dummy_o total_nb_top50_o total_nb_top75_o total_nb_top90_o total_nb_top99_o total_nb_top50_to_75_o total_nb_top75_to_90_o total_nb_top90_to_99_o total_nb_top99_to_100_o firm_age (sum) single_owner multiple_owner two_owners three_owners four_owners five_owners two_four_owners five_more_owners six_more_owners occurences relocation_in_2012 relocation_in_2018 relocation_out_2012 relocation_out_2018 creation_2012 creation_2018 max_occurences death_2012 death_2018, by(gdenr_2018)
		gen year = `y'
	save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_final_18_nosingle.dta", replace
}

* Append all years (gdenr_2018)
use "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_count_mun_final_18_nosingle.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_final_18_nosingle.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_2003_count_owners_mun_18_nosingle.dta", replace

********************************************************************************
* Gen data on nominal capital

*** Corrections of capital (gdenr)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {

use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'_movements.dta", clear
drop if max_occurences == 0

rename capital capital_orig
gen capital = capital_orig
replace capital = ustrlower(capital)
replace capital = ustrregexra(capital,"[a-z]","")
replace capital = usubinstr(capital, ")", "",.)
replace capital = usubinstr(capital, "í", "",.)
replace capital = usubinstr(capital, "'", "",.)
replace capital = usubinstr(capital, "-", "",.)
replace capital = usubinstr(capital, ",", "",.)
replace capital = usubinstr(capital, "!", "",.)
replace capital = usubinstr(capital, "`", "",.)
replace capital = usubinstr(capital, "�", "",.)
replace capital = usubinstr(capital, "/", "",.)
replace capital = usubinstr(capital, "?", "",.)
replace capital = usubinstr(capital, ":", "",.)
replace capital = usubinstr(capital, "è", "",.)
replace capital = usubinstr(capital, "ó", "",.)
replace capital = usubinstr(capital, "ö", "",.)
replace capital = usubinstr(capital, "ü", "",.)
replace capital = usubinstr(capital, "ä", "",.)
replace capital = usubinstr(capital, "   ", " ",.)
replace capital = usubinstr(capital, "  ", " ",.)
replace capital = usubinstr(capital, " ", "",.)
replace capital = ustrtrim(capital)
destring capital, replace force 
replace capital = 0.5 if capital_orig == "0.S"
replace capital = 0.5 if capital_orig == "0.5)."
replace capital = 0.05 if capital_orig == "0:05"
replace capital = 0.3 if capital_orig == "0:3"
replace capital = 0.3 if capital_orig == "0.3) 5"
replace capital = 0.1 if capital_orig == "0-1"
replace capital = 0.1 if capital_orig == "0.'"
replace capital = 0.051 if capital_orig == "0.05l"
replace capital = 32.66 if capital_orig == "32.6602).5"
replace capital = 0.1 if capital_orig == "/0,1"
replace capital = 145.0 if capital_orig == "145:0"
replace capital = 0.1 if capital_orig == "0.1) 5"
replace capital = 0.98 if capital_orig == "O,98"
replace capital = 0.06 if capital_orig == "0.C6"
replace capital = . if capital_orig == "Liesta140,1"
replace capital = . if capital_orig == ",4,0"
replace capital = . if capital_orig == "41,"
replace capital = 0.1 if capital_orig == "0.i"
replace capital = 0.05 if capital_orig == "0..05"
replace capital = 0.05 if capital_orig == "0.0.5"
replace capital = 0.05 if capital_orig == "0.05) Pr.-"
replace capital = 0.1 if capital_orig == "0.1) Vr."
replace capital = 0.1 if capital_orig == "0.1."
replace capital = 0.15 if capital_orig == "0.15."
replace capital = 0.1 if capital_orig == "0.25)' t..."

drop if missing(capital) | capital == .
		
		gen tag = 1
		egen ndistinct = total(tag), by(gdenr)
		drop tag
		egen total_capital = sum(capital), by(gdenr)
		gen avg_capital = total_capital/ndistinct
		egen min_capital = min(capital), by(gdenr)
		egen max_capital = max(capital), by(gdenr)
		gen log_avg_capital = ln(avg_capital)
		
		summarize capital, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_c = r(p`p') 
		}
		
		gen top50_dummy_c = 1 if capital>p50_thresh_c
		gen top75_dummy_c = 1 if capital>p75_thresh_c
		gen top90_dummy_c = 1 if capital>p90_thresh_c
		gen top99_dummy_c = 1 if capital>p99_thresh_c
		*gen bottom50_dummy_c = 1 if capital<=p50_thresh_c
		
		egen total_nb_top50_c = sum(top50_dummy_c), by(gdenr)
		egen total_nb_top75_c = sum(top75_dummy_c), by(gdenr)
		egen total_nb_top90_c = sum(top90_dummy_c), by(gdenr)
		egen total_nb_top99_c = sum(top99_dummy_c), by(gdenr)
		
		gen total_nb_top50_to_75_c = total_nb_top50_c - total_nb_top75_c
		gen total_nb_top75_to_90_c = total_nb_top75_c - total_nb_top90_c
		gen total_nb_top90_to_99_c = total_nb_top90_c - total_nb_top99_c
		gen total_nb_top99_to_100_c = total_nb_top99_c
		
		collapse (mean) ndistinct total_capital avg_capital min_capital max_capital log_avg_capital p1_thresh_c p5_thresh_c p10_thresh_c p25_thresh_c p50_thresh_c p75_thresh_c p90_thresh_c p99_thresh_c top50_dummy_c top75_dummy_c top90_dummy_c top99_dummy_c total_nb_top50_c total_nb_top75_c total_nb_top90_c total_nb_top99_c total_nb_top50_to_75_c total_nb_top75_to_90_c total_nb_top90_to_99_c total_nb_top99_to_100_c, by(gdenr)
		gen year = `y'

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_nosingle.dta", replace
}

* Append all years (gdenr)
use "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_capital_wstat_mun_nosingle.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_nosingle.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_2003_capital_wstat_mun_nosingle.dta", replace

*** Corrections of capital (gdenr_2012)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {

use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'_movements.dta", clear
drop if max_occurences == 0

rename capital capital_orig
gen capital = capital_orig
replace capital = ustrlower(capital)
replace capital = ustrregexra(capital,"[a-z]","")
replace capital = usubinstr(capital, ")", "",.)
replace capital = usubinstr(capital, "í", "",.)
replace capital = usubinstr(capital, "'", "",.)
replace capital = usubinstr(capital, "-", "",.)
replace capital = usubinstr(capital, ",", "",.)
replace capital = usubinstr(capital, "!", "",.)
replace capital = usubinstr(capital, "`", "",.)
replace capital = usubinstr(capital, "�", "",.)
replace capital = usubinstr(capital, "/", "",.)
replace capital = usubinstr(capital, "?", "",.)
replace capital = usubinstr(capital, ":", "",.)
replace capital = usubinstr(capital, "è", "",.)
replace capital = usubinstr(capital, "ó", "",.)
replace capital = usubinstr(capital, "ö", "",.)
replace capital = usubinstr(capital, "ü", "",.)
replace capital = usubinstr(capital, "ä", "",.)
replace capital = usubinstr(capital, "   ", " ",.)
replace capital = usubinstr(capital, "  ", " ",.)
replace capital = usubinstr(capital, " ", "",.)
replace capital = ustrtrim(capital)
destring capital, replace force 
replace capital = 0.5 if capital_orig == "0.S"
replace capital = 0.5 if capital_orig == "0.5)."
replace capital = 0.05 if capital_orig == "0:05"
replace capital = 0.3 if capital_orig == "0:3"
replace capital = 0.3 if capital_orig == "0.3) 5"
replace capital = 0.1 if capital_orig == "0-1"
replace capital = 0.1 if capital_orig == "0.'"
replace capital = 0.051 if capital_orig == "0.05l"
replace capital = 32.66 if capital_orig == "32.6602).5"
replace capital = 0.1 if capital_orig == "/0,1"
replace capital = 145.0 if capital_orig == "145:0"
replace capital = 0.1 if capital_orig == "0.1) 5"
replace capital = 0.98 if capital_orig == "O,98"
replace capital = 0.06 if capital_orig == "0.C6"
replace capital = . if capital_orig == "Liesta140,1"
replace capital = . if capital_orig == ",4,0"
replace capital = . if capital_orig == "41,"
replace capital = 0.1 if capital_orig == "0.i"
replace capital = 0.05 if capital_orig == "0..05"
replace capital = 0.05 if capital_orig == "0.0.5"
replace capital = 0.05 if capital_orig == "0.05) Pr.-"
replace capital = 0.1 if capital_orig == "0.1) Vr."
replace capital = 0.1 if capital_orig == "0.1."
replace capital = 0.15 if capital_orig == "0.15."
replace capital = 0.1 if capital_orig == "0.25)' t..."

drop if missing(capital) | capital == .
		
		gen tag = 1
		egen ndistinct = total(tag), by(gdenr_2012)
		drop tag
		egen total_capital = sum(capital), by(gdenr_2012)
		gen avg_capital = total_capital/ndistinct
		egen min_capital = min(capital), by(gdenr_2012)
		egen max_capital = max(capital), by(gdenr_2012)
		gen log_avg_capital = ln(avg_capital)
		
		summarize capital, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_c = r(p`p') 
		}
		
		gen top50_dummy_c = 1 if capital>p50_thresh_c
		gen top75_dummy_c = 1 if capital>p75_thresh_c
		gen top90_dummy_c = 1 if capital>p90_thresh_c
		gen top99_dummy_c = 1 if capital>p99_thresh_c
		*gen bottom50_dummy_c = 1 if capital<=p50_thresh_c
		
		egen total_nb_top50_c = sum(top50_dummy_c), by(gdenr_2012)
		egen total_nb_top75_c = sum(top75_dummy_c), by(gdenr_2012)
		egen total_nb_top90_c = sum(top90_dummy_c), by(gdenr_2012)
		egen total_nb_top99_c = sum(top99_dummy_c), by(gdenr_2012)
		
		gen total_nb_top50_to_75_c = total_nb_top50_c - total_nb_top75_c
		gen total_nb_top75_to_90_c = total_nb_top75_c - total_nb_top90_c
		gen total_nb_top90_to_99_c = total_nb_top90_c - total_nb_top99_c
		gen total_nb_top99_to_100_c = total_nb_top99_c
		
		collapse (mean) ndistinct total_capital avg_capital min_capital max_capital log_avg_capital p1_thresh_c p5_thresh_c p10_thresh_c p25_thresh_c p50_thresh_c p75_thresh_c p90_thresh_c p99_thresh_c top50_dummy_c top75_dummy_c top90_dummy_c top99_dummy_c total_nb_top50_c total_nb_top75_c total_nb_top90_c total_nb_top99_c total_nb_top50_to_75_c total_nb_top75_to_90_c total_nb_top90_to_99_c total_nb_top99_to_100_c, by(gdenr_2012)
		gen year = `y'

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_12_nosingle.dta", replace
}

* Append all years (gdenr_2012)
use "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_capital_wstat_mun_12_nosingle.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_12_nosingle.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_2003_capital_wstat_mun_12_nosingle.dta", replace

*** Corrections of capital (gdenr_2018)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {

use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'_movements.dta", clear
drop if max_occurences == 0

rename capital capital_orig
gen capital = capital_orig
replace capital = ustrlower(capital)
replace capital = ustrregexra(capital,"[a-z]","")
replace capital = usubinstr(capital, ")", "",.)
replace capital = usubinstr(capital, "í", "",.)
replace capital = usubinstr(capital, "'", "",.)
replace capital = usubinstr(capital, "-", "",.)
replace capital = usubinstr(capital, ",", "",.)
replace capital = usubinstr(capital, "!", "",.)
replace capital = usubinstr(capital, "`", "",.)
replace capital = usubinstr(capital, "�", "",.)
replace capital = usubinstr(capital, "/", "",.)
replace capital = usubinstr(capital, "?", "",.)
replace capital = usubinstr(capital, ":", "",.)
replace capital = usubinstr(capital, "è", "",.)
replace capital = usubinstr(capital, "ó", "",.)
replace capital = usubinstr(capital, "ö", "",.)
replace capital = usubinstr(capital, "ü", "",.)
replace capital = usubinstr(capital, "ä", "",.)
replace capital = usubinstr(capital, "   ", " ",.)
replace capital = usubinstr(capital, "  ", " ",.)
replace capital = usubinstr(capital, " ", "",.)
replace capital = ustrtrim(capital)
destring capital, replace force 
replace capital = 0.5 if capital_orig == "0.S"
replace capital = 0.5 if capital_orig == "0.5)."
replace capital = 0.05 if capital_orig == "0:05"
replace capital = 0.3 if capital_orig == "0:3"
replace capital = 0.3 if capital_orig == "0.3) 5"
replace capital = 0.1 if capital_orig == "0-1"
replace capital = 0.1 if capital_orig == "0.'"
replace capital = 0.051 if capital_orig == "0.05l"
replace capital = 32.66 if capital_orig == "32.6602).5"
replace capital = 0.1 if capital_orig == "/0,1"
replace capital = 145.0 if capital_orig == "145:0"
replace capital = 0.1 if capital_orig == "0.1) 5"
replace capital = 0.98 if capital_orig == "O,98"
replace capital = 0.06 if capital_orig == "0.C6"
replace capital = . if capital_orig == "Liesta140,1"
replace capital = . if capital_orig == ",4,0"
replace capital = . if capital_orig == "41,"
replace capital = 0.1 if capital_orig == "0.i"
replace capital = 0.05 if capital_orig == "0..05"
replace capital = 0.05 if capital_orig == "0.0.5"
replace capital = 0.05 if capital_orig == "0.05) Pr.-"
replace capital = 0.1 if capital_orig == "0.1) Vr."
replace capital = 0.1 if capital_orig == "0.1."
replace capital = 0.15 if capital_orig == "0.15."
replace capital = 0.1 if capital_orig == "0.25)' t..."

drop if missing(capital) | capital == .
		
		gen tag = 1
		egen ndistinct = total(tag), by(gdenr_2018)
		drop tag
		egen total_capital = sum(capital), by(gdenr_2018)
		gen avg_capital = total_capital/ndistinct
		egen min_capital = min(capital), by(gdenr_2018)
		egen max_capital = max(capital), by(gdenr_2018)
		gen log_avg_capital = ln(avg_capital)
		
		summarize capital, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_c = r(p`p') 
		}
		
		gen top50_dummy_c = 1 if capital>p50_thresh_c
		gen top75_dummy_c = 1 if capital>p75_thresh_c
		gen top90_dummy_c = 1 if capital>p90_thresh_c
		gen top99_dummy_c = 1 if capital>p99_thresh_c
		*gen bottom50_dummy_c = 1 if capital<=p50_thresh_c
		
		egen total_nb_top50_c = sum(top50_dummy_c), by(gdenr_2018)
		egen total_nb_top75_c = sum(top75_dummy_c), by(gdenr_2018)
		egen total_nb_top90_c = sum(top90_dummy_c), by(gdenr_2018)
		egen total_nb_top99_c = sum(top99_dummy_c), by(gdenr_2018)
		
		gen total_nb_top50_to_75_c = total_nb_top50_c - total_nb_top75_c
		gen total_nb_top75_to_90_c = total_nb_top75_c - total_nb_top90_c
		gen total_nb_top90_to_99_c = total_nb_top90_c - total_nb_top99_c
		gen total_nb_top99_to_100_c = total_nb_top99_c
		
		collapse (mean) ndistinct total_capital avg_capital min_capital max_capital log_avg_capital p1_thresh_c p5_thresh_c p10_thresh_c p25_thresh_c p50_thresh_c p75_thresh_c p90_thresh_c p99_thresh_c top50_dummy_c top75_dummy_c top90_dummy_c top99_dummy_c total_nb_top50_c total_nb_top75_c total_nb_top90_c total_nb_top99_c total_nb_top50_to_75_c total_nb_top75_to_90_c total_nb_top90_to_99_c total_nb_top99_to_100_c, by(gdenr_2018)
		gen year = `y'

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_18_nosingle.dta", replace
}

* Append all years (gdenr_2018)
use "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_capital_wstat_mun_18_nosingle.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_18_nosingle.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_2003_capital_wstat_mun_18_nosingle.dta", replace
