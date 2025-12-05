clear
set more off
version 18

*** Set paths

global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"


/* Firm names changes

import excel using "$path\01_Raw_data\17_SIX\TOC_historische Kurs Gesellschafts und Instrument Stammdaten.xlsx", sheet("Name Changes") firstrow  clear

keep GKnumber Longnameold Longnamenew 
rename GKnumber GK
rename Longnameold oldname
rename Longnamenew newname
sort GK
save "$path\02_Processed_data\17_SIX\SIX_names_changes.dta", replace

* Firm names today

import excel using "$path\01_Raw_data\17_SIX\TOC_historische Kurs Gesellschafts und Instrument Stammdaten.xlsx", sheet("all price history") firstrow  clear

foreach var of varlist _all {
    tostring `var', replace
	replace `var' = "" if `var' == "."
}

replace P = "2024.09" if P != ""
replace Q = "2024.08" if Q != ""
replace R = "2024.07" if R != ""
replace S = "2024.06" if S != ""
replace T = "2024.05" if T != ""
replace U = "2024.04" if U != ""
replace V = "2024.03" if V != ""
replace W = "2024.02" if W != ""
replace X = "2024.01" if X != ""
replace Y = "2023.12" if Y != ""
replace Z = "2023.11" if Z != ""
replace AA = "2023.10" if AA != ""
replace AB = "2023.09" if AB != ""
replace AC = "2023.08" if AC != ""
replace AD = "2023.07" if AD != ""
replace AE = "2023.06" if AE != ""
replace AF = "2023.05" if AF != ""
replace AG = "2023.04" if AG != ""
replace AH = "2023.03" if AH != ""
replace AI = "2023.02" if AI != ""
replace AJ = "2023.01" if AJ != ""
replace AK = "2022.12" if AK != ""
replace AL = "2022.11" if AL != ""
replace AM = "2022.10" if AM != ""
replace AN = "2022.09" if AN != ""
replace AO = "2022.08" if AO != ""
replace AP = "2022.07" if AP != ""
replace AQ = "2022.06" if AQ != ""
replace AR = "2022.05" if AR != ""
replace AS = "2022.04" if AS != ""
replace AT = "2022.03" if AT != ""
replace AU = "2022.02" if AU != ""
replace AV = "2022.01" if AV != ""
replace AW = "2021.12" if AW != ""
replace AX = "2021.11" if AX != ""
replace AY = "2021.10" if AY != ""
replace AZ = "2021.09" if AZ != ""
replace BA = "2021.08" if BA != ""
replace BB = "2021.07" if BB != ""
replace BC = "2021.06" if BC != ""
replace BD = "2021.05" if BD != ""
replace BE = "2021.04" if BE != ""
replace BF = "2021.03" if BF != ""
replace BG = "2021.02" if BG != ""
replace BH = "2021.01" if BH != ""
replace BI = "2020.12" if BI != ""
replace BJ = "2020.11" if BJ != ""
replace BK = "2020.10" if BK != ""
replace BL = "2020.09" if BL != ""
replace BM = "2020.08" if BM != ""
replace BN = "2020.07" if BN != ""
replace BO = "2020.06" if BO != ""
replace BP = "2020.05" if BP != ""
replace BQ = "2020.04" if BQ != ""
replace BR = "2020.03" if BR != ""
replace BS = "2020.02" if BS != ""
replace BT = "2020.01" if BT != ""
replace BU = "2019.12" if BU != ""
replace BV = "2019.11" if BV != ""
replace BW = "2019.10" if BW != ""
replace BX = "2019.09" if BX != ""
replace BY = "2019.08" if BY != ""
replace BZ = "2019.07" if BZ != ""
replace CA = "2019.06" if CA != ""
replace CB = "2019.05" if CB != ""
replace CC = "2019.04" if CC != ""
replace CD = "2019.03" if CD != ""
replace CE = "2019.02" if CE != ""
replace CF = "2019.01" if CF != ""
replace CG = "2018.12" if CG != ""
replace CH = "2018.11" if CH != ""
replace CI = "2018.10" if CI != ""
replace CJ = "2018.09" if CJ != ""
replace CK = "2018.08" if CK != ""
replace CL = "2018.07" if CL != ""
replace CM = "2018.06" if CM != ""
replace CN = "2018.05" if CN != ""
replace CO = "2018.04" if CO != ""
replace CP = "2018.03" if CP != ""
replace CQ = "2018.02" if CQ != ""
replace CR = "2018.01" if CR != ""
replace CS = "2017.12" if CS != ""
replace CT = "2017.11" if CT != ""
replace CU = "2017.10" if CU != ""
replace CV = "2017.09" if CV != ""
replace CW = "2017.08" if CW != ""
replace CX = "2017.07" if CX != ""
replace CY = "2017.06" if CY != ""
replace CZ = "2017.05" if CZ != ""
replace DA = "2017.04" if DA != ""
replace DB = "2017.03" if DB != ""
replace DC = "2017.02" if DC != ""
replace DD = "2017.01" if DD != ""
replace DE = "2016.12" if DE != ""
replace DF = "2016.11" if DF != ""
replace DG = "2016.10" if DG != ""
replace DH = "2016.09" if DH != ""
replace DI = "2016.08" if DI != ""
replace DJ = "2016.07" if DJ != ""
replace DK = "2016.06" if DK != ""
replace DL = "2016.05" if DL != ""
replace DM = "2016.04" if DM != ""
replace DN = "2016.03" if DN != ""
replace DO = "2016.02" if DO != ""
replace DP = "2016.01" if DP != ""
replace DQ = "2015.12" if DQ != ""
replace DR = "2015.11" if DR != ""
replace DS = "2015.10" if DS != ""
replace DT = "2015.09" if DT != ""
replace DU = "2015.08" if DU != ""
replace DV = "2015.07" if DV != ""
replace DW = "2015.06" if DW != ""
replace DX = "2015.05" if DX != ""
replace DY = "2015.04" if DY != ""
replace DZ = "2015.03" if DZ != ""
replace EA = "2015.02" if EA != ""
replace EB = "2015.01" if EB != ""
replace EC = "2014.12" if EC != ""
replace ED = "2014.11" if ED != ""
replace EE = "2014.10" if EE != ""
replace EF = "2014.09" if EF != ""
replace EG = "2014.08" if EG != ""
replace EH = "2014.07" if EH != ""
replace EI = "2014.06" if EI != ""
replace EJ = "2014.05" if EJ != ""
replace EK = "2014.04" if EK != ""
replace EL = "2014.03" if EL != ""
replace EM = "2014.02" if EM != ""
replace EN = "2014.01" if EN != ""
replace EO = "2013.12" if EO != ""
replace EP = "2013.11" if EP != ""
replace EQ = "2013.10" if EQ != ""
replace ER = "2013.09" if ER != ""
replace ES = "2013.08" if ES != ""
replace ET = "2013.07" if ET != ""
replace EU = "2013.06" if EU != ""
replace EV = "2013.05" if EV != ""
replace EW = "2013.04" if EW != ""
replace EX = "2013.03" if EX != ""
replace EY = "2013.02" if EY != ""
replace EZ = "2013.01" if EZ != ""
replace FA = "2012.12" if FA != ""
replace FB = "2012.11" if FB != ""
replace FC = "2012.10" if FC != ""
replace FD = "2012.09" if FD != ""
replace FE = "2012.08" if FE != ""
replace FF = "2012.07" if FF != ""
replace FG = "2012.06" if FG != ""
replace FH = "2012.05" if FH != ""
replace FI = "2012.04" if FI != ""
replace FJ = "2012.03" if FJ != ""
replace FK = "2012.02" if FK != ""
replace FL = "2012.01" if FL != ""
replace FM = "2011.12" if FM != ""
replace FN = "2011.11" if FN != ""
replace FO = "2011.10" if FO != ""
replace FP = "2011.09" if FP != ""
replace FQ = "2011.08" if FQ != ""
replace FR = "2011.07" if FR != ""
replace FS = "2011.06" if FS != ""
replace FT = "2011.05" if FT != ""
replace FU = "2011.04" if FU != ""
replace FV = "2011.03" if FV != ""
replace FW = "2011.02" if FW != ""
replace FX = "2011.01" if FX != ""
replace FY = "2010.12" if FY != ""
replace FZ = "2010.11" if FZ != ""
replace GA = "2010.10" if GA != ""
replace GB = "2010.09" if GB != ""
replace GC = "2010.08" if GC != ""
replace GD = "2010.07" if GD != ""
replace GE = "2010.06" if GE != ""
replace GF = "2010.05" if GF != ""
replace GG = "2010.04" if GG != ""
replace GH = "2010.03" if GH != ""
replace GI = "2010.02" if GI != ""
replace GJ = "2010.01" if GJ != ""
replace GK = "2009.12" if GK != ""
replace GL = "2009.11" if GL != ""
replace GM = "2009.10" if GM != ""
replace GN = "2009.09" if GN != ""
replace GO = "2009.08" if GO != ""
replace GP = "2009.07" if GP != ""
replace GQ = "2009.06" if GQ != ""
replace GR = "2009.05" if GR != ""
replace GS = "2009.04" if GS != ""
replace GT = "2009.03" if GT != ""
replace GU = "2009.02" if GU != ""
replace GV = "2009.01" if GV != ""
replace GW = "2008.12" if GW != ""
replace GX = "2008.11" if GX != ""
replace GY = "2008.10" if GY != ""
replace GZ = "2008.09" if GZ != ""
replace HA = "2008.08" if HA != ""
replace HB = "2008.07" if HB != ""
replace HC = "2008.06" if HC != ""
replace HD = "2008.05" if HD != ""
replace HE = "2008.04" if HE != ""
replace HF = "2008.03" if HF != ""
replace HG = "2008.02" if HG != ""
replace HH = "2008.01" if HH != ""
replace HI = "2007.12" if HI != ""
replace HJ = "2007.11" if HJ != ""
replace HK = "2007.10" if HK != ""
replace HL = "2007.09" if HL != ""
replace HM = "2007.08" if HM != ""
replace HN = "2007.07" if HN != ""
replace HO = "2007.06" if HO != ""
replace HP = "2007.05" if HP != ""
replace HQ = "2007.04" if HQ != ""
replace HR = "2007.03" if HR != ""
replace HS = "2007.02" if HS != ""
replace HT = "2007.01" if HT != ""
replace HU = "2006.12" if HU != ""
replace HV = "2006.11" if HV != ""
replace HW = "2006.10" if HW != ""
replace HX = "2006.09" if HX != ""
replace HY = "2006.08" if HY != ""
replace HZ = "2006.07" if HZ != ""
replace IA = "2006.06" if IA != ""
replace IB = "2006.05" if IB != ""
replace IC = "2006.04" if IC != ""
replace ID = "2006.03" if ID != ""
replace IE = "2006.02" if IE != ""
replace IF = "2006.01" if IF != ""
replace IG = "2005.12" if IG != ""
replace IH = "2005.11" if IH != ""
replace II = "2005.10" if II != ""
replace IJ = "2005.09" if IJ != ""
replace IK = "2005.08" if IK != ""
replace IL = "2005.07" if IL != ""
replace IM = "2005.06" if IM != ""
replace IN = "2005.05" if IN != ""
replace IO = "2005.04" if IO != ""
replace IP = "2005.03" if IP != ""
replace IQ = "2005.02" if IQ != ""
replace IR = "2005.01" if IR != ""
replace IS = "2004.12" if IS != ""
replace IT = "2004.11" if IT != ""
replace IU = "2004.10" if IU != ""
replace IV = "2004.09" if IV != ""
replace IW = "2004.08" if IW != ""
replace IX = "2004.07" if IX != ""
replace IY = "2004.06" if IY != ""
replace IZ = "2004.05" if IZ != ""
replace JA = "2004.04" if JA != ""
replace JB = "2004.03" if JB != ""
replace JC = "2004.02" if JC != ""
replace JD = "2004.01" if JD != ""
replace JE = "2003.12" if JE != ""
replace JF = "2003.11" if JF != ""
replace JG = "2003.10" if JG != ""
replace JH = "2003.09" if JH != ""
replace JI = "2003.08" if JI != ""
replace JJ = "2003.07" if JJ != ""
replace JK = "2003.06" if JK != ""
replace JL = "2003.05" if JL != ""
replace JM = "2003.04" if JM != ""
replace JN = "2003.03" if JN != ""
replace JO = "2003.02" if JO != ""
replace JP = "2003.01" if JP != ""
replace JQ = "2002.12" if JQ != ""
replace JR = "2002.11" if JR != ""
replace JS = "2002.10" if JS != ""
replace JT = "2002.09" if JT != ""
replace JU = "2002.08" if JU != ""
replace JV = "2002.07" if JV != ""
replace JW = "2002.06" if JW != ""
replace JX = "2002.05" if JX != ""
replace JY = "2002.04" if JY != ""
replace JZ = "2002.03" if JZ != ""
replace KA = "2002.02" if KA != ""
replace KB = "2002.01" if KB != ""
replace KC = "2001.12" if KC != ""
replace KD = "2001.11" if KD != ""
replace KE = "2001.10" if KE != ""
replace KF = "2001.09" if KF != ""
replace KG = "2001.08" if KG != ""
replace KH = "2001.07" if KH != ""
replace KI = "2001.06" if KI != ""
replace KJ = "2001.05" if KJ != ""
replace KK = "2001.04" if KK != ""
replace KL = "2001.03" if KL != ""
replace KM = "2001.02" if KM != ""
replace KN = "2001.01" if KN != ""
replace KO = "2000.12" if KO != ""
replace KP = "2000.11" if KP != ""
replace KQ = "2000.10" if KQ != ""
replace KR = "2000.09" if KR != ""
replace KS = "2000.08" if KS != ""
replace KT = "2000.07" if KT != ""
replace KU = "2000.06" if KU != ""
replace KV = "2000.05" if KV != ""
replace KW = "2000.04" if KW != ""
replace KX = "2000.03" if KX != ""
replace KY = "2000.02" if KY != ""
replace KZ = "2000.01" if KZ != ""
replace LA = "1999.12" if LA != ""
replace LB = "1999.11" if LB != ""
replace LC = "1999.10" if LC != ""
replace LD = "1999.09" if LD != ""
replace LE = "1999.08" if LE != ""
replace LF = "1999.07" if LF != ""
replace LG = "1999.06" if LG != ""
replace LH = "1999.05" if LH != ""
replace LI = "1999.04" if LI != ""
replace LJ = "1999.03" if LJ != ""
replace LK = "1999.02" if LK != ""
replace LL = "1999.01" if LL != ""
replace LM = "1998.12" if LM != ""
replace LN = "1998.11" if LN != ""
replace LO = "1998.10" if LO != ""
replace LP = "1998.09" if LP != ""
replace LQ = "1998.08" if LQ != ""
replace LR = "1998.07" if LR != ""
replace LS = "1998.06" if LS != ""
replace LT = "1998.05" if LT != ""
replace LU = "1998.04" if LU != ""
replace LV = "1998.03" if LV != ""
replace LW = "1998.02" if LW != ""
replace LX = "1998.01" if LX != ""
replace LY = "1997.12" if LY != ""
replace LZ = "1997.11" if LZ != ""
replace MA = "1997.10" if MA != ""
replace MB = "1997.09" if MB != ""
replace MC = "1997.08" if MC != ""
replace MD = "1997.07" if MD != ""
replace ME = "1997.06" if ME != ""
replace MF = "1997.05" if MF != ""
replace MG = "1997.04" if MG != ""
replace MH = "1997.03" if MH != ""
replace MI = "1997.02" if MI != ""
replace MJ = "1997.01" if MJ != ""
replace MK = "1996.12" if MK != ""
replace ML = "1996.11" if ML != ""
replace MM = "1996.10" if MM != ""
replace MN = "1996.09" if MN != ""
replace MO = "1996.08" if MO != ""
replace MP = "1996.07" if MP != ""
replace MQ = "1996.06" if MQ != ""
replace MR = "1996.05" if MR != ""
replace MS = "1996.04" if MS != ""
replace MT = "1996.03" if MT != ""
replace MU = "1996.02" if MU != ""
replace MV = "1996.01" if MV != ""
replace MW = "1995.12" if MW != ""
replace MX = "1995.11" if MX != ""
replace MY = "1995.10" if MY != ""
replace MZ = "1995.09" if MZ != ""
replace NA = "1995.08" if NA != ""
replace NB = "1995.07" if NB != ""
replace NC = "1995.06" if NC != ""
replace ND = "1995.05" if ND != ""
replace NE = "1995.04" if NE != ""
replace NF = "1995.03" if NF != ""
replace NG = "1995.02" if NG != ""
replace NH = "1995.01" if NH != ""
replace NI = "1994.12" if NI != ""
replace NJ = "1994.11" if NJ != ""
replace NK = "1994.10" if NK != ""
replace NL = "1994.09" if NL != ""
replace NM = "1994.08" if NM != ""
replace NN = "1994.07" if NN != ""
replace NO = "1994.06" if NO != ""
replace NP = "1994.05" if NP != ""
replace NQ = "1994.04" if NQ != ""
replace NR = "1994.03" if NR != ""
replace NS = "1994.02" if NS != ""
replace NT = "1994.01" if NT != ""
replace NU = "1993.12" if NU != ""
replace NV = "1993.11" if NV != ""
replace NW = "1993.10" if NW != ""
replace NX = "1993.09" if NX != ""
replace NY = "1993.08" if NY != ""
replace NZ = "1993.07" if NZ != ""
replace OA = "1993.06" if OA != ""
replace OB = "1993.05" if OB != ""
replace OC = "1993.04" if OC != ""
replace OD = "1993.03" if OD != ""
replace OE = "1993.02" if OE != ""
replace OF = "1993.01" if OF != ""
replace OG = "1992.12" if OG != ""
replace OH = "1992.11" if OH != ""
replace OI = "1992.10" if OI != ""
replace OJ = "1992.09" if OJ != ""
replace OK = "1992.08" if OK != ""
replace OL = "1992.07" if OL != ""
replace OM = "1992.06" if OM != ""
replace ON = "1992.05" if ON != ""
replace OO = "1992.04" if OO != ""
replace OP = "1992.03" if OP != ""
replace OQ = "1992.02" if OQ != ""
replace OR = "1992.01" if OR != ""
replace OS = "1991.12" if OS != ""
replace OT = "1991.11" if OT != ""
replace OU = "1991.10" if OU != ""
replace OV = "1991.09" if OV != ""
replace OW = "1991.08" if OW != ""
replace OX = "1991.07" if OX != ""
replace OY = "1991.06" if OY != ""
replace OZ = "1991.05" if OZ != ""
replace PA = "1991.04" if PA != ""
replace PB = "1991.03" if PB != ""
replace PC = "1991.02" if PC != ""
replace PD = "1991.01" if PD != ""
replace PE = "1990.12" if PE != ""
replace PF = "1990.11" if PF != ""
replace PG = "1990.10" if PG != ""
replace PH = "1990.09" if PH != ""
replace PI = "1990.08" if PI != ""
replace PJ = "1990.07" if PJ != ""
replace PK = "1990.06" if PK != ""
replace PL = "1990.05" if PL != ""
replace PM = "1990.04" if PM != ""
replace PN = "1990.03" if PN != ""
replace PO = "1990.02" if PO != ""
replace PP = "1990.01" if PP != ""
replace PQ = "1989.12" if PQ != ""
replace PR = "1989.11" if PR != ""
replace PS = "1989.10" if PS != ""
replace PT = "1989.09" if PT != ""
replace PU = "1989.08" if PU != ""
replace PV = "1989.07" if PV != ""
replace PW = "1989.06" if PW != ""
replace PX = "1989.05" if PX != ""
replace PY = "1989.04" if PY != ""
replace PZ = "1989.03" if PZ != ""
replace QA = "1989.02" if QA != ""
replace QB = "1989.01" if QB != ""
replace QC = "1988.12" if QC != ""
replace QD = "1988.11" if QD != ""
replace QE = "1988.10" if QE != ""
replace QF = "1988.09" if QF != ""
replace QG = "1988.08" if QG != ""
replace QH = "1988.07" if QH != ""
replace QI = "1988.06" if QI != ""
replace QJ = "1988.05" if QJ != ""
replace QK = "1988.04" if QK != ""
replace QL = "1988.03" if QL != ""
replace QM = "1988.02" if QM != ""
replace QN = "1988.01" if QN != ""
replace QO = "1987.12" if QO != ""
replace QP = "1987.11" if QP != ""
replace QQ = "1987.10" if QQ != ""
replace QR = "1987.09" if QR != ""
replace QS = "1987.08" if QS != ""
replace QT = "1987.07" if QT != ""
replace QU = "1987.06" if QU != ""
replace QV = "1987.05" if QV != ""
replace QW = "1987.04" if QW != ""
replace QX = "1987.03" if QX != ""
replace QY = "1987.02" if QY != ""
replace QZ = "1987.01" if QZ != ""

foreach var of varlist P-QZ {
    gen `var'_y = real(substr(`var', 1, 4))  			// Extract year
	gen `var'_m = real(substr(`var', 6, 2)) 			// Extract month
	gen `var'_md = ym(`var'_y, `var'_m) 	// Create Stata monthly date
	format `var'_md %tm            				// Format as monthly date
	drop `var'_y `var'_m `var'
}
egen date_min = rowmin(*_md)
format date_min %tm  
egen date_max = rowmax(*_md)
format date_max %tm 
drop *_md

keep F J L M date_min date_max
rename J firmname
rename M domicile
rename L namechange
rename F GK
destring GK, replace
sort GK

duplicates report GK
duplicates tag GK, gen(dupl1)
duplicates tag GK firmname, gen(dupl2)
*br if dupl1 != dupl2
duplicates drop GK, force

merge 1:m GK using "$path\02_Processed_data\17_SIX\SIX_names_changes.dta"
drop _merge
sort GK

preserve 
keep GK domicile firmname date_*
save "$path\02_Processed_data\17_SIX\SIX_names_1.dta", replace
restore

preserve 
keep GK domicile newname date_*
rename newname firmname
save "$path\02_Processed_data\17_SIX\SIX_names_2.dta", replace
restore

preserve 
keep GK domicile oldname date_*
rename oldname firmname
save "$path\02_Processed_data\17_SIX\SIX_names_3.dta", replace
restore

* append names
use "$path\02_Processed_data\17_SIX\SIX_names_1.dta", clear
forv v = 2(1)3 {
	append using "$path\02_Processed_data\17_SIX\SIX_names_`v'.dta"
	erase "$path\02_Processed_data\17_SIX\SIX_names_`v'.dta"
}
erase "$path\02_Processed_data\17_SIX\SIX_names_1.dta"
erase "$path\02_Processed_data\17_SIX\SIX_names_changes.dta"

duplicates drop domicile firmname, force
drop if firmname == ""

sort firmname
save  "$path\02_Processed_data\17_SIX\SIX_names.dta", replace

keep if domicile == "CH"
drop domicile
sort GK firmname
save  "$path\02_Processed_data\17_SIX\SIX_names_CH.dta", replace

* Coding Bisnode:
gen nr = _n
set obs `=_N + 1'
replace nr = 0 if nr == .
sort nr
gen duns1 = .
gen duns2 = .
gen duns3 = .
gen vormduns1 = .
gen vormduns2 = .
gen vormduns3 = .
gen comment = "''duns1'' is the primary ID field, ''vormduns'' is meant as ID field for firms named as 'vorm.' in brackets" if nr == 0

order nr GK date_* firmname duns* vorm* comment
export excel "$path\02_Processed_data\17_SIX\SIX_names_CH_coding1994-2018.xlsx", firstrow(var) replace

* Coding Sugarcube:
drop duns* vorm* comment
gen ID1 = .
gen ID2 = .
gen ID3 = .
gen vormID1 = .
gen vormID2 = .
gen vormID3 = .
gen comment = "''ID1'' is the primary ID field, ''vormID'' is meant as ID field for firm IDs named as 'vorm.' in brackets" if nr == 0

order nr GK date_* firmname ID* vorm* comment
export excel "$path\02_Processed_data\17_SIX\SIX_names_CH_coding1934-2003.xlsx", firstrow(var) replace

*/

**** Prepare Datasets with SIX IDs for merge with Bisnode and Sugarcube
*** Bisnode

/*
Steps:
1) import both manual codings, merge and compare
2) Create dataset with DUNS from Bisnode and year-range of when they were quoted at SIX
*/


** 1) Read-in coded SIX data with Bisnode DUNS (1994-2018) and compare
* merge
foreach coder in ED JU { 

	import excel "$path\02_Processed_data\17_SIX\SIX_names_CH_coding1994-2018_in_`coder'.xlsx", firstrow clear  
	
// Stata loses precision in the duns variables as they are too long for the automatic storage type. It does not display values correctly.
// This is a storage issue. If I transform the imported data into "double", I can recover the correct DUNS.
	foreach col in duns1 duns2 duns3 vormduns1 vormduns2 vormduns3 {
		gen double `col'_fixed =  `col'
		drop `col'
		rename `col'_fixed `col'
		format `col' %15.0f
	}
	
	drop if nr == 0 // first raw contains information for coders
	
	foreach col in duns1 duns2 duns3 vormduns1 vormduns2 vormduns3 {		
		preserve
		keep if `col' != .
		gen double duns = `col'
		format duns %15.0f
		gen origin = "`col'"
		gen coder = "`coder'"
		drop duns1 duns2 duns3 vormduns1 vormduns2 vormduns3
		save "$path\02_Processed_data\17_SIX\\`coder'_`col'.dta", replace
		restore
	}	

	use "$path\02_Processed_data\17_SIX\\`coder'_duns1.dta", clear
	foreach col in duns2 duns3 vormduns1 vormduns2 vormduns3 {
		append using "$path\02_Processed_data\17_SIX\\`coder'_`col'.dta"
		erase "$path\02_Processed_data\17_SIX\\`coder'_`col'.dta"
	}
	rename comment comment_`coder'
	duplicates tag GK duns firmname, gen(dupl)
	duplicates drop GK duns firmname, force
	drop dupl
	sort nr
	save "$path\02_Processed_data\17_SIX\\`coder'_codings.dta", replace
	erase "$path\02_Processed_data\17_SIX\\`coder'_duns1.dta"
}
drop date_min date_max

* compare and correct
merge 1:1 nr GK firmname duns using "$path\02_Processed_data\17_SIX\ED_codings.dta"
gen check = 1 if _merge != 3
drop _merge
sort nr GK firmname duns
order check nr GK firmname duns comment_*

/* 

    Result                      Number of obs
    -----------------------------------------
    Not matched                           116
        from master                        59  (_merge==1)
        from using                         57  (_merge==2)

    Matched                               629  (_merge==3)
    -----------------------------------------

Coding looks very reasonable. Some differences occure because there are various firmnames in SIX and coders have attributed 
the same DUNS but to different SIX name variants (of the same company according to SIX ID - called GK)	

There is the problem that some holdings have various AGs related to them, which are found in Bisnode/Sugarcube. It is not always clear 
what belongs to the holding or which entry is itself the holding (e.g., ABB Group AG and ABB Holding). 

Coders have made Zefix inquiries to find "best" match. It happens that some AGs belonging to a holding are included and sometimes not.
E.g., UBS has many sub-entities that are the holding themselfs. Here, coders have only included UBS AG and UBS Group.

This coding is very difficult. Decision: stick to the coding of the coders.

*/

erase "$path\02_Processed_data\17_SIX\ED_codings.dta"
erase "$path\02_Processed_data\17_SIX\JU_codings.dta"

** 2)  Create dataset with DUNS (from Bisnode) connected to SIX firm with the years of quotation
duplicates report GK duns
duplicates drop GK duns, force
duplicates tag duns, gen(dupl)  // some duplicate duns accross SIX-ID (GK) as mergers changed SIX-ID
duplicates drop duns, force
gen SIX_min = year(date_min)
gen SIX_max = year(date_max)
gen SIX = 1
keep duns SIX*
sort duns
label var duns "ID bisnode"
label var SIX_min "start of quotation @SIX"
label var SIX_max "end of quotation @SIX"
order duns SIX SIX_*
save "$path\02_Processed_data\17_SIX\DUNS_related_SIX_1994-2018.dta", replace


** Sugarcube

/*
Steps:
1) import both manual codings, merge and compare
2) Create dataset with DUNS from Bisnode and year-range of when they were quoted at SIX
*/


** 1) Read-in coded SIX data with Sugarcube Pabel ID (1934-2003) and compare
* merge
foreach coder in ED JU { 

	import excel "$path\02_Processed_data\17_SIX\SIX_names_CH_coding1934-2003_in_`coder'.xlsx", firstrow clear  
	destring ID1, replace
	
	foreach col in ID1 ID2 ID3 vormID1 vormID2 vormID3 {
		
		gen double `col'_fixed =  `col'							// there are strings...!
		drop `col'
		rename `col'_fixed `col'
		format `col' %15.0f
	}
	drop if nr == 0 // first raw contains information for coders
	
	foreach col in ID1 ID2 ID3 vormID1 vormID2 vormID3 {		
		preserve
		keep if `col' != .
		gen double ID = `col'
		format ID %15.0f
		gen origin = "`col'"
		gen coder = "`coder'"
		drop ID1 ID2 ID3 vormID1 vormID2 vormID3
		save "$path\02_Processed_data\17_SIX\\`coder'_`col'.dta", replace
		restore
	}	

	use "$path\02_Processed_data\17_SIX\\`coder'_ID1.dta", clear
	foreach col in ID2 ID3 vormID1 vormID2 vormID3 {
		append using "$path\02_Processed_data\17_SIX\\`coder'_`col'.dta"
		erase "$path\02_Processed_data\17_SIX\\`coder'_`col'.dta"
	}
	rename comment comment_`coder'
	duplicates tag GK ID firmname, gen(dupl)
	duplicates drop GK ID firmname, force
	drop dupl
	sort nr
	save "$path\02_Processed_data\17_SIX\\`coder'_codings.dta", replace
	erase "$path\02_Processed_data\17_SIX\\`coder'_ID1.dta"
}
drop date_min date_max

* compare and correct
merge 1:1 nr GK firmname ID using "$path\02_Processed_data\17_SIX\ED_codings.dta"
gen check = 1 if _merge != 3
replace coder = "ED&JU" if _merge == 3
drop _merge M
sort nr GK firmname ID
order check nr GK firmname ID comment_*

/* 

    Result                      Number of obs
    -----------------------------------------
    Not matched                           370
        from master                       227  (_merge==1)
        from using                        143  (_merge==2)

    Matched                               706  (_merge==3)
    -----------------------------------------


Coding rather difficult as company names have changed over time. What is belonging to the SIX-quoted company? 

There is the problem that some holdings have various AGs related to them, which are found in Bisnode/Sugarcube. It is not always clear 
what belongs to the holding or which entry is itself the holding (e.g., ABB Group AG and ABB Holding). 

Coders have made Zefix inquiries to find "best" match. It happens that some AGs belonging to a holding are included and sometimes not.
E.g., UBS has many sub-entities that are the holding themselfs. Here, coders have only included UBS AG and UBS Group.

This coding is very difficult. Decision: stick to the coding of the coders and include all ID's that coders potentially link to company.

*/

erase "$path\02_Processed_data\17_SIX\ED_codings.dta"
erase "$path\02_Processed_data\17_SIX\JU_codings.dta"

** 2)  Create dataset with ID (from Sugarcube) connected to SIX firm with the years of quotation
duplicates report GK ID
duplicates drop GK ID, force
duplicates tag ID, gen(dupl)  // some duplicate IDs accross SIX-ID (GK) as mergers changed SIX-ID
duplicates drop ID, force
gen SIX_min = year(date_min)
gen SIX_max = year(date_max)
gen SIX = 1
keep ID origin SIX*
sort ID
label var ID "ID sugarcube"
label var SIX_min "start of quotation @SIX"
label var SIX_max "end of quotation @SIX"
order ID SIX SIX_*
rename ID UCID_1934_2003 // this is the name of the firm panel ID
sort UCID_1934_2003
merge 1:m UCID_1934_2003 using "$path\02_Processed_data\17_SIX\Sugarcube_Firms_Key_ID_old.dta"  // coders worked with first panel-ID (before last corrections)

/*
    Result                      Number of obs
    -----------------------------------------
    Not matched                     4,206,201
        from master                         0  (_merge==1)
        from using                  4,206,201  (_merge==2)

    Matched                            10,514  (_merge==3)
    -----------------------------------------
*/

keep if _merge == 3
drop _merge period UCID_1934_2003
order CID year
sort CID year
save "$path\02_Processed_data\17_SIX\CID_year_related_SIX_1934-2003.dta", replace














