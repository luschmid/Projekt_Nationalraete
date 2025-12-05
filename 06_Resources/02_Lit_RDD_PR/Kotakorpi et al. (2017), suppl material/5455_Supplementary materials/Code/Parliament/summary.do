/*SUMMARY STATISTICS*/
local sumvars female age /*yob*/ incumbent prevelected voteshare  nseats nvotes south KOK KESK SDP OTHER ///
    found avgtulo avgpotulo next1_3tulo next5_7tulo next9_11tulo last1_3tulo

order `sumvars'
preserve
keep elected `sumvars'

sort elected
quietly bys elected: outreg2 `sumvars' using $OUTPATH/summary, ///
  tex sum(log) eqkeep (N mean sd) label addn("Incomes in 2011 euros per annum.") ///
  title("Table 1. Elected and defeated candidates in the 1970-2007 parliamentary elections") replace

restore


