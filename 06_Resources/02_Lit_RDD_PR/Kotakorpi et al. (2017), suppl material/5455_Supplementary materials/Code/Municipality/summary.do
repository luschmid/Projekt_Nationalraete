* SUMMARY STATISTICS - Municipal elections

local sumvars female age /*yob*/ incumbent voteshare KOK KESK SDP nseats nvotes ///
    found avgtulo avgpotulo next1_3tulo next5_7tulo next9_11tulo last1_3tulo

order `sumvars'
sort elected
quietly bys elected: outreg2 `sumvars' using $OUTPATH/summary, ///
  tex sum(log) eqkeep (N mean sd) label addn("Incomes in 2011 euros per annum.") ///
  title("Table 1. Elected and defeated candidates in the 1996-2008 municipal elections") replace
