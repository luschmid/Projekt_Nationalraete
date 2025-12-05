version 18
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"

use "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_corr.dta", clear

*** Prepare dataset for coding of SIX quotation

keep duns firma rechtsform filialeindikator gruendungsjahr
tab rechtsform
/*
. tab rechtsform

                             Rechtsform |      Freq.     Percent        Cum.
----------------------------------------+-----------------------------------
                     Aktiengesellschaft |    362,600       28.71       28.71
                                Anstalt |          3        0.00       28.71
                    Duns Support Record |          1        0.00       28.71
                  Einfache Gesellschaft |      2,489        0.20       28.90
                            Einzelfirma |    492,828       39.01       67.92
                FL - Aktiengesellschaft |          3        0.00       67.92
                         Genossenschaft |     15,668        1.24       69.16
  Gesellschaft mit beschränkter Haftung |    270,555       21.42       90.58
            Haupt von Gemeinderschaften |         55        0.00       90.58
Investmentgesellschaft mit variablem .. |         14        0.00       90.58
                  Kollektivgesellschaft |     41,076        3.25       93.83
                  Kommanditgesellschaft |      5,076        0.40       94.24
Kommanditgesellschaft für kollektive .. |         16        0.00       94.24
                   Rechtsform unbekannt |         53        0.00       94.24
                               Stiftung |     29,050        2.30       96.54
                                 Verein |     14,312        1.13       97.67
                     Zweigniederlassung |     16,705        1.32       99.00
Zweigniederlassung (Hauptsitz im Ausl.. |     11,420        0.90       99.90
      Öffentlich-rechtliche Institution |      1,262        0.10      100.00
----------------------------------------+-----------------------------------
                                  Total |  1,263,186      100.00
*/

keep if rechtsform == "Aktiengesellschaft" | rechtsform == "Öffentlich-rechtliche Institution"

tab filialeindikator
/*
. tab filialeindikator

    Filiale |
  Indikator |      Freq.     Percent        Cum.
------------+-----------------------------------
          N |    363,862      100.00      100.00
------------+-----------------------------------
      Total |    363,862      100.00
*/
drop filialeindikator
sort firma
order duns rechtsform gruendungsjahr firma
save "$path\02_Processed_data\17_SIX\Bisnode_firms_AGs.dta", replace
export excel "$path\02_Processed_data\17_SIX\Bisnode_firms_AGs.xlsx", firstrow(var) replace