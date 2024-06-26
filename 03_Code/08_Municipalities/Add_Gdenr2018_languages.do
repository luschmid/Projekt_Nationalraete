clear 	
cap log close
set more off	  

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"


use "$path\01_Raw_data\08_Municipalities\BfS_Municipalities_Language.dta", clear

preserve
use "$path\02_Processed_data\08_Municipalities\Gde_1960_2018.dta"
keep if year== 2017
keep gdenr gdename gdenr_2018
rename gdenr gem_cd
duplicates report gdenr
sort gem_cd
save "$path\02_Processed_data\08_Municipalities\Gde_2017_temp.dta", replace
restore

merge 1:1 gem_cd using "$path\02_Processed_data\08_Municipalities\Gde_2015_temp.dta"
replace gdenr_2018 = gem_cd if _merge==1 // manually checked
replace gem_name = gdename if _merge == 2
replace language = "Deutsches Sprachgebiet" if gem_name == "Illnau-Effretikon" & gdenr_2018 == 296
replace language = "Deutsches Sprachgebiet" if gem_name == "Kyburg" & gdenr_2018 == 296
replace language = "Deutsches Sprachgebiet" if gem_name == "Niederösch" & gdenr_2018 == 405
replace language = "Deutsches Sprachgebiet" if gem_name == "Oberösch" & gdenr_2018 == 405
replace language = "Deutsches Sprachgebiet" if gem_name == "Bangerten" & gdenr_2018 == 310
replace language = "Deutsches Sprachgebiet" if gem_name == "Tägertschi" & gdenr_2018 == 616
replace language = "Deutsches Sprachgebiet" if gem_name == "Hermiswil" & gdenr_2018 == 988
replace language = "Französisches Sprachgebiet" if gem_name == "Bussy (FR)" & gdenr_2018 == 2054
replace language = "Französisches Sprachgebiet" if gem_name == "Châbles" & gdenr_2018 == 2055
replace language = "Französisches Sprachgebiet" if gem_name == "Cheyres" & gdenr_2018 == 2055
replace language = "Französisches Sprachgebiet" if gem_name == "Domdidier" & gdenr_2018 == 2053
replace language = "Französisches Sprachgebiet" if gem_name == "Dompierre (FR)" & gdenr_2018 == 2053
replace language = "Französisches Sprachgebiet" if gem_name == "Estavayer-le-Lac" & gdenr_2018 == 2054
replace language = "Französisches Sprachgebiet" if gem_name == "Léchelles" & gdenr_2018 == 2053
replace language = "Französisches Sprachgebiet" if gem_name == "Morens (FR)" & gdenr_2018 == 2054
replace language = "Französisches Sprachgebiet" if gem_name == "Murist" & gdenr_2018 == 2054
replace language = "Französisches Sprachgebiet" if gem_name == "Rueyres-les-Prés" & gdenr_2018 == 2054
replace language = "Französisches Sprachgebiet" if gem_name == "Russy" & gdenr_2018 == 2053
replace language = "Französisches Sprachgebiet" if gem_name == "Villeneuve (FR)" & gdenr_2018 == 2044
replace language = "Französisches Sprachgebiet" if gem_name == "Vuissens" & gdenr_2018 == 2054
replace language = "Französisches Sprachgebiet" if gem_name == "Vernay" & gdenr_2018 == 2054
replace language = "Französisches Sprachgebiet" if gem_name == "Autafond" & gdenr_2018 == 2175
replace language = "Französisches Sprachgebiet" if gem_name == "Chésopelloz" & gdenr_2018 == 2183
replace language = "Französisches Sprachgebiet" if gem_name == "Corpataux-Magnedens" & gdenr_2018 == 2236
replace language = "Französisches Sprachgebiet" if gem_name == "Farvagny" & gdenr_2018 == 2236
replace language = "Französisches Sprachgebiet" if gem_name == "Rossens (FR)" & gdenr_2018 == 2236
replace language = "Französisches Sprachgebiet" if gem_name == "Le Glèbe" & gdenr_2018 == 2236
replace language = "Französisches Sprachgebiet" if gem_name == "Vuisternens-en-Ogoz" & gdenr_2018 == 2236
replace language = "Französisches Sprachgebiet" if gem_name == "Barberêche" & gdenr_2018 == 2254
replace language = "Französisches Sprachgebiet" if gem_name == "Courlevon" & gdenr_2018 == 2275
replace language = "Deutsches Sprachgebiet" if gem_name == "Jeuss" & gdenr_2018 == 2275
replace language = "Deutsches Sprachgebiet" if gem_name == "Lurtigen" & gdenr_2018 == 2275
replace language = "Deutsches Sprachgebiet" if gem_name == "Salvenach" & gdenr_2018 == 2275
replace language = "Französisches Sprachgebiet" if gem_name == "Villarepos" & gdenr_2018 == 2254
replace language = "Französisches Sprachgebiet" if gem_name == "Bas-Vully" & gdenr_2018 == 2284
replace language = "Französisches Sprachgebiet" if gem_name == "Haut-Vully" & gdenr_2018 == 2284
replace language = "Deutsches Sprachgebiet" if gem_name == "Wallenried" & gdenr_2018 == 2254
replace language = "Deutsches Sprachgebiet" if gem_name == "Oberschrot" & gdenr_2018 == 2299
replace language = "Deutsches Sprachgebiet" if gem_name == "Zumholz" & gdenr_2018 == 2299
replace language = "Rätoromanisches Sprachgebiet" if gem_name == "Bivio" & gdenr_2018 == 3543
replace language = "Rätoromanisches Sprachgebiet" if gem_name == "Cunter" & gdenr_2018 == 3543
replace language = "Rätoromanisches Sprachgebiet" if gem_name == "Marmorera" & gdenr_2018 == 3543
replace language = "Rätoromanisches Sprachgebiet" if gem_name == "Mulegns" & gdenr_2018 == 3543
replace language = "Rätoromanisches Sprachgebiet" if gem_name == "Riom-Parsonz" & gdenr_2018 == 3543
replace language = "Rätoromanisches Sprachgebiet" if gem_name == "Salouf" & gdenr_2018 == 3543
replace language = "Rätoromanisches Sprachgebiet" if gem_name == "Savognin" & gdenr_2018 == 3543
replace language = "Rätoromanisches Sprachgebiet" if gem_name == "Sur" & gdenr_2018 == 3543
replace language = "Rätoromanisches Sprachgebiet" if gem_name == "Tinizong-Rona" & gdenr_2018 == 3543
replace language = "Deutsches Sprachgebiet" if gem_name == "Obersaxen" & gdenr_2018 == 3988
replace language = "Deutsches Sprachgebiet" if gem_name == "Mundaun" & gdenr_2018 == 3988
replace language = "Italienisches Sprachgebiet" if gem_name == "Leggia" & gdenr_2018 == 3832
replace language = "Italienisches Sprachgebiet" if gem_name == "Verdabbio" & gdenr_2018 == 3832
replace language = "Deutsches Sprachgebiet" if gem_name == "Saas" & gdenr_2018 == 3871
replace language = "Deutsches Sprachgebiet" if gem_name == "St. Antönien" & gdenr_2018 == 3891
replace language = "Italienisches Sprachgebiet" if gem_name == "Sobrio" & gdenr_2018 == 5072
replace language = "Italienisches Sprachgebiet" if gem_name == "Gresso" & gdenr_2018 == 5136
replace language = "Italienisches Sprachgebiet" if gem_name == "Mosogno" & gdenr_2018 == 5136
replace language = "Italienisches Sprachgebiet" if gem_name == "Vergeletto" & gdenr_2018 == 5136
replace language = "Italienisches Sprachgebiet" if gem_name == "Isorno" & gdenr_2018 == 5136
replace language = "Französisches Sprachgebiet" if gem_name == "Brenles" & gdenr_2018 == 5675
replace language = "Französisches Sprachgebiet" if gem_name == "Chesalles-sur-Moudon" & gdenr_2018 == 5675
replace language = "Französisches Sprachgebiet" if gem_name == "Cremin" & gdenr_2018 == 5675
replace language = "Französisches Sprachgebiet" if gem_name == "Forel-sur-Lucens" & gdenr_2018 == 5675
replace language = "Französisches Sprachgebiet" if gem_name == "Sarzens" & gdenr_2018 == 5675
replace language = "Französisches Sprachgebiet" if gem_name == "Corcelles-sur-Chavornay" & gdenr_2018 == 5749
replace language = "Französisches Sprachgebiet" if gem_name == "Carrouge (VD)" & gdenr_2018 == 5806
replace language = "Französisches Sprachgebiet" if gem_name == "Ferlens (VD)" & gdenr_2018 == 5806
replace language = "Französisches Sprachgebiet" if gem_name == "Mézières (VD)" & gdenr_2018 == 5806
replace language = "Französisches Sprachgebiet" if gem_name == "Essert-Pittet" & gdenr_2018 == 5749
replace language = "Deutsches Sprachgebiet" if gem_name == "Blitzingen" & gdenr_2018 == 6077
replace language = "Deutsches Sprachgebiet" if gem_name == "Niederwald" & gdenr_2018 == 6077
replace language = "Deutsches Sprachgebiet" if gem_name == "Grafschaft" & gdenr_2018 == 6077
replace language = "Deutsches Sprachgebiet" if gem_name == "Münster-Geschinen" & gdenr_2018 == 6077
replace language = "Deutsches Sprachgebiet" if gem_name == "Reckingen-Gluringen" & gdenr_2018 == 6077
replace language = "Französisches Sprachgebiet" if gem_name == "Les Agettes" & gdenr_2018 == 6266
replace language = "Französisches Sprachgebiet" if gem_name == "Chermignon" & gdenr_2018 == 6253
replace language = "Französisches Sprachgebiet" if gem_name == "Mollens (VS)" & gdenr_2018 == 6253
replace language = "Französisches Sprachgebiet" if gem_name == "Montana" & gdenr_2018 == 6253
replace language = "Französisches Sprachgebiet" if gem_name == "Randogne" & gdenr_2018 == 6253
replace language = "Französisches Sprachgebiet" if gem_name == "Brot-Dessous" & gdenr_2018 == 6413

drop _merge gdename

sort gdenr_2018
save "$path\02_Processed_data\08_Municipalities\BfS_Municipalities_Language_gdenr2018.dta", replace
erase "$path\02_Processed_data\08_Municipalities\Gde_2017_temp.dta"
