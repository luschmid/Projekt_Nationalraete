

cd "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data"
use nationalraete_1931_2015_Geo.dta, replace


keep ID birthyear firstname name sex municipality municipalityno E_CNTR_w N_CNTR_w

drop if municipalityno==.

cd "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data\17_Municipality_GdeNo_PLZ_Geo\data"
export delimited nationalraete_1931_2015_Export_RecordLinkage.csv, replace  delimiter(";")
