global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte\01_Raw_data\02_Elections_1971_2015\01_Candidates_2019_Election"

* AG
import exc using "$path\AG\Kanton Aargau NR Wahlen 2019.xlsx", ///
	sh(KKanton) cellrange(A1:M497) first clear
desc, s
duplicates report AmtlName AmtlVorname
tab bisher

* AI
import exc using "$path\AI\Daten aufbereitet 20190914.xlsx", ///
	sh(Tabelle1) cellrange(A1:E5) first clear
desc, s
tab Bisher

* AR
import exc using "$path\AR\Daten aufbereitet 20190914.xlsx", ///
	sh(Tabelle1) cellrange(A1:F3) first clear
desc, s
tab Bisher

* BL
import exc using "$path\BL\Daten aufbereitet 20190914.xlsx", ///
	sh(Tabelle1) cellrange(A1:G136) first clear
duplicates drop KandNr Name, force
desc, s
tab Bisher

* BS
import delim using "$path\BS\BS_NR-Wahl.csv", delim(";") clear
duplicates drop kandidatennr vorname name, force
desc, s
tab bisher

* BE
import exc using "$path\BE\Daten aufbereitet 20190912.xlsx", ///
	sh(reportKandidatenListe-A) cellrange(A1:K755) first clear
duplicates drop Name Vorname GebJahr, force
desc, s
tab Bisher

* FR
import exc using "$path\FR\Daten aufbereitet 20190908.xlsx", ///
	sh(de_Prop.BereinigteWahlvorschläg) cellrange(A1:M161) first clear
desc, s
duplicates drop Name Vorname, force
desc, s
tab Bisher

* GE
import exc using "$path\GE\Daten aufbereitet 20190911.xlsx", ///
	sh(Tabelle1) cellrange(A1:G177) first clear
duplicates tag NomPrénom, g(dups)
list NoListe Nom NomPrénom Sexe if dups>0
desc, s
tab Bisher 
* 176 candidates according to official website

* GL
import exc using "$path\GL\Daten aufbereitet 20190913.xlsx", ///
	sh(Tabelle1) cellrange(A1:D3) first clear
desc, s
tab Bisher

* GR
import exc using "$path\GR\Daten aufbereitet 20190908.xlsx", ///
	sh(Tabelle1) cellrange(A1:I101) first clear
duplicates report Kandidat
desc, s
tab Bisher

* JU
import exc using "$path\JU\Daten aufbereitet 20190913.xlsx", ///
	sh(Tabelle1) cellrange(A1:G35) first clear
desc, s
tab Bisher

* LU
import exc using "$path\LU\Kanton Luzern NR Wahlen 2019.xlsx", ///
	sh(Tabelle1) cellrange(A1:L253) first clear
desc, s
tab Bisher

* NE
import exc using "$path\NE\Daten aufbereitet 20190911.xlsx", ///
	sh(Tabelle1) cellrange(A1:G47) first clear
duplicates report Nom Prénom
desc, s
tab Bisher

* NW
import exc using "$path\NW\Daten aufbereitet 20190913.xlsx", ///
	sh(Tabelle1) cellrange(A1:C3) first clear
desc, s
tab Bisher

* OW
import exc using "$path\OW\Kanton Obwalden NR Wahlen 2019.xlsx", ///
	sh(Tabelle1) cellrange(A3:E11) clear
rename A Anrede
rename B Kandidat
rename C Adresse
rename D PLZ
rename E Wohnort
drop if Anrede == ""
desc, s

* SH
	
* SZ
import exc using "$path\SZ\Kandidaten NR-Wahlen 2019.xlsx", ///
	sh(WABSTI Datenexport) cellrange(A1:J85) first clear
desc, s
duplicates report Nachname Vorname
tab LiNr

* SO
import exc using "$path\SO\Daten aufbereitet 20190913.xlsx", ///
	sh(Tabelle1) cellrange(A1:J167) first clear
desc, s
tab bisher

* SG
import exc using "$path\SG\Daten bearbeitet 20190909.xlsx", ///
	sh(Listen und Kandidaten) cellrange(A1:I283) first clear
duplicates drop ListenundKandidierennummer, force
desc, s
tab bisher

* TI
import exc using "$path\TI\Daten aufbereitet 20190910.xlsx", ///
	sh(Tabelle1) cellrange(A1:E157) first clear
duplicates drop CognomeNome Lista, force
desc, s
tab Bisher

* TG
import exc using "$path\TG\Kanton ThurgauNR Wahlen 2019.xlsx", ///
	sh(Table 1) cellrange(A1:J139) first clear
duplicates drop KandNr Name Vorname, force
desc, s
tab Bisher

* UR
import exc using "$path\UR\Daten aufbereitet 20190913.xlsx", ///
	sh(Tabelle1) cellrange(A1:C5) first clear
desc, s
tab Bisher

* VD
import exc using "$path\VD\Daten aufbereitet 20190910.xlsx", ///
	sh(Export candidats) cellrange(A1:I375) first clear
duplicates report Nomusuel Prénomusuel Numérocandidat
tab Bisher

* VS

* ZG
import exc using "$path\ZG\Daten aufbereitet 20190914.xlsx", ///
	sh(Tabelle1) cellrange(A1:H76) first clear
duplicates report KandNr Name
tab H

* ZH
import exc using "$path\ZH\Provisorisch_AlleKand_NR2019_20190816_Online.xlsx", ///
	sh(WABSTI Datenexport) cellrange(A1:M1023) first clear
desc, s
duplicates drop LiNr KNr, force
desc, s
tab Bisher
