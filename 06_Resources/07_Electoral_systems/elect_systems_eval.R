

setwd("C:/Users/julia/OneDrive/Desktop/Research Assistant/Democratic Electoral Systems around the world")

library(readr)
library(car)

# Daten einlesen

elect_systems_data <- read_csv("es_data-v3.csv")

str(elect_systems_data)

elect_systems_data <- elect_systems_data[, -c(28, 29, 30, 31, 32, 33)]

# Numerische Variablen für Lesbarkeit recodieren

elect_systems_data$legislative_elect_system <- 
  recode(elect_systems_data$legislative_type, "1='Majoritarian'; 2='Proportional'; 3='Mixed'; 
         -99='Not available'")


elect_systems_data$electoral_rule <- 
  recode(elect_systems_data$elecrule, "1='Single-Member-District-Plurality';
         2='Two-Round-System'; 3='Alternative Vote'; 4='Borda Count'; 5='Block Vote'; 6='Party Block Vote';
         7='Limited Vote'; 8='Single-Nontransferable Vote'; 9='List Proportional Representation';
         10='Single-Transferable Vote'; 11='Mixed-Member-Proportional'; 12='Mixed Parallel'; 
         -99='Not available'")

# Anmerkung: Tier1 ist immer auf nationaler Ebene und 2-4 sind immer kleinere Wahleinheiten
# (Aber es ist ohnehin nur die nationale Ebene wirklih gut abgedeckt - trotzdem habe ich tier 2 mal auch recodiert)

elect_systems_data$PR_Seat_calculation_national <- 
  recode(elect_systems_data$tier1_formula, "c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)='Not applicable';
         12='Hare Quota'; 13='Hare Quota (with largest remainder)'; 
         14='Hare Quota (highest average remainder)'; 15='Hagenbach_Bischoff Quota'; 
         16='Hagenbach-Bischoff Quota (with largest remainder)';
         17='Hagenbach-Bischoff Quota (with highest average remainder)'; 18='Droop Quota'; 
         19='Droop Quota (with largest remainder)';
         20='Droop Quota (with highest average remainder)'; 21='Imperiali Quota'; 
         22='Imperiali Quota (with largest remainder)'; 23='imperiali Quota (with highest average remainder)';
         24='Reinforced Imperiali Quota'; 25='d`Hondt'; 26='Sainte-Laguë'; 27='Modified Sainte-Laguë';
         28='Not applicable'; -99='Not available'; -88='Not available'")


elect_systems_data$PR_Seat_calculation_regional <- 
  recode(elect_systems_data$tier2_formula, "c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)='Not applicable';
         12='Hare Quota'; 13='Hare Quota (with largest remainder)'; 
         14='Hare Quota (highest average remainder)'; 15='Hagenbach_Bischoff Quota'; 
         16='Hagenbach-Bischoff Quota (with largest remainder)';
         17='Hagenbach-Bischoff Quota (with highest average remainder)'; 18='Droop Quota'; 
         19='Droop Quota (with largest remainder)';
         20='Droop Quota (with highest average remainder)'; 21='Imperiali Quota'; 
         22='Imperiali Quota (with largest remainder)'; 23='imperiali Quota (with highest average remainder)';
         24='Reinforced Imperiali Quota'; 25='d`Hondt'; 26='Sainte-Laguë'; 27='Modified Sainte-Laguë';
         28='Not applicable'; -99='Not available'; -88='Not available'")


Nice_data <- elect_systems_data[, c(1, 2, 3, 4, 36, 37, 38, 39, 30, 31, 32)]









