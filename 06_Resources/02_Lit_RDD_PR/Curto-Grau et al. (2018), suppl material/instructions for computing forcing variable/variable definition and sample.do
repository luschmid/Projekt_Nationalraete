************************
** Variable definition
************************
vloc0 = initial number of votes for party x in local elections
sloc0 = initial number of seats of party x in the city council
sval = share of valid votes in local elections
pol_party = political party running 
left = 1 if the political party is left-wing
right = 1 if the political party is right-wing
mun_code = municipality code
stotal = total amount of seats distributed in the city council
vtotal = total amount of votes cast in the local election
seatsleft0 = initial number of seats that belong to left-wing parties in the city council
seatsright0 = initial number of seats that belong to right-wing parties in the city council
alpha_l = votes that party x (which is left-wing) has as a share of the total amount of votes of the left-wing bloc in the municipality
alpha_r = votes that party y (which is right-wing) has as a share of the total amount of votes of the right-wing bloc in the municipality
vblank = number of blank ballots in local elections
ul_right = 1 if the upper-level government is right-wing

** To compute the distance using method2, the following variables are needed:
//1. Generate a share of right-wing votes, left-wing votes and abstention
*note voteleft0 = initial amount of left-wing votes in the city council, voteright0 = initial amount of right-wing votes in the city council
gen ro_l=voteleft0/census
gen ro_r=voteright0/census
gen ro_a=(census-voters)/census
* note: voters = total amount of votes cast in local elections in municipality X , and census = population census
gen fi_l=ro_l/(ro_a+ro_l)
gen fi_r=ro_r/(ro_a+ro_r)




************************
** Sample restrictions
************************

1. We dropped municipalities where 2 parties have exactly the same amount of votes because the Stata command v2seats fails)
2. We dropped municipalities where all parties are right-wing or all parties are left-wing

These 2 situations happen rarely in our sample.
