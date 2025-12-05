

*******************************************************************
** Paper: Curto-Grau, M., Solé-Ollé, A., and Sorribas-Navarro, P. 2018. "Does electoral competition curb party favoritism?", AEJ:AE
** Code to create the forcing variable based on Hypothesis 2 (part of the votes lost by the regional incumbent's bloc go to abstention and part to the opposition)
*******************************************************************

// Note 1 : the allocation of seats is done using program v2seats. It can be found on http://ideas.repec.org/c/boc/bocode/s456973.html
// Note 2 : in Spain the minimum threshold to obtain representation in the city council is 5% of valid votes
// Download Stata module "electool":  https://ideas.repec.org/c/boc/bocode/s456973.html


use "file with electoral results", clear
** with the initial distribution of votes amongst parties, generate the allocation of seats using the Stata program v2seats
qui v2seats vloc if sval>0.05, party(pol_party) formula(dhondt) district(mun_code) size(stotal) ///
	save("new file seats.dta") preserve(mun_code pol_party)

** merge initial database ("file with electoral results") with the database with the allocation of seats ("new file seats.dta"), then:
use "file with electoral results and seats_method2", clear	
** calculate the initial number of seats that belong to right and left-wing parties (named "seatsleft0" and "seatsright0")
** define threshold to have a majority of seats
gen seatmaj=floor((stotal/2)+1)
** using the previous threshold define whether the municipality has a majority of left or right-wing seats
gen majleft0=seatsleft0>=seatmaj
gen majright0=seatsright0>=seatmaj
*********************************
** computing comparison numbers
*********************************
gen cl=vloc0/sloc if left==1
gen cr=vloc0/sloc if right==1
gen cl_1=vloc0/(sloc+1) if (right==0 & sval>0.05)
gen cr_1=vloc0/(sloc+1) if (left==0  & sval>0.05)
sort mun_code
by mun_code: egen clmax=max(cl_1)
by mun_code: egen crmax=max(cr_1)
replace clmax=0 if clmax==.
replace crmax=0 if crmax==.
** in case that either clmax or crmax =0, then: 
ge id=0
replace id=1 if (clmax==0 | crmax==0)
replace cl_1=vloc0/(sloc+1) if (right==0 & id==1)
replace cr_1=vloc0/(sloc+1) if (left==0 & id==1)
sort mun_code
drop clmax crmax
by mun_code: egen clmax=max(cl_1)
by mun_code: egen crmax=max(cr_1)
replace clmax=0 if clmax==.
replace crmax=0 if crmax==.
** Note: up until this point the procedure is the same as in "method1"
*********************************
** compute the distance if there's a left-wing majority in the city council (note that deltai=votes that each left-w party loses, delta=total votes lost by the left, etai=votes that each right-w party wins)
*********************************
gen x1=(cl-crmax)*sloc0 if majleft0==1
gen sloc_rrr=sloc0 if cr_1==crmax
by codiine: egen sloc_rr=max(sloc_rrr)
gen alpha_rrr=alpha_r if cr_1==crmax & sloc0==sloc_rr
by codiine: egen alpha_rr=max(alpha_rrr)
gen delta_00=((x1/alpha_l)/(1+((sloc0/(sloc_rr+1))*(alpha_rr/alpha_l)*fi_r)))+1   if majleft0==1
by codiine: egen delta0=min(delta_00) if majleft0==1
gen delta1=-delta0 if majleft0==1
gen delta1i=delta1*alpha_l if majleft0==1
gen eta1i=delta0*alpha_r*fi_r
*********************************
** compute the distance if there's a right-wing majority in the city council(note that betai=votes that each right-w party loses, beta=total votes lost by the right, mui=votes that each left-wing party wins)
*********************************
gen y1=(cr-clmax)*sloc0 if majright0==1
gen sloc_lll=sloc0 if cl_1==clmax
by codiine: egen sloc_ll=max(sloc_lll)
gen alpha_lll=alpha_l if cl_1==clmax & sloc0==sloc_ll
by codiine: egen alpha_ll=max(alpha_lll)
gen beta_00=((y1/alpha_r)/(1+((sloc0/(sloc_ll+1))*(alpha_ll/alpha_r)*fi_l)))+1   if majright0==1
by codiine: egen beta0=min(beta_00)  if majright0==1
gen beta1=-beta0 if majright0==1
gen beta1i=beta1*alpha_r if majright0==1
gen mu1i=beta0*alpha_l*fi_l
** round betai deltai, etai and mui 
replace delta1i=floor(delta1i)
replace beta1i=floor(beta1i)
replace eta1i=ceil(eta1i)
replace mu1i=ceil(mu1i)
** check that x,y,etai and mui are positive and delta1 and beta1 negative
sum x1 y1 delta1 beta1 eta1i mu1i
*********************************
** new distribution of votes
*********************************
gen vloc1=vloc
replace vloc1=vloc+delta1i if delta1i!=. & majleft0==1 & left==1
replace vloc1=vloc+beta1i if  beta1i!=. & majright0==1 & right==1
replace vloc1=vloc+eta1i if eta1i!=. & majleft0==1 & right==1
replace vloc1=vloc+mu1i if  mu1i!=. & majright0==1 & left==1
** if, due to rounding delta and beta, we end up with a negative amount of votes:
replace vloc1=0 if vloc1<0
** if two parties have the same amount of votes (different than 0), v2seats cannot compute correctly (thus, drop this observations if they exist). Then:
sort mun_code
by codiine: egen vtot1=sum(vloc1)
gen sval1=vloc1/(vtot1+vblanc)
save "file with electoral results and seats1_method2", replace
*********************************
** new distribution of seats
*********************************
** once delta, beta, mu and eta are added, with the new allocation of votes generate the new allocation of seats
qui v2seats vloc if sval>0.05, party(pol_party) formula(dhondt) district(mun_code) size(stotal)  ///
	save("new file seats.dta") preserve(mun_code pol_party)
** merge previous database ("file with electoral results and seats1") with the database with the new allocation of seats ("new file seats.dta"), then:
use "file with electoral results and seats1_method2", clear	
*********************************
** Checks
*********************************
** calculate the new number of total seats that belong to right and left-wing parties (named "seatsleft1" and "seatsright1")
** perform the following final checks (note that seats1 is the new allocation of seats):
// 1. The initial amount of total seats has to remain unchanged
sort mun_code
by mun_code: egen seats0=sum(sloc0)
by mun_code: egen seats1=sum(sloc1)
ge check_seats0=stotal-seats0
ge check_seats1=stotal-seats1
tab check_seats0	
tab check_seats1
// 2. How many seats have switched? (it should only be one to ensure the minimum distance is computed)
ge s_lost1=sloc1-sloc0
ge s_lost_left1=seatleft1-seatleft0
ge s_lost_right1=seatright1-seatright0
tab s_lost1
tab s_lost_left1
tab s_lost_right1
// 3. Does the party that initially had a majority keep it?
gen majleft1=seatleft1>=seatmaj
gen majright1=seatright1>=seatmaj

** FINAL note: If there has been a majority switch (compared to the initial situation), 
// then the process stops here. If not, the procedure is iterated until all municipalities that had initially a left-wing majority (majleft0==1) have finally 
// a right-wing majority, and the municipalities that initially had a right-wing majority (majright0==1) have finally a left-wing majority
tab majleft1 if majleft0==1
tab majright1 if majright0==1

*********************************
** Computing the FINAL DISTANCE (forcing variable, labelled "distance2")
*********************************
replace delta1i=0 if delta1i==.
replace beta1i=0 if beta1i==.
by mun_code: egen tdelta1=sum(delta1i)
by mun_code: egen tbeta1=sum(beta1i)
gen delta=abs(tdelta1)/vtotal if majleft0==1
gen beta=abs(tbeta1)/vtotal if majright0==1
** delta and beta cannot be >1 as the votes lost cannot be greater than the initial amount of votes (if this happens it's due to rounding. Such observations have to be dropped)
sum delta beta

** If the upper-level government is right-wing and the city council has a left-wing majority (i.e. they are likely to be unaligned)
gen distance2=-delta if majleft0==1 & ul_right==1
** If the upper-level government is left-wing and the city council has a left-wing majority (i.e. they are likely to be aligned)
replace distance2=delta if majleft0==1 &  ul_right==0
** If the upper-level government is left-wing and the city council has a right-wing majority (i.e. they are likely to be unaligned)
replace distance2=-beta if majleft0==0 &  ul_right==0
** If the upper-level government is right-wing and the city council has a right-wing majority (i.e. they are likely to be aligned)
replace distance2=beta if majleft0==0 & ul_right==1
