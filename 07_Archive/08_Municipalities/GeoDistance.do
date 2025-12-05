
******************************************************
* (0) Set directories
******************************************************

	clear
	set more off

	capture cd "D:\SchmidLu\Dropbox\Projekt Nationalräte\Verschiedenes\GeoDistance"
	if _rc==0{
		global hauptpfad "D:\SchmidLu\Dropbox\Projekt Nationalräte\Verschiedenes\GeoDistance"

	}
	


	******************************************************
	* (A) READ-IN DATA
	******************************************************

		use "$hauptpfad\municipality_coordinates.dta", clear
		
 		levelsof gdenr_2018, local(levels) 
		foreach i of local levels {		
		use "$hauptpfad\municipality_coordinates.dta", clear
		gen gdenr_2018_destination=`i'
		if `i'==1{
		save "$hauptpfad\temp.dta", replace
		}
		else{
		append using "$hauptpfad\temp.dta"
		save "$hauptpfad\temp.dta", replace
		}
		}
		
		rename gdenr_2018 gdenr_orig
		rename X_CNTR x_orig
		rename Y_CNTR y_orig
		rename gdenr_2018_destination gdenr_dest
		rename GMDNAME name_orig
		
		preserve
		use "$hauptpfad\municipality_coordinates.dta", clear
		rename gdenr_2018 gdenr_dest
		rename X_CNTR x_dest
		rename Y_CNTR y_dest
		rename GMDNAME name_dest		
		save "$hauptpfad\temp2.dta", replace
		restore
		
		merge m:1 gdenr_dest using "$hauptpfad\temp2.dta", gen(merge_orig_dest)
		drop merge_orig_dest
		
		sort gdenr_orig gdenr_dest
		
		//gen gdenr_2018_origin= gdenr_2018 
		//gen gdenr_2018_destination= gdenr_2018 
		
		//joinby gdenr_2018_origin  gdenr_2018_destination using "$hauptpfad/municipality_coordinates.dta"
		//rename gdenr_2018 gdenr_2018_destination
		save "$hauptpfad\municipality_coordinates_origin_destination.dta", replace
		erase "$hauptpfad\temp.dta"
		erase "$hauptpfad\temp2.dta"

		
		
		******************************************************
		* (B) CALCULATE REAL AND PSEUDO-DISTANCES
		******************************************************
		
		use "$hauptpfad\municipality_coordinates_origin_destination.dta",clear
		
		
		g x2_orig = x_orig * x_orig
		g y2_orig = y_orig * y_orig
		g x2_dest = x_dest * x_dest
		g y2_dest = y_dest * y_dest
		
		g xsr_orig = sqrt(x_orig)
		g ysr_orig = sqrt(y_orig)
		g xsr_dest = sqrt(x_dest)
		g ysr_dest = sqrt(y_dest)

		gen real_dist=sqrt((x_orig-x_dest)^2+(y_orig-y_dest)^2)
		
		gen pseudo_dist1=abs(x_orig-x_dest)+abs(y_orig-y_dest)

		gen pseudo_dist2=abs(xsr_orig-xsr_dest)+abs(ysr_orig-ysr_dest)
		
		gen x_y_ratio=abs(x_orig-x_dest)/abs(y_orig-y_dest)
		sum x_y_ratio, d
		
		

		cor real_dist pseudo_* 
		
		cor real_dist pseudo_*  if x_y_ratio==. & gdenr_orig!=gdenr_dest
		
		cor real_dist pseudo_*  if x_y_ratio<1.42704 
		cor real_dist pseudo_*  if x_y_ratio>1.42704 
		
		
		cor real_dist pseudo_*  if inrange(x_y_ratio,0.5,1.5)
		cor real_dist pseudo_*  if !inrange(x_y_ratio,0.5,1.5)		
