	* (iv) Explore missing values
	
		foreach var of varlist nationalratswahl-merge_cantonaltranslator{
		display "`var'" 
		capture  noisily count if `var'!=""
		capture  noisily count if `var'!=.	
		}
		