
* This code produces a lookup table for the ksmirnovk test for 2 to 10 groups

cd "YOUR FOLDER"
clear all

do "ksmirnovk.do"

/** Reproduce lookup values from Table 1 in Kiefer 1959, and extend lookup table to more groups */
/** WARNING TAKES A LONG TIME TO RUN **/

clear matrix

* start with empty matrix
forvalues x = 0(0.01)5.00 {
	mat res = nullmat(res) \ `x'
}

* Collect results into matrix
forvalues k = 2(1)10 {
	local stop = 0
	forvalues x = 0(0.01)3.4 {
		if (`stop'==0) {
			mata: Phi = cdf_ksmirnovk(`x', `k'-1, theroots)
			mata: st_local("Phi",strofreal(Phi))
			mat res`k' = nullmat(res`k') \ `Phi'
			di "k=`k', x=`x' ::: func = `Phi'"	
			if (`Phi'>0.99999) {
				local stop = 1
			}
		}
		else {
			mat res`k' = nullmat(res`k') \ .
			di "k=`k', x=`x' ::: func = not computed"
		}
	}
	mat colnames res`k' = "k`k'"
	mat res = res , res`k'
	di "FINISHED k=`k'!!!!"
}
putexcel set "LookupTable.xlsx" , modify sheet("LookupTable")
putexcel A1=mat(res)



