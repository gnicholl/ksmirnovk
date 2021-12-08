* Program to run K-sample Kolmogorov-Smirnov test

//// SPECIFY DIRECTORY /////////////////////////////////////////////////////////
global mydir "YOUR FOLDER"
////////////////////////////////////////////////////////////////////////////////

//// INSTALL PACKAGES //////////////////////////////////////////////////////////
*ssc install integrate.pkg
////////////////////////////////////////////////////////////////////////////////

//// SOURCES ///////////////////////////////////////////////////////////////////
* K-Sample Analogues of the Kolmogorov-Smirnov and Cramer-V. Mises Tests
*  by J. Kiefer
*  The Annals of Mathematical Statistics, Jun 1959, Vol 30, No 2 pp. 420-447
*  https://www.jstor.org/stable/2237091
* Special Functions: An introduction to the classical functions of mathematical physics 
*  by Temme, Nico M. (1996)  pp. 228-231
* chebroots adapted from matlab code here: https://github.com/chebfun/chebfun
////////////////////////////////////////////////////////////////////////////////

//// DESCRIPTION ///////////////////////////////////////////////////////////////
* syntax: ksmirnovk y , by(x)
  * y: variable of interest
  * x: categorical variable specifying k groups numbered 1,2,...,k

* ksmirnovk produces a p-value for testing the null hypothesis
 * that the distribution of y is the same in all k groups of x

* It is a non-parametric test, and is the k-sample analogue of the
 * ksmirnov command in Stata. for k=2 groups, ksmirnovk is theoretically
 * identical to ksmirnov, though I include more terms in
 * the approximation of the p-value, so output may differ slightly

* As a side effect, it also produces a plot of the empirical cdfs of each group

* I compute the test based on Kiefer (1959)
  * - the test statistic requires the computation of each group's empirical cdf,
  *   as well as the empirical cdf of the pooled data. The statistic is based on
  *   a sum of max (squared) differences between the group and pooled cdfs.
  * - To compute a p-value, we need to calculate the cdf of the test statistic
      * - Given k groups...
	  * - if k=2 (standard ksmirnov case) or k=4, there is an easy, closed-form expression
	  * - otherwise, it is a numercially intensive procedure.
	  *      ->the function is an infinite series, and each term depends on finding
	  *        the root of a bessel function, and evaluating a different bessel function at that root
	  *      ->bessel functions themselves do not have a closed form solution, and instead can
      *        be represented by an infinite sum or a difference of two integrals. Here I
      *        use the integral representation (see Temme 1996) and use 200 quadrature points
	  * - To save time, I pre-compute the first 30 roots of the bessel function for cases
	  *    k = 1,3,5,6,7,8,9,10. I found that 30 roots was enough to reproduce
      *    Kiefer (1959)'s results (which go to 6 decimal places). 
	  * - for k>10, the program has to compute the bessel roots on the fly--not recommended,
      *     it is slow and has NOT been tested

* Notes:
  * - I could likely speed up the program by reducing the number of roots used
  *   and/or reducing the number of qudrature draws without sacrificing much
  *   (if any) precision. Perhaps one day I will be motivated to experiment with that...
  * - I should probably use temporary variables to avoid potential conflicts
  *   with already-existing variables
////////////////////////////////////////////////////////////////////////////////

* Read in bessel roots generated from genroots.do
mata:
rootsdata = xl()
rootsdata.load_book("$mydir\besselroots.xlsx")
rootsdata.set_sheet("besselroots")
theroots = rootsdata.get_number((2,31),(1,7))
rootsdata.close_book()
end

capture program drop ksmirnovk
capture mata mata drop ecdf_groups()
capture mata mata drop myecdf_grp()
capture mata mata drop computeecdfs()
capture mata mata drop besselmat()
capture mata mata drop offset_diag()
capture mata mata drop chebroots_rec()
capture mata mata drop chebroots()
capture mata mata drop besselJ_int1()
capture mata mata drop besselJ_int2()
capture mata mata drop besselJ()
capture mata mata drop besselK_int()
capture mata mata drop besselK()
capture mata mata drop cdf_ksmirnovk()
capture mata mata drop nroots()


program define ksmirnovk, eclass
	syntax varlist , BY(string)
	
	
	* Calculate number of groups in by-variable
	quietly count if `by' == .
	local miss = r(N)
	tempvar nvals
	by `by', sort: gen `nvals' = (_n==1)
	quietly count if `nvals'
	local ngroups = r(N)
	if `miss' > 0 {
		local ngroups = `ngroups' - 1
	}
	
	* Sort values and funnel missings to bottom
	qui gen new_`varlist' = `varlist'
	qui gen new_`by' = `by'
	qui replace new_`varlist' = . if `by'==.
	qui replace new_`by' = . if `varlist'==.
	sort new_`varlist'
	drop new_`varlist' new_`by'
	
	* Generate ecdfs by group
	forvalues g = 1(1)`ngroups' {
		capture drop ecdf_`by'_`g'
		qui gen ecdf_`by'_`g' = .
		
		capture drop n_`by'_`g'
		qui gen n_`by'_`g' = .
	}
	capture drop ecdf_overall
	qui gen ecdf_overall = .
	mata: computeecdfs("`varlist'", "`by'",`ngroups')
	
	* Calculate test statistic for each x
	local genstring = ""
	
	forvalues g = 1(1)`ngroups' {
		if ("`g'" == "1") {
			local genstring = "n_`by'_1*(ecdf_`by'_1 - ecdf_overall)^2"
		}
		else {
			local genstring = "`genstring' + n_`by'_`g'*(ecdf_`by'_`g' - ecdf_overall)^2"
		}
	}
	capture drop ksktest
	qui gen ksktest = `genstring'
	
	* Calculate test statistics
	quietly summ ksktest
	local maxksk = r(max)
	local meanksk = r(mean)
	local maxrootT = sqrt(`maxksk')
	
	* Compute p-values with mata, store in local macros
	local pvalue = 0
	local pvaluecvm = 0
	qui mata: cdf_ksmirnovk(`maxrootT',`ngroups'-1,theroots)
	
	* Output
	display "-----"
	display "Number of groups = `ngroups'"
	display "-----"
	display "KS = `maxrootT'"
	display "Pvalue (KS) = `pvalue'"
	display "-----"
	
	* Produce fancy graph
	local ecdflist = ""
	forvalue i = 1(1)`ngroups' {
		local ecdflist = "`ecdflist' ecdf_`by'_`i'"
	}
	line `ecdflist' `varlist', sort
	
	* Set post-estimation values
	ereturn scalar KS = `maxrootT'
	ereturn scalar pKS = `pvalue'
	
end

mata:
	
	// Calculate ecdf for each group sample and for pooled sample (output as row vec)
	// Assumes data is sorted so it can stop once it reaches high enough value
	// x - value we are evaluating ecdf at
	// mydata - Nx2 matrix, contains data in first col, group in second col
	// grouptotals - row vector containing sample sizes, length 
	function ecdf_groups(numeric scalar x, mydata, grouptotals) {
	
		// Number of groups
		ng = cols(grouptotals)
		
		// Sample sizes in each group, last 
		ns = grouptotals, rowsum(grouptotals)
		
		// Will contain number of values <= x
		belowx = J(1,ng+1,0)
		
		i = 1
		while (mydata[i,1] <= x) {
			g = mydata[i,2]
			belowx[g] = belowx[g] + 1
			belowx[ng+1] = belowx[ng+1] + 1
			i++
			if (i > rows(mydata)){
				break
			}
		} 
		
		return(belowx:/ns)
	}
	
	// Compute ecdfs from stata data and output stata variables
	function computeecdfs(string scalar studyvar, string scalar groupvar, numeric scalar ngroups) {
	
		// Copy data from stata
		data = st_data(., (studyvar,groupvar))
		
		// Sort data by study variable (already done in stata part)
		//data = sort(data,1)
		
		// Trim off missing values (missing values are sent to bottom from sorting)
		nomiss = (rows(data)-sum(rowmissing(data):>0))
		data = data[1..nomiss,]
	
		// Calculate sample sizes in each group
		gtotals = J(1,ngroups,0)
		for (i=1; i<= rows(data); i++) {
			g = data[i,2]
			gtotals[g] = gtotals[g]+1
		}
	
		// Calculate ecdf's for every x
		ecdfs = ecdf_groups(data[1,1], data, gtotals)
		for (i=2; i<=rows(data); i++) {
			ecdfs = ecdfs \ ecdf_groups(data[i,1],data,gtotals)
		}
	
		// Create vectors with sample sizes
		gtotals = J(rows(data),1,gtotals)
	
		// Generate stata variables
		for (grp=1; grp<=ngroups; grp++) {
			st_store(1::nomiss, "ecdf"+"_"+groupvar+"_"+strofreal(grp), ecdfs[,grp])
			st_store(1::nomiss, "n"+"_"+groupvar+"_"+strofreal(grp), gtotals[,grp])
		}
		st_store(1::nomiss, "ecdf_overall", ecdfs[,ngroups+1])
	
	}
	
	// Calculate bessel function from matrix to matrix
	function besselmat(in,params) {
		out = in
		for (i=1; i<=rows(in); i++) {
			for (j=1; j<=cols(in); j++) {
				out[i,j] = besselJ(in[i,j],params)
			}
		}
		return(out)
	}
	
	// Function to create a matrix with diagonal elements
	// formed from vector v. off is an offset number (+ve
	// or -ve integer) which offsets the diagonal up or down
	function offset_diag (v, off) {
		dim = length(v) + abs(off)
		result = J(dim,dim, 0)
		if (off>=0) {
			for (i=1; i<=length(v); i++) {
				result[i,i+off] = v[i]
			}
		}
		else {
			for (i=1; i<=length(v); i++) {
				result[i-off,i] = v[i]
			}
		}
		return(result)
	}
	
   function chebroots_rec ( f , a , b, N, x, V, params ) {

        // get the centre and half-width of the interval [a,b]
        m = (a + b) / 2
        h = (b - a) / 2

        // Compute the coefficients
        fx = (*f)( m :+ h:*x, params )
		help = (*f)(x, params)
		nf = norm(help,.)
        c = qrsolve(V, fx)

        // Trim off residual trailing coefficients
        n = N+1
        while (abs(c[n]) < (2^-52)*100*nf) {
			n = n - 1
		}
      
		if (n <= 2) { // Trivial case		
			r = - c[1] / c[2]
			return(r)
		}
        else if (N - n < 3) { // Recurse?
			rleft = chebroots_rec( f , a , m, N, x, V, params )
			rright = chebroots_rec( f , m , b, N, x, V, params )
			r = rleft \ rright
			return(r)
		}
        else { // Otherwise not
            // normalize the coefficients
            c = c[1..n] / c[n]

            // Construct the colleague matrix M
            gamma = J(n-2,1,1) / 2
            M = offset_diag(gamma,1) + offset_diag(gamma,-1)
			for (i=1; i<=n-1; i++) {
				M[n-1,i] = M[n-1,i] - 0.5*c[i]
			}
			M[1,2] = 1

            // Get the eigenvalues of M
            lambda = eigenvalues(M)'
			
            // reduce to values with negligible imaginary part
            r = select(lambda, abs(Im(lambda)) :< 10*(2^-52)*abs(lambda))

            // reduce to roots inside the interval
            r = Re(select(r, abs(r) :<= 1+1e4*(2^-52) ))
			r = rowmin( (rowmax( (r,J(length(r),1,-1)) ), J(length(r),1,1))  )
            // map back to the interval [a,b]
			r = m :+ h:*r
            return(r)

        } 
    }

	function chebroots( f , a , b , N, params ) {
		// CHEBROOTS approximates the roots of f in [a,b].
		//   R=CHEBROOTS(F,A,B,N) approximates the real roots R of F(X)
		//   in the interval [A,B] by computing a Chebyshev interpolation
		//   of degree at least N and finding its roots using a
		//   colleague matrix.

    // Compute the nodes and the Vandermonde-like matrix
	x = cos( pi() * (0::N)/N )
	V = J(N+1,1,1), x
    for (k=3; k<=N+1; k++) {
		V = V, 2*x:*V[1..N+1,k-1] - V[1..N+1,k-2]
	}
	
    // Run helper
	return(chebroots_rec ( f , a , b, N, x, V, params ))
	}
	
	/////// BESSEL FUNCTION ////////////////////////////////////////////////////
	
	// Compute bessel J function using quadrature (currently always use 200 points)
      // formula comes from Temme (1996)
	function besselJ_int1(t, params) {
		return(  cos(params[2]*t - params[1]*sin(t))  )
	}
	function besselJ_int2(t, params) {
		return(  exp(-params[1]*sinh(t) - params[2]*t)  )
	}
	function besselJ (z,v) {
		z = abs(z)
		return(   ((1/pi()):*integrate(&besselJ_int1(), 0, pi(),200,(z,v))) - ((sin(v*pi())/pi()):*(integrate(&besselJ_int2(),0,.,200,(z,v))))    )
	}
	
	////////////////////////////////////////////////////////////////////////////
	
	
	/////// Asymptotic cdfs ////////////////////////////////////////////////////
	
	// Compute asymptotic ksmirnov cdf given h (= num of groups-1)
	// Based on results by Kiefer (1959)
	function cdf_ksmirnovk(x, h, theroots) {
		if (x<=0) { // Trivial-> outside of range has 0 probability
			result = 0
		}
		else if (x>=5) { // If >=5, result = 1...true for h up to 9 => i.e. 10-group comparison
			result = 1
		}
		else {
			if (h == 1) { // Easier formula from Keifer
				outside = ((2*pi())^(1/2))/x
				runsum = 0
				for (i=1; i<=40; i++) {
					runsum = runsum + exp(   -((2*i-1)^2)*pi()*pi()/(8*x^2)  )
				}
				result = outside*runsum
			}
			else if (h==3) { // Easier formula from Keifer
				outside = ((2*pi()^5)^(1/2))/(x^3)
				runsum = 0
				for (i=1; i<=40; i++) {
					runsum = runsum + (i^2)*exp(   -(i^2)*pi()*pi()/(2*x^2)  )
				}
				result = outside*runsum
			}
			else if (h<=9) { // Use pre-calculated roots
				if (h==2) {
					gammas = theroots[,1]
				}
				else {
					gammas = theroots[,h-2]
				}

				outside = 4/(gamma(h/2)*(2^(h/2))*(x^h))
				runsum = 0
				for (i=1; i<=length(gammas); i++) {
					runsum = runsum + (  ((gammas[i]^(h-2)) * exp(-(gammas[i]^2)/(2*x^2))) / (besselJ(gammas[i],h/2))^2  )
				}
				result = outside*runsum
			}
			else { // Calculate roots on the spot
				gammas = nroots(&besselmat(),(h-2)/2,10)
				outside = 4/(gamma(h/2)*(2^(h/2))*(x^h))
				runsum = 0
				for (i=1; i<=length(gammas); i++) {
					runsum = runsum + (  ((gammas[i]^(h-2)) * exp(-(gammas[i]^2)/(2*x^2))) / (besselJ(gammas[i],h/2))^2  )
				}
				result = outside*runsum
			}
		}
		st_local("pvalue",strofreal(1-result))
		return(result)
	}
	
	////////////////////////////////////////////////////////////////////////////
	
	// Get first n roots of a function by incrementing by 20
	function nroots(f,params,n) {
		roots = chebroots( f , 0.3 , 20.3 , 100, params )
		startpoint = 20.3
		while(length(roots) < n) {
			roots = roots \ chebroots( f , startpoint , startpoint+20 , 100, params )
			startpoint = startpoint+20
		}
		if ( abs((*f)(0,params)) < 1e-9) {
			roots = 0 \ roots
		}
		return(sort(roots,1))
	}
end