
* Generate first 30 roots of the bessel function with
*   v = 0, 1, 1.5, 2, 2.5, 3, 3.5
* and save them to specified csv file
* this code takes a while to run!! which is why I have pre-computed them rather
*  then calculate them on the fly in ksmirnovk program

* Easier formulas are provided by Keifer (1959) in the case of 2 groups 4 groups (h=1 or h=3)
* Thus, we need only pre-compute bessel roots in the case of 3 or 5+ groups
* The resulting csv file contains roots for k = 3,5,6,7,8,9,10 respectively
*   (In Keifer's paper, this corresponds to h = 2,4,5,6,7,8,9)

cd "YOUR FOLDER"

global csvfile "besselroots.csv"

* Need to install numerical integration package for mata
	* run 'ssc install integrate.pkg'

clear all

set obs 30
gen roots0 = .
gen roots1 = .
gen roots1_5 = .
gen roots2 = .
gen roots2_5 = .
gen roots3 = .
gen roots3_5 = .

mata:

// Compute bessel J function using quadrature (currently always use 200 points)
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

// helper function for chebroots, see below
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

x = nroots(&besselmat(), 0, 30)
st_store(., "roots0", x[1::30])

x = nroots(&besselmat(), 1, 30)
st_store(., "roots1", x[1::30])

x = nroots(&besselmat(), 1.5, 30)
st_store(., "roots1_5", x[1::30])

x = nroots(&besselmat(), 2, 30)
st_store(., "roots2", x[1::30])

x = nroots(&besselmat(), 2.5, 30)
st_store(., "roots2_5", x[1::30])

x = nroots(&besselmat(), 3, 30)
st_store(., "roots3", x[1::30])

x = nroots(&besselmat(), 3.5, 30)
st_store(., "roots3_5", x[1::30])

end

export delimited using "$csvfile"
