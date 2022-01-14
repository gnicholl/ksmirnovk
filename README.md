# ksmirnovk
Implementation of k-sample version Kolmogorov-Smirnov test in Stata. Read more discussion on my blog post: https://gnicholl.github.io/post/ksmirnovk/

# Required package
ssc install integrate.pkg

# Description of program
* syntax: ksmirnovk y , by(x)
  * y: variable of interest
  * x: categorical variable specifying k groups numbered 1,2,...,k

* ksmirnovk produces a p-value for testing the null hypothesis
  that the distribution of y is the same in all k groups of x
  * It is a non-parametric test, and is the k-sample analogue of the
    ksmirnov command in Stata. 
  * For k=2 groups, ksmirnovk is theoretically identical to ksmirnov, though I include more terms in
    the approximation of the p-value, so output may differ slightly
  * As a side effect, it also produces a plot of the empirical cdfs of each group

* I compute the test based on Kiefer (1959)
  * the test statistic requires the computation of each group's empirical cdf,
    as well as the empirical cdf of the pooled data. The statistic is based on
    a sum of max (squared) differences between the group and pooled cdfs.
  * To compute a p-value, we need to calculate the cdf of the test statistic
    * Given k groups,
	    * if k=2 (standard ksmirnov case) or k=4, there is an easy, closed-form expression
	    * otherwise, it is a numercially intensive procedure.
	      *  the function is an infinite series, and each term depends on finding
	         the root of a bessel function, and evaluating a different bessel function at that root
	      * bessel functions themselves do not have a closed form solution, and instead can
          be represented by an infinite sum or a difference of two integrals. Here I
          use the integral representation (see Temme 1996) and use 200 quadrature points
	  * To save time, I pre-compute the first 30 roots of the bessel function for cases
	    k = 1,3,5,6,7,8,9,10. I found that 30 roots was enough to reproduce
      Kiefer (1959)'s results (which go to 6 decimal places) (see LookupTable.xlsx). 
	  * for k>10, the program has to compute the bessel roots on the fly--not recommended,
      it is slow and has NOT been tested

# Sources
* K-Sample Analogues of the Kolmogorov-Smirnov and Cramer-V. Mises Test
  * by J. Kiefer
  *  The Annals of Mathematical Statistics, Jun 1959, Vol 30, No 2 pp. 420-447
  *  https://www.jstor.org/stable/2237091
* Special Functions: An introduction to the classical functions of mathematical physics 
  *  by Temme, Nico M. (1996)  pp. 228-231
* chebroots adapted from matlab code here: https://github.com/chebfun/chebfun
