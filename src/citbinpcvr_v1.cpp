#include <R.h>
#include <Rmath.h>
#include <vector>
#include <algorithm>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <iostream>
#include <random>       // std::default_random_engine
#include "logisticfunc.h"
#include <Rcpp.h>
#include "maxElementWithNan.h"

using namespace Rcpp;
using namespace std;

/*
L: matrix of continuous instrumental variables
G: matrix of candidate causal mediators
T: matrix of 0/1 variables
Programmer: Joshua Millstein
*/

// conduct permutations individually for each test so an intersection-union type approach can be applied to permutation-FDR

// [[Rcpp::export]]
void citbinpcvr( Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, 
    int &maxit, int &permit, int &boots, int &nrow, int &ncol, int &ncolc,
	Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, Rcpp::NumericVector pval3, Rcpp::NumericVector pval4, Rcpp::NumericVector pval3nc, Rcpp::IntegerVector perm_index, int &rseed)
{
	unsigned seed = rseed;
	int rw, brw, cl, i, j, rind, df, df1, df2, npos, nperm, dncol, perm, firstloop;
	int *bootind, *nposperm;
	double rss2, rss3, rss5, F, pv, pvalind, pvp, tmp, rhs;
	double *designmat, *phenovec, *pindep;
	bool aa, bb, cc, converged, permute;
	const int posno = 20;
	vector<vector<double> > LL;
	vector<vector<double> > CC;
	vector<double> gpred;
	vector<double> gresid;
	gsl_matrix *cov, *X;
	gsl_vector *Gm, *Gp, *c;

	bootind = new int[nrow];
	nposperm = new int[boots];
	designmat = new double[nrow * (ncol + ncolc + 2)];
	phenovec = new double[nrow];
	pindep = new double[boots];
	
	for(i = 0; i < boots; i++){
		nposperm[ i ] = 0;
	}	
	firstloop = permit;

	LL.resize( nrow );
	CC.resize( nrow );
	GetRNGstate();

	for(rw = 0; rw < nrow; rw++) {
		LL[rw].resize( ncol );	
		CC[rw].resize( ncolc );
	}

	for(cl = 0; cl < ncol; cl++) {
		for(rw = 0; rw < nrow; rw++) {
			LL[rw][cl] = L[rw + nrow * cl];
		}
	}
	for(cl = 0; cl < ncolc; cl++) {
		for(rw = 0; rw < nrow; rw++) {
			CC[rw][cl] = C[rw + nrow * cl];
		}
	}

	// begin permutation loop	
	for(perm = 0; perm < (boots + 1); perm++) {

		permute = perm > 0;
		for(rw = 0; rw < nrow; rw++) {
			bootind[ rw ]  = ( permute ) ? perm_index[ (perm - 1) * nrow + rw ] : rw;
		}

		// fit model T ~ C + L, to test T~L|C
		// create design matrix with no missing values
		dncol = 1 + ncolc + ncol;                               // intercept + multiple L variable
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ] ;

			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ brw ][ cl ]  != -9999 ) ? aa : 0;
			}
			for(cl = 0; cl < ncolc; cl++) {
                  aa = ( CC[ rw ][ cl ]  != -9999 ) ? aa : 0;
			}
			if( aa ){
				phenovec[ rind ] = T[ rw ];
				designmat[ rind * dncol  ] = 1;      // intercept
				for(cl = 0; cl < ncolc; cl++) {
                  designmat[ rind * dncol + 1 + cl  ]  = CC[ rw ][ cl ];
				}
				for(cl = 0; cl < ncol; cl++) {
                  designmat[ rind * dncol + ncolc + 1 + cl  ]  = LL[ brw ][ cl ];
				}
				rind++;
			} // end if aa
		} // end for rw

		df = ncol;
		converged = logisticReg( pv, phenovec, designmat, rind, dncol, df );
		if(!converged)Rcpp::warning("Cannot Converge when doing regression for calculating P-value.");
		pv = ( converged ) ? pv : std::numeric_limits<double>::quiet_NaN();
		pval1[perm] = pv;  // pval for T ~ L|C, 9 if it did not converge, p1

		// fit model T ~ C + L + G
		dncol = 1 + ncolc + ncol + 1;
		rind = 0;

		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ] ;
			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
			}
			for(cl = 0; cl < ncolc; cl++) {
                  aa = ( CC[ rw ][ cl ]  != -9999 ) ? aa : 0;
			}
			aa = ( G[ brw ] != -9999 ) ? aa : 0;
			if( aa ){
				phenovec[ rind ] = T[ rw ];
				designmat[ rind * dncol  ] = 1;      // intercept
				for(cl = 0; cl < ncolc; cl++) {
                  designmat[ rind * dncol + 1 + cl  ]  = CC[ rw ][ cl ];
				}
				for(cl = 0; cl < ncol; cl++) {
                  designmat[ rind * dncol + ncolc + 1 + cl  ]  = LL[ rw ][ cl ];
				}
				designmat[ rind * dncol + ncolc + 1 + ncol  ] = G[ brw ];
				rind++;
			} // end if aa
		} // end for rw
		df = 1;
		converged = logisticReg( pv, phenovec, designmat, rind, dncol, df );
		if(!converged)Rcpp::warning("Cannot Converge when doing regression for calculating P-value.");
		pv = ( converged ) ? pv : std::numeric_limits<double>::quiet_NaN();
		pval2[perm]  = pv;  // pval for T ~ G|L,C, 9 if it did not converge, p2
		
		// fit model G ~ T
		dncol = 2;
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ] ;
			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ brw ][ cl ]  != -9999 ) ? aa : 0;
			}
			aa = ( G[ rw ] != -9999 ) ? aa : 0;
			if( aa ){
				rind++;
			} // end if aa
		} // end for rw
		
		X = gsl_matrix_alloc(rind, dncol);
		Gm = gsl_vector_alloc (rind);
		
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ] ;
			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ brw ][ cl ]  != -9999 ) ? aa : 0;
			}
			aa = ( G[ rw ] != -9999 ) ? aa : 0;
			if( aa ){
				gsl_matrix_set(X, rind, 0, 1.0);      // intercept
				gsl_matrix_set(X, rind, 1, T[ rw ]); 
				gsl_vector_set(Gm, rind, G[rw]);
				rind++;
			} // end if aa
		} // end for rw
		
		c = gsl_vector_alloc(dncol);
		cov = gsl_matrix_alloc(dncol, dncol);
		gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(rind, dncol);
		gsl_multifit_linear(X, Gm, c, cov, &rss2, work);
		gsl_multifit_linear_free(work);
		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_vector_free(c);
	
		// fit model G ~ L + T
		dncol = 1 + ncol + 1;
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ] ;
			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ brw ][ cl ]  != -9999 ) ? aa : 0;
			}
			aa = ( G[ rw ] != -9999 ) ? aa : 0;
			if( aa ){
				rind++;
			} // end if aa
		} // end for rw
		
		X = gsl_matrix_alloc(rind, dncol);
		Gm = gsl_vector_alloc (rind);
		
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ] ;
			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ brw ][ cl ]  != -9999 ) ? aa : 0;
			}
			aa = ( G[ rw ] != -9999 ) ? aa : 0;
			if( aa ){
				gsl_matrix_set(X, rind, 0, 1.0);      // intercept
				for(cl = 0; cl < ncol; cl++) {
                  gsl_matrix_set(X, rind, cl + 1, LL[ brw ][ cl ] );
		     	}
				gsl_matrix_set(X, rind, 1 + ncol, T[ rw ]); 
				gsl_vector_set(Gm, rind, G[rw]);
				rind++;
			} // end if aa
		} // end for rw
		
		c = gsl_vector_alloc(dncol);
		cov = gsl_matrix_alloc(dncol, dncol);
		work = gsl_multifit_linear_alloc(rind, dncol);
		gsl_multifit_linear(X, Gm, c, cov, &rss3, work);
		gsl_multifit_linear_free(work);
		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_vector_free(c);
		df1 = ncol;
		df2 = rind - dncol;
		F = df2*(rss2-rss3)/(rss3*df1);
		pv = gsl_cdf_fdist_Q(F, df1, df2);
		pval3[perm]  = pv; // pval for G ~ L|T, p3
		
		// non-centrality parameter
		dncol = 1 + ncolc + 1 + ncol;
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
		  brw = bootind[ rw ] ;
		  aa = 1;
		  aa = ( T[ rw ] != -9999 ) ? aa : 0;
		  for(cl = 0; cl < ncol; cl++) {
		    aa = ( LL[ brw ][ cl ]  != -9999 ) ? aa : 0;
		  }
		  for(cl = 0; cl < ncolc; cl++) {
		    aa = ( CC[ rw ][ cl ]  != -9999 ) ? aa : 0;
		  }
		  aa = ( G[ rw ] != -9999 ) ? aa : 0;
		  if( aa ){
		    phenovec[ rind ] = T[ rw ];
		    designmat[ rind * dncol  ] = 1;      // intercept
		    for(cl = 0; cl < ncolc; cl++) {
		      designmat[ rind * dncol + 1 + cl  ]  = CC[ rw ][ cl ];
		    }
		    designmat[ rind * dncol + ncolc + 1  ] = G[ rw ];
		    for(cl = 0; cl < ncol; cl++) {
		      designmat[ rind * dncol + ncolc + 2 + cl  ]  = LL[ brw ][ cl ];
		    }
		    rind++;
		  } // end if aa
		} // end for rw
		
		df = ncol;
		converged = logisticReg( pv, phenovec, designmat, rind, dncol, df );
		if(!converged)Rcpp::warning("Cannot Converge when doing regression for calculating P-value.");
		pval3nc[perm] = ( converged ) ? pv : std::numeric_limits<double>::quiet_NaN();    // p-value for T ~ L|G,C, p-value for computing non-centrality parameter
		
		// fit model T ~ C + G + L
		if( perm == 0 ){
			dncol = 1 + ncolc + 1 + ncol;
			rind = 0;
			for(rw = 0; rw < nrow; rw++) {
				brw = bootind[ rw ] ;
				aa = 1;
				aa = ( T[ rw ] != -9999 ) ? aa : 0;
				for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				for(cl = 0; cl < ncolc; cl++) {
                  aa = ( CC[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				aa = ( G[ brw ] != -9999 ) ? aa : 0;
				if( aa ){
					phenovec[ rind ] = T[ rw ];
					designmat[ rind * dncol  ] = 1;      // intercept
					for(cl = 0; cl < ncolc; cl++) {
                  		designmat[ rind * dncol + 1 + cl  ]  = CC[ rw ][ cl ];
					}
					designmat[ rind * dncol + ncolc + 1  ] = G[ brw ];
					for(cl = 0; cl < ncol; cl++) {
                 		designmat[ rind * dncol + ncolc + 2 + cl  ]  = LL[ rw ][ cl ];
					}
					rind++;
				} // end if aa
			} // end for rw
		
			df = ncol;
			converged = logisticReg( pv, phenovec, designmat, rind, dncol, df );
			if(!converged)Rcpp::warning("Cannot Converge when doing regression for calculating P-value.");
			pvalind = ( converged ) ? pv : std::numeric_limits<double>::quiet_NaN();    // p-value for T ~ L|G,C
		
			// fit model G ~ L
			dncol = 1 + ncol;
			rind = 0;
			for(rw = 0; rw < nrow; rw++) {
				aa = 1;
				aa = ( G[ rw ] != -9999 ) ? aa : 0;
				for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				if( aa ){
					rind++;
				} // end if aa
			} // end for rw
		
			X = gsl_matrix_alloc(rind, dncol);
			Gm = gsl_vector_alloc (rind);
		
			rind = 0;
			for(rw = 0; rw < nrow; rw++) {
				aa = 1;
				aa = ( G[ rw ] != -9999 ) ? aa : 0;
				for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				if( aa ){
					gsl_matrix_set(X, rind, 0, 1.0);      // intercept
					for(cl = 0; cl < ncol; cl++) {
                 	gsl_matrix_set(X, rind, cl + 1, LL[ rw ][ cl ] );
		     		}
					gsl_vector_set(Gm, rind, G[rw]);
					rind++;
				} // end if aa
			} // end for rw
		
			c = gsl_vector_alloc(dncol);
			cov = gsl_matrix_alloc(dncol, dncol);
			work = gsl_multifit_linear_alloc(rind, dncol);
			gsl_multifit_linear(X, Gm, c, cov, &rss5, work);
			gsl_multifit_linear_free(work);
			gsl_matrix_free(cov);
			
			// residuals for G ~ L
			for(rw = 0; rw < nrow; rw++) {
				aa = 1;
				aa = ( G[ rw ] != -9999 ) ? aa : 0;
				for(cl = 0; cl < ncol; cl++) {
                 aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				if( aa ){
					rhs = gsl_vector_get(c, 0); // intercept
					for(cl = 0; cl < ncol; cl++) {
                  rhs += gsl_vector_get(c, (cl+1)) * LL[ rw ][ cl ];
		     		}
					gpred.push_back(rhs);
					tmp = gsl_vector_get(Gm, rw) - rhs;
					gresid.push_back(tmp);
				} 
				else {
					gpred.push_back(-9999);
					gresid.push_back(-9999);
				}
			} // end rw loop
			gsl_vector_free (c);
		} // end if perm == 0
			
		// Conduct an initial set of permutations
		if( perm > 0 ){
			// bootstrap the residuals like the other tests, but since the outcome is not bootstrapped, this test is conducted under the null.
			// compute G* based on marginal L effects and permuted residuals
			Gp = gsl_vector_alloc(nrow);
			for(rw = 0; rw < nrow; rw++) {
				brw = bootind[ rw ];
				aa = 1;
				aa = ( gpred[rw] != -9999 ) ? aa : 0;
				aa = ( gresid[brw] != -9999 ) ? aa : 0;
				if( aa ){
					gsl_vector_set(Gp, rw, gpred[rw] + gresid[brw] );
				}
				else {
					gsl_vector_set(Gp, rw, -9999 );
				}
			}
			
			// Recompute p-value for T ~ L|G,C based on G*
			// fit model T ~ C + G* + L to test L 
			dncol = 1 + ncolc + 1 + ncol;
			rind = 0;
			for(rw = 0; rw < nrow; rw++) {
				aa = 1;
				aa = ( T[ rw ] != -9999 ) ? aa : 0;
				for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				for(cl = 0; cl < ncolc; cl++) {
                  aa = ( CC[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				aa = ( gsl_vector_get(Gp, rw ) != -9999 ) ? aa : 0;
				if( aa ){
					phenovec[ rind ] = T[ rw ];
					designmat[ rind * dncol  ] = 1;      // intercept
					for(cl = 0; cl < ncolc; cl++) {
                  		designmat[ rind * dncol + 1 + cl  ]  = CC[ rw ][ cl ];
					}
					designmat[ rind * dncol + ncolc + 1  ] = gsl_vector_get(Gp, rw );
					for(cl = 0; cl < ncol; cl++) {
                  		designmat[ rind * dncol + ncolc + 2 + cl  ]  = LL[ rw ][ cl ];
					}
					rind++;
				} // end if aa
			} // end for rw
		
			df = ncol;
			converged = logisticReg( pvp, phenovec, designmat, rind, dncol, df );
			if(!converged)Rcpp::warning("Cannot Converge when doing regression for calculating P-value.");
			pindep[ perm - 1 ] = ( converged ) ? pvp : std::numeric_limits<double>::quiet_NaN();    // p-value for T ~ L|G*,C

		} // end if perm > 0
	} // End perm loop
		
	npos = 0;
	for(i = 0; i < firstloop; i++){
		// randomly permute residuals
		
		shuffle( gresid.begin(), gresid.end(), std::default_random_engine(seed) );

		for(rw = 0; rw < nrow; rw++) {
			brw = bootind[ rw ];
			aa = 1;
			aa = ( gpred[rw] != -9999 ) ? aa : 0;
			aa = ( gresid[brw] != -9999 ) ? aa : 0;
			if( aa ){
				gsl_vector_set(Gp, rw, gpred[rw] + gresid[rw] );
			}
			else {
				gsl_vector_set(Gp, rw, -9999 );
			}
		} // end rw loop
		// fit model T ~ C + G* + L to test L 
		dncol = 1 + ncolc + 1 + ncol;
		rind = 0;
		for(rw = 0; rw < nrow; rw++) {
			aa = 1;
			aa = ( T[ rw ] != -9999 ) ? aa : 0;
			for(cl = 0; cl < ncol; cl++) {
          		aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
			}
			for(cl = 0; cl < ncolc; cl++) {
           	aa = ( CC[ rw ][ cl ]  != -9999 ) ? aa : 0;
			}
			aa = ( gsl_vector_get(Gp, rw ) != -9999 ) ? aa : 0;
			if( aa ){
				phenovec[ rind ] = T[ rw ];
				designmat[ rind * dncol  ] = 1;      // intercept
				for(cl = 0; cl < ncolc; cl++) {
                  	designmat[ rind * dncol + 1 + cl  ]  = CC[ rw ][ cl ];
				}
				designmat[ rind * dncol + ncolc + 1  ] = gsl_vector_get(Gp, rw );
				for(cl = 0; cl < ncol; cl++) {
                  	designmat[ rind * dncol + ncolc + 2 + cl  ]  = LL[ rw ][ cl ];
				}
				rind++;
			} // end if aa
		} // end for rw
			
		df = ncol;
		converged = logisticReg( pvp, phenovec, designmat, rind, dncol, df );
		if(!converged)Rcpp::warning("Cannot Converge when doing regression for calculating P-value.");
		pvp = ( converged ) ? pvp : std::numeric_limits<double>::quiet_NaN();    // p-value for T ~ L|G*,C
		if( pvp > pvalind ) npos++;
			
		for(j = 0; j < boots; j++){
			if( pvp > pindep[ j ] ) nposperm[ j ] += 1;
		}
	} // end initial i permutation loop
		
	nperm = firstloop;
			
	// Conduct additional permutations if there is some indication of statistical significance
	if( boots == 0 ){
		aa = npos < posno;
		bb = nperm < maxit;
				
		while(aa && bb) {
	
			// randomly permute residuals
            
            shuffle( gresid.begin(), gresid.end(), std::default_random_engine(seed) );

			for(rw = 0; rw < nrow; rw++) {
				brw = bootind[ rw ];
				aa = 1;
				aa = ( gpred[rw] != -9999 ) ? aa : 0;
				aa = ( gresid[brw] != -9999 ) ? aa : 0;
				if( aa ){
					gsl_vector_set(Gp, rw, gpred[rw] + gresid[rw] );
				}
				else {
					gsl_vector_set(Gp, rw, -9999 );
				}
			} // end rw loop
			// fit model T ~ C + G* + L to test L 
			dncol = 1 + ncolc + 1 + ncol;
			rind = 0;
			for(rw = 0; rw < nrow; rw++) {
				aa = 1;
				aa = ( T[ rw ] != -9999 ) ? aa : 0;
				for(cl = 0; cl < ncol; cl++) {
                  aa = ( LL[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				for(cl = 0; cl < ncolc; cl++) {
           		aa = ( CC[ rw ][ cl ]  != -9999 ) ? aa : 0;
				}
				aa = ( gsl_vector_get(Gp, rw ) != -9999 ) ? aa : 0;
				if( aa ){
					phenovec[ rind ] = T[ rw ];
					designmat[ rind * dncol  ] = 1;      // intercept
					for(cl = 0; cl < ncolc; cl++) {
                  		designmat[ rind * dncol + 1 + cl  ]  = CC[ rw ][ cl ];
					}
					designmat[ rind * dncol + ncolc + 1  ] = gsl_vector_get(Gp, rw );
					for(cl = 0; cl < ncol; cl++) {
                  		designmat[ rind * dncol + ncolc + 2 + cl  ]  = LL[ rw ][ cl ];
					}
					rind++;
				} // end if aa
			} // end for rw
			
			df = ncol;
			converged = logisticReg( pvp, phenovec, designmat, rind, dncol, df );
			if(!converged)Rcpp::warning("Cannot Converge when doing regression for calculating P-value.");
			pvp = ( converged ) ? pvp : std::numeric_limits<double>::quiet_NaN();    // p-value for T ~ L|G*,C
			if( pvp > pvalind ) npos++;
			
			for(j = 0; j < boots; j++){
				if( pvp > pindep[ j ] ) nposperm[ j ] += 1;
			}
				
			aa = npos < posno;
			cc = nperm < ( maxit - 1 );
			nperm++;
		} // end 'while' permutation loop
	} // end if boots == 0 
	
	for(perm = 0; perm < (boots + 1); perm++) {
		pv = ( perm == 0 ) ? 1.0 * npos / nperm : 1.0 * nposperm[ perm - 1 ] / nperm;
		
		// To avoid a p-value = 0, make 0 p-values = 1/nperm
		pval4[perm] = ( pv > 0 ) ? pv : 1.0 * 1.0 / nperm;   // pval for L ind T|G,C
		
	} // End perm loop

	gresid.clear();
	gpred.clear();
	LL.clear();
	gsl_vector_free (Gm);
	gsl_vector_free (Gp);
	PutRNGstate();
	
	delete [] bootind;
	delete [] nposperm;
	delete [] designmat;
	delete [] phenovec;
	delete [] pindep;

} // End citconlog3pcvr