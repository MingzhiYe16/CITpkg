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
#include <random> // std::default_random_engine
#include "linearfunc.h"
#include "maxElementWithNan.h"

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

/*
C: matrix of continuous adjustment covariates for T
L: matrix of continuous instrumental variables
G: matrix of candidate causal mediators
T: matrix of 0/1 variables
CG: matrix of continuous covariates for T
Programmer: Joshua Millstein
*/

// [[Rcpp::export]]
void citconlog2cvr_linear(Rcpp::NumericVector L, Rcpp::NumericVector G, Rcpp::NumericVector T, Rcpp::NumericVector C, 
                          Rcpp::NumericVector CG, int &nrow, int &ncol, int &ncolc, int &ncolct,
                          Rcpp::NumericVector pval, Rcpp::NumericVector pval1, Rcpp::NumericVector pval2, 
                          Rcpp::NumericVector pval3, Rcpp::NumericVector pval4, int &maxit, int &rseed)
{
    unsigned seed = rseed;
    int rw, cl, i, rind, df, nobs, ip, npos, nperm, nmiss, stride;
    double rss5, tmp, rhs, maxp, testval;
    double pv = 9.0;  // Need to initialize value, otherwise conditional jump
    double pvp = 9.0; // Need to initialize value, otherwise conditional jump
    double *designmat, *phenovec;
    bool aa, bb, cc, dd, converged;
    const int firstloop = 2000;
    const int posno = 20;
    const double alpha = .1;
    vector<vector<double>> LL;
    vector<vector<double>> CC;
    vector<vector<double>> CGG;
    vector<double> pvec;
    vector<double> gpred;
    vector<double> gresid;
    gsl_matrix *Lm, *cov, *X, *Cm, *CGm;
    gsl_vector *Gm, *Tm, *Gp, *c;

    designmat = new double[nrow * (ncol + ncolc + ncolct + 2)];
    phenovec = new double[nrow];

    LL.resize(nrow);
    CC.resize(nrow);
    CGG.resize(nrow);

    GetRNGstate();

    for (rw = 0; rw < nrow; rw++)
    {
        LL[rw].resize(ncol);
        CC[rw].resize(ncolc);
        CGG[rw].resize(ncolct);
    }

    for (cl = 0; cl < ncol; cl++)
    {
        for (rw = 0; rw < nrow; rw++)
        {
            LL[rw][cl] = L[rw + nrow * cl];
        }
    }
    for (cl = 0; cl < ncolc; cl++)
    {
        for (rw = 0; rw < nrow; rw++)
        {
            CC[rw][cl] = C[rw + nrow * cl];
        }
    }
    for (cl = 0; cl < ncolct; cl++)
    {
        for (rw = 0; rw < nrow; rw++)
        {
            CGG[rw][cl] = CG[rw + nrow * cl];
        }
    }

    // create analysis vectors w/no missing data
    nobs = 0;
    for (rw = 0; rw < nrow; rw++)
    {
        nmiss = 0;
        for (cl = 0; cl < ncol; cl++)
        {
            if (LL[rw][cl] == -9999)
            {
                nmiss++;
            }
        }
        dd = 1;
        for (cl = 0; cl < ncolc; cl++)
        {
            dd = (CC[rw][cl] != -9999) ? dd : 0;
        }
        for (cl = 0; cl < ncolct; cl++)
        {
            dd = (CGG[rw][cl] != -9999) ? dd : 0;
        }
        aa = nmiss == 0;
        bb = G[rw] != -9999;
        cc = T[rw] != -9999;
        if (aa && bb && cc && dd)
        {
            nobs++;
        }
    } // End for rw

    if (ncolc > 0)
	{
        Cm = gsl_matrix_alloc(nobs, ncolc);
    }
    if (ncolct > 0)
    {
        CGm = gsl_matrix_alloc(nobs, ncolct);
    }
    Lm = gsl_matrix_alloc(nobs, ncol);
    Gm = gsl_vector_alloc(nobs);
    Tm = gsl_vector_alloc(nobs);
    rind = 0;
    for (rw = 0; rw < nrow; rw++)
    {
        nmiss = 0;
        for (cl = 0; cl < ncol; cl++)
        {
            if (LL[rw][cl] == -9999)
            {
                nmiss++;
            }
        }
        aa = nmiss == 0;
        bb = G[rw] != -9999;
        cc = T[rw] != -9999;
        dd = 1;
        for (cl = 0; cl < ncolc; cl++)
        {
            dd = (CC[rw][cl] != -9999) ? dd : 0;
        }
        for (cl = 0; cl < ncolct; cl++)
        {
            dd = (CGG[rw][cl] != -9999) ? dd : 0;
        }

        if (aa && bb && cc && dd)
        {
            for (cl = 0; cl < ncol; cl++)
            {
                gsl_matrix_set(Lm, rind, cl, LL[rw][cl]);
            }
            for (cl = 0; cl < ncolc; cl++)
            {
                gsl_matrix_set(Cm, rind, cl, CC[rw][cl]);
            }
            for (cl = 0; cl < ncolct; cl++)
            {
                gsl_matrix_set(CGm, rind, cl, CGG[rw][cl]);
            }
            gsl_vector_set(Gm, rind, G[rw]);
            gsl_vector_set(Tm, rind, T[rw]);
            rind++;
        }
    }

	// fit model T ~ C + L
	ip = 1 + ncolc + ncol; // intercept + covariates + multiple L variable
	for (rw = 0; rw < nobs; rw++)
	{
		phenovec[rw] = gsl_vector_get(Tm, rw);
		designmat[rw * ip] = 1; // intercept
		for (cl = 0; cl < ncolc; cl++)
		{
			designmat[rw * ip + 1 + cl] = gsl_matrix_get(Cm, rw, cl);
		}
		for (cl = 0; cl < ncol; cl++)
		{
			designmat[rw * ip + 1 + ncolc + cl] = gsl_matrix_get(Lm, rw, cl);
		}
	}
	df = ncol;
	converged = linearRegCompare(pv, phenovec, designmat, nobs, ip, df);
	if (!converged)
		Rcpp::Rcout << "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
	pv = (converged) ? pv : std::numeric_limits<double>::quiet_NaN();
	pvec.push_back(pv); // pval for T ~ C + L, 9 if it did not converge, p1

	// fit model T ~ C + L + G
	stride = ip + 1;
	for (rw = 0; rw < nobs; rw++)
	{
		designmat[rw * stride] = 1; // intercept
		for (cl = 0; cl < ncolc; cl++)
		{
			designmat[rw * stride + 1 + cl] = gsl_matrix_get(Cm, rw, cl);
		}
		for (cl = 0; cl < ncol; cl++)
		{
			designmat[rw * stride + 1 + ncolc + cl] = gsl_matrix_get(Lm, rw, cl);
		}
		designmat[rw * stride + 1 + ncolc + ncol] = gsl_vector_get(Gm, rw);
	}

	df = 1;
	converged = linearRegCompare(pv, phenovec, designmat, nobs, stride, df);
	if (!converged)
		Rcpp::Rcout << "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
	pv = (converged) ? pv : std::numeric_limits<double>::quiet_NaN();
	pvec.push_back(pv); // pval for T ~ G|L,C, 9 if it did not converge, p2

	// fit model G ~ L + T + CG, L could be reduced
	gsl_multifit_linear_workspace *work;
	rind=0;
	int dncol=ncolct + ncol + 1 + 1; // CG + L + T + Intercept
	for (rw = 0; rw < nrow; rw++)
	{
		aa = 1;
		aa = (G[rw] != -9999) ? aa : 0;
		for (cl = 0; cl < ncol; cl++)
		{
			aa = (LL[rw][cl] != -9999) ? aa : 0;
		}
		for (cl = 0; cl < ncolct; cl++)
		{
			aa = (CGG[rw][cl] != -9999) ? aa : 0;
		}
		if (aa)
		{
			phenovec[rind] = G[rw];
			designmat[rind * dncol] = 1; // intercept
			designmat[rind * dncol + 1] = T[rw]; // T
			for (cl = 0; cl < ncolct; cl++) // CG
			{
				designmat[rind * dncol + 2 + cl] = CGG[rw][cl];
			}
			for (cl = 0; cl < ncol; cl++)
			{
				designmat[rind * dncol + ncolct + 2 + cl] = LL[rw][cl];
			}
			rind++;
		} // end if aa
	} // end for rw

	df=ncol;
	converged = linearRegCompare(pv, phenovec, designmat, rind, dncol, df);
	if (!converged)
		Rcpp::Rcout << "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
	pv = (converged) ? pv : std::numeric_limits<double>::quiet_NaN();
	pvec.push_back(pv); // pval for G ~ L + T + CG, L could be reduced, p3



	// fit model T ~ C + G + L to test L
	for (rw = 0; rw < nobs; rw++)
	{
		designmat[rw * stride] = 1; // intercept
		for (cl = 0; cl < ncolc; cl++)
		{
			designmat[rw * stride + 1 + cl] = gsl_matrix_get(Cm, rw, cl);
		}
		designmat[rw * stride + ncolc + 1] = gsl_vector_get(Gm, rw);
		for (cl = 0; cl < ncol; cl++)
		{
			designmat[rw * stride + ncolc + 2 + cl] = gsl_matrix_get(Lm, rw, cl);
		}
	}

	df = ncol;
	converged = linearRegCompare(pv, phenovec, designmat, nobs, stride, df);
	if (!converged)
		Rcpp::Rcout << "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
	pv = (converged) ? pv : std::numeric_limits<double>::quiet_NaN(); // p-value for T ~ L|G + C

	// fit model G ~ L
	X = gsl_matrix_alloc(nobs, ncol + 1);
	for (rw = 0; rw < nobs; rw++)
	{
		gsl_matrix_set(X, rw, 0, 1.0); // intercept
		for (cl = 0; cl < ncol; cl++)
		{
			gsl_matrix_set(X, rw, cl + 1, gsl_matrix_get(Lm, rw, cl));
		}
	}
	c = gsl_vector_alloc(ncol + 1);
	cov = gsl_matrix_alloc(ncol + 1, ncol + 1);
	work = gsl_multifit_linear_alloc(nobs, ncol + 1);
	gsl_multifit_linear(X, Gm, c, cov, &rss5, work);
	gsl_multifit_linear_free(work);
	gsl_matrix_free(cov);

	// residuals for G ~ L
	for (rw = 0; rw < nobs; rw++)
	{
		rhs = 0;
		for (cl = 0; cl < (ncol + 1); cl++)
		{
			rhs += gsl_vector_get(c, cl) * gsl_matrix_get(X, rw, cl);
		}

		gpred.push_back(rhs);
		tmp = gsl_vector_get(Gm, rw) - rhs;
		gresid.push_back(tmp);
	}
	gsl_vector_free(c);

	// Conduct an initial set of permutations

	Gp = gsl_vector_alloc(nobs);
	npos = 0;
	for (i = 0; i < firstloop; i++)
	{
		// randomly permute residuals

		shuffle(gresid.begin(), gresid.end(), std::default_random_engine(seed));
		seed += 1;

		// compute G* based on marginal L effects and permuted residuals
		for (rw = 0; rw < nobs; rw++)
		{
			gsl_vector_set(Gp, rw, gpred[rw] + gresid[rw]);
		}

		// Recompute p-value for T ~ L|G based on G*
		// fit model T ~ C + G* + L to test L
		stride = 1 + ncolc + 1 + ncol;
		for (rw = 0; rw < nobs; rw++)
		{
			designmat[rw * stride] = 1; // intercept
			for (cl = 0; cl < ncolc; cl++)
			{
				designmat[rw * stride + 1 + cl] = gsl_matrix_get(Cm, rw, cl);
			}
			designmat[rw * stride + 1 + ncolc] = gsl_vector_get(Gp, rw);
			for (cl = 0; cl < ncol; cl++)
			{
				designmat[rw * stride + 2 + ncolc + cl] = gsl_matrix_get(Lm, rw, cl);
			}
		}

		df = ncol;
		converged = linearRegCompare(pvp, phenovec, designmat, nobs, stride, df);
		if (!converged)
			Rcpp::Rcout << "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
		pvp = (converged) ? pvp : std::numeric_limits<double>::quiet_NaN(); // p-value for T ~ L|G*
		if (pvp > pv)
			npos++;

	} // end initial permutation loop

	// Conduct additional permutations if there is some indication of statistical significance
	maxp = maxElementWithNan(pvec);
	nperm = firstloop;
	aa = npos < posno;
	bb = maxp < alpha;
	cc = nperm < maxit;
	testval = (double)(npos + 1) / nperm;
	dd = maxp < testval; // check that other component p-values are small

	if (aa && bb && cc && dd)
	{
		while (aa && cc)
		{

			// randomly permute residuals

			shuffle(gresid.begin(), gresid.end(), std::default_random_engine(seed));
			seed += 1;

			// compute G* based on marginal L effects and permuted residuals
			for (rw = 0; rw < nobs; rw++)
			{
				gsl_vector_set(Gp, rw, gpred[rw] + gresid[rw]);
			}

			// Recompute p-value for T ~ L|G based on G*
			// fit model T ~ C + G* + L to test L
			stride = 1 + ncolc + 1 + ncol; // ip + 1;
			for (rw = 0; rw < nobs; rw++)
			{
				designmat[rw * stride] = 1; // intercept
				for (cl = 0; cl < ncolc; cl++)
				{
					designmat[rw * stride + 1 + cl] = gsl_matrix_get(Cm, rw, cl);
				}
				designmat[rw * stride + 1 + ncolc] = gsl_vector_get(Gp, rw);
				for (cl = 0; cl < ncol; cl++)
				{
					designmat[rw * stride + 2 + ncolc + cl] = gsl_matrix_get(Lm, rw, cl);
				}
			}

			df = ncol;
			converged = linearRegCompare(pvp, phenovec, designmat, nobs, stride, df);
			if (!converged)
				Rcpp::Rcout << "Warning: Cannot Converge when doing regression for calculating P-value." << std::endl;
			pvp = (converged) ? pvp : std::numeric_limits<double>::quiet_NaN(); // p-value for T ~ L|G*
			if (pvp > pv)
				npos++;

			aa = npos < posno;
			cc = nperm < (maxit - 1);
			nperm++;
		} // end 'while' permutation loop
	} // End if
	pv = 1.0 * npos / nperm;

	pvec.push_back(pv); // pval for L ind T|G

	maxp = maxElementWithNan(pvec);
	pval[0] = maxp;
	pval1[0] = pvec[0]; // pval for T ~ L
	pval2[0] = pvec[1]; // pval for T ~ G|L
	pval3[0] = pvec[2]; // pval for G ~ L|T
	pval4[0] = pvec[3]; // pval for L ind T|G

	pvec.clear();
	gresid.clear();
	gpred.clear();
	gsl_matrix_free(Lm);
	gsl_vector_free(Gm);
	gsl_vector_free(Tm);
	if(ncolc > 0){gsl_matrix_free(Cm);}
	if(ncolct > 0){gsl_matrix_free(CGm);}

	gsl_vector_free(Gp);
	gsl_matrix_free(X);

	delete[] designmat;
	delete[] phenovec;

	PutRNGstate();
	LL.clear();
	CC.clear();
	CGG.clear();

} // End citconlog2_linear
