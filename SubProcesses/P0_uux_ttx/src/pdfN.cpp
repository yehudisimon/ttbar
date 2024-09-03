// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Implements the PDF module and the PDF fit (needed for Mellin space PDFs).

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <complex>
#include "process.h"

#include "gsl/gsl_blas.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_monte_vegas.h"
#include "gsl/gsl_multifit_nlinear.h"
#include "gsl/gsl_sf.h"
#include "gsl/gsl_types.h"

using namespace std;

#define NFLVR 5   // Number of active flavours.
#define NDATA 10000 // Number of data points for the fit

void SetLHAPDF(const PDF *, double&, double*, double*, double&);

extern double muF, xmin, xmax;
extern string namePDF;

// function used to sample the PDFs.
double sampling(double xmin2, double xmax2 ,const size_t n, int i) {
    return pow(xmin2, 1.0 - (1.0 + (double)i) / (xmax2 + (double)n));
}


// Struct for GSL fit.
struct data {
  size_t n;
  double *y;
  double *t;
  double *sigma;
  double xr;
};

int expb_f(const gsl_vector *x, void *data, gsl_vector *f) {
    size_t n = ((struct data *)data)->n;
    double *y = ((struct data *)data)->y;
    double *sigma = ((struct data *)data)->sigma;
    double xr = ((struct data *)data)->xr;

    double A[8];
    for (size_t i0 = 0; i0 < 8; i0++) {
        A[i0] = gsl_vector_get(x, i0);
    }

    for (size_t i0 = 0; i0 < n; i0++) {
      double t = sampling(xr,xmax,n,i0);
      double Yi = A[0] * pow(t, A[1]) * pow(1.0 - t, A[2])* (1.0 + A[3] * sqrt(t) + A[4] * t + A[5] * pow(t, 1.5)+ A[6] * pow(t, 2.0) + A[7] * pow(t, 2.5));
      gsl_vector_set(f, i0, (Yi - y[i0]) /sigma[i0]);
    }
	
	return GSL_SUCCESS;
}

int expb_df(const gsl_vector *x, void *data, gsl_matrix *J) {
    size_t n     = ((struct data *)data)->n;
    double xr    = ((struct data *)data)->xr;
    double *sigma = ((struct data *)data)->sigma;
    double A[8];
    for (size_t i0 = 0; i0 < 8; i0++) {
        A[i0] = gsl_vector_get(x, i0);
    }

    for (size_t i0 = 0; i0 < n; i0++) {
        double t = sampling(xr,xmax,n,i0);
	double ssigma=sigma[i0];
	double e = A[0] * pow(t, A[1]) * pow(1.0 - t, A[2])* (1.0 + A[3] * sqrt(t) + A[4] * t + A[5] * pow(t, 1.5)+ A[6] * pow(t, 2.0) + A[7] * pow(t, 2.5))/ssigma;
        gsl_matrix_set(J, i0, 0, e / A[0]);
        gsl_matrix_set(J, i0, 1, e * log(t));
        gsl_matrix_set(J, i0, 2, e * log(1.0 - t));

    	e = A[0] * pow(t, A[1]) * pow(1.0 - t, A[2])/ssigma;
        gsl_matrix_set(J, i0, 3, e * sqrt(t));
        gsl_matrix_set(J, i0, 4, e * t);
        gsl_matrix_set(J, i0, 5, e * pow(t, 1.5));
        gsl_matrix_set(J, i0, 6, e * pow(t, 2.0));
        gsl_matrix_set(J, i0, 7, e * pow(t, 2.5));
    }

    return GSL_SUCCESS;
}

void
callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double rcond;

  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);
}

void Fit(double A[8], double E[8], int flag, double xr, double weight_valence, double weight_sea, double weight_gluon) {
    // Defines the function to minimize.
    const size_t n = NDATA;
    const size_t p = 8;
    
    double y[NDATA], sigma[NDATA], t[NDATA];


    PDF* F=mkPDF(namePDF,0);
    struct data d = { n, y, t, sigma, xr };
    // First guess for the parameters.
    double x_init[8] = {1.0, -1.4, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    if (flag == 1 || flag == 2) {
        x_init[1] = -0.6;
    }

    gsl_vector_view x = gsl_vector_view_array(x_init, p);
    gsl_vector_view wts = gsl_vector_view_array(sigma, n);
    gsl_multifit_nlinear_fdf f;
    f.f = &expb_f;
    f.df = &expb_df;
    f.fvv=NULL;
    //f.fdf = &expb_fdf; //Old form
    f.n = n;
    f.p = p;
    f.params = &d;

    // PDFs to fit.
    double gamma_weight = -1.6;
    for (size_t i0 = 0; i0 < n; i0++)
      {
	double tc = sampling(xr,xmax,n,i0);
	double q[2][5];
	double g;
	double qq;
	
	SetLHAPDF(F,tc,q[0],q[1],g);
	switch (flag) {
	case 0:
	  qq = g;
	  gamma_weight = weight_gluon;
	  break;
	case 1:
	  qq = q[0][0]; //d
	  gamma_weight = weight_valence;
	  break;
	case 2:
	  qq = q[0][1]; //u
	  gamma_weight = weight_valence;
	  break;
	case 3:
	  qq = q[1][0]; //dbar
	  gamma_weight = weight_sea;
	  break;
	case 4:
	  qq = q[1][2]; //sbar or s
	  gamma_weight = weight_sea;
	  break;
	case 5:
	  qq = q[1][4]; //bbar or b
	  gamma_weight = weight_sea;
	  break;
	case 6:
	  qq = q[1][1]; //ubar
	  gamma_weight = weight_sea;
	  break;
	case 7:
	  qq = q[1][3]; //cbar or c
	  gamma_weight = weight_sea;
	  break;
	default:
	  fprintf(stderr,
		  "error: while retrieving PDF: unkown flag %d\n", flag);
	  exit(1);
	}

        y[i0] = qq;
	t[i0]=tc;
        sigma[i0]= pow(tc,gamma_weight); // DeFlorian-like pdf-weights t^(-1.6); DeJonathan t^(-1)
    }


    delete F;
    
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();


    w=gsl_multifit_nlinear_alloc(T,&fdf_params,n,p);
    gsl_multifit_nlinear_winit(&x.vector, &wts.vector,&f,w);

    /* compute initial cost function */
    gsl_vector *f2;
    double chisq0, chisq;
    f2 = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f2, f2, &chisq0);

    /* solve the system with a maximum of 100000 iterations */
    int info;
    const double xtol = 1e-18, gtol=1e-18, ftol=0.0;
    int status;
    status= gsl_multifit_nlinear_driver(1000, xtol, gtol, ftol, NULL, NULL, &info, w);

    gsl_matrix *Jac;
    Jac=gsl_multifit_nlinear_jac(w);
    
    
    gsl_matrix *covar;
    covar= gsl_matrix_alloc(p, p);
    gsl_multifit_nlinear_covar(Jac, 0.0, covar);
    gsl_blas_ddot(f2, f2, &chisq);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, sqrt(chisq) / sqrt(dof));

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
    
    // Gets parameters.
    for (size_t i0 = 0; i0 < 8; i0++) {
        A[i0] = 0.0;
    }
    for (size_t i0 = 0; i0 < p; i0++) {
        A[i0] = gsl_vector_get(w->x, i0);
    }
    for (size_t i0 = 0; i0 < p; i0++) {
        E[i0] = c * gsl_matrix_get(covar, i0, i0);
    }

    //    gsl_multifit_fdfsolver_free(s); //Old
    gsl_multifit_nlinear_free(w);
    gsl_matrix_free(covar);
}

void pdfFit(double &A1MIN, double A[8][8], double weight_valence, double weight_sea, double weight_gluon) {
    const int nf = NFLVR;
    //fprintf(stderr,"Performing PDF fit with %d flavors, Q^2 = mu_F^2 = %g\n and weights: valence: x^%g, sea: x^%g, gluon: x^%g and xmin = %g \n Fit function: f = A0 * x^A1 * (1 - x)^A2 * ( 1 + A3 * x^(1/2) + A4 * x + A5 * x^(3/2) + A6 * x^2 + A7 * x^(5/2) )\n", nf, muF*muF, weight_valence, weight_sea, weight_gluon, xmin);

    for (int i0 = 0; i0 < 8; i0++) {
        double err[8];
        switch (i0) {
        case 0:
	  //fprintf(stderr, "Fitting gluon PDF...");
	  break;
        case 1:
	  //fprintf(stderr, "Fitting valence down quark PDF...");
	  break;
        case 2:
	  //fprintf(stderr, "Fitting valence up quark PDF...");
	  break;
        case 3:
	  //fprintf(stderr, "Fitting sea down quark PDF...");
	  break;
        case 4:
	  //fprintf(stderr, "Fitting strange quark PDF...");
	  break;
        case 5:
	  //fprintf(stderr, "Fitting bottom quark PDF...");
	  break;
        case 6:
	  //fprintf(stderr, "Fitting sea up quark PDF...");
	  break;
        case 7:
	  //fprintf(stderr, "Fitting charm quark PDF...");
	  break;
        }
	
        Fit(A[i0], err, i0, xmin, weight_valence,  weight_sea,  weight_gluon); 
        //for (int i1 = 0; i1 < 8; i1++) {
          //fprintf(stderr, "A%i =  %.5f # +-%.5f\n",i1, A[i0][i1], err[i1]);
        //}
        if (A[i0][1] < A1MIN) {
            A1MIN = A[i0][1];
        }
    }
}
