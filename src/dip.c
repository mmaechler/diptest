/*   ALGORITHM AS 217 APPL. STATIST. (1985) VOL.34, NO.3 

  @article{HarP85,
     author = {P. M. Hartigan},
     title = {Computation of the Dip Statistic to Test for Unimodality},
     year = 1985,
     journal = {Applied Statistics},
     pages = {320--325},
     volume = 34 }
  @article{HarJH85,
     author = {J. A. Hartigan and P. M. Hartigan},
     title = {The Dip Test of Unimodality},
     year = 1985,
     journal = {Ann. of Statistics},
     pages = {70--84},
     volume = 13 }

  Does the dip calculation for an ordered vector X using the 
  greatest convex minorant and the least concave majorant, skipping 
  through the data using the change points of these distributions. 

  It returns the dip statistic 'DIP' and the modal interval (XL, XU).
				===                          ======

   dip.f -- translated by f2c (version of 22 July 1992  22:54:52).

   Pretty--Edited by 	Martin Maechler <maechler@stat.math.ethz.ch>
   			Seminar fuer Statistik, ETH 8092 Zurich	 SWITZERLAND

   $Id: dip.c,v 1.10 2000/12/12 21:56:06 mm Exp $
*/
#include <stdio.h> /*--- for debugging only ---*/

/* Subroutine */ 
int diptst (float *x, long int *n, float *dip, float *xl, float *xu, long int *ifault, long int *gcm, long int *lcm, long int *mn, long int *mj, long int *debug)
{
    /* Initialized data */

    static float zero = (float)0.;
    static float one  = (float)1.;

    /* Local variables */
    long low, high,  gcmi, gcmi1, gcmix,  lcm1, lcmiv, lcmiv1, 
	mnj, mnmnj, mjk, mjmjk,   ic, icv, icva, icx, icxa,
	ig, ih, iv, ix, j, jb, je, jk, jr, k, kb, ke, kr, N1;
    float fn, dip_l, dip_u, dipnew, a, b, d, dx, t, temp, C;

    long N = *n;

    /* Parameter adjustments, so I can do "as with index 1" : x[1]..x[N] */
    --mj;    --mn;
    --lcm;   --gcm;
    --x;
    N1 = N - 1;

/*-------- Function Body ------------------------------ */

    *ifault = 1;    if (N <= 0) { return 0; }
    *ifault = 0;

/* Check that X is sorted --- if not, return with  ifault = 2*/

    *ifault = 2;    for (k = 2; k <= N; ++k) if (x[k] < x[k - 1]) return 0;
    *ifault = 0;

/* Check for all values of X identical, */
/*     and for 1 <= N < 4. */

    if (N < 4 || x[N] == x[1]) {	
      *xl = x[1];  *xu = x[N];   *dip = zero;      return 0;
    }

/* LOW contains the index of the current estimate  of the lower end
   of the modal interval, HIGH contains the index for the upper end. 
*/

    low = 1;    high = N; /*-- IDEA:  *xl = x[low];    *xu = x[high]; --*/
    fn = (float) (N);
    *dip = one / fn;

    if(*debug) printf( "'dip': starting with dip = %5g\n", *dip);

/* Establish the indices   mn[1..N]  over which combination is necessary
   for the convex MINORANT (GCM) fit.
*/
    mn[1] = 1;
    for (j = 2; j <= N; ++j) {
	mn[j] = j - 1;
	while(1) {
	  mnj = mn[j];
	  mnmnj = mn[mnj];
	  a = (float) (mnj - mnmnj);
	  b = (float) (j - mnj);
	  if (mnj == 1 ||
	      (x[j] - x[mnj]) * a < (x[mnj] - x[mnmnj]) * b) break;
	  mn[j] = mnmnj;
	}
    }

/* Establish the indices   mj[1..N]  over which combination is necessary
   for the concave MAJORANT (LCM) fit. 
*/
    mj[N] = N;
    for (jk = 1; jk <= N1; ++jk) {
	k = N - jk;
	mj[k] = k + 1;
	while(1) {
	  mjk = mj[k];
	  mjmjk = mj[mjk];
	  a = (float) (mjk - mjmjk);
	  b = (float) (k - mjk);
	  if (mjk == N ||
	      (x[k] - x[mjk]) * a < (x[mjk] - x[mjmjk]) * b) break;
	  mj[k] = mjmjk;
	}
    }

/* ----------------------- Start the cycling. ------------------------------- */
LOOP_Start:

/* Collect the change points for the GCM from HIGH to LOW. */

    if(*debug) printf( "'dip':LOOP-BEGIN: low, high = %5ld,%5ld\n", low,high);

    ic = 1; gcm[1] = high;
    while(gcm[ic] > low) {
      gcmi1 = gcm[ic];
      ++ic;
      gcm[ic] = mn[gcmi1];
    }
    icx = ic;

/* Collect the change points for the LCM from LOW to HIGH. */

    ic = 1; lcm[1] = low;
    while(lcm[ic] < high) {
      lcm1 = lcm[ic];
      ++ic;
      lcm[ic] = mj[lcm1];
    }
    icv = ic;

/*     ICX, IX, IG are counters for the convex  minorant, */
/*     ICV, IV, IH are counters for the concave majorant. */

    ig = icx;
    ih = icv;

/*     Find the largest distance greater than 'DIP' between the GCM and */
/*     the LCM from LOW to HIGH. */

    ix = icx - 1;    iv = 2;    d = zero;
    if (icx != 2 || icv != 2) {
      while(1) { /* gcm[ix] != lcm[iv]  (after first loop) */
	gcmix = gcm[ix];
	lcmiv = lcm[iv];
	if (gcmix > lcmiv) {	
	  lcmiv = lcm[iv];
	  gcmi  = gcm[ix];
	  gcmi1 = gcm[ix + 1];
	  a = (float) (lcmiv - gcmi1 + 1);
	  b = (float) (gcmi - gcmi1);
	  dx = a / fn - (x[lcmiv] - x[gcmi1]) * b / (fn * (x[gcmi] - x[gcmi1]));
	  ++iv;
	  if (dx >= d) {
	    d = dx;
	    ig = ix + 1;
	    ih = iv - 1;
	  }
	} else {
	  /*     If the next point of either the GCM or LCM is from the LCM, */
	  /*     calculate the distance here. */

	  lcmiv1 = lcm[iv - 1];
	  a = (float) (lcmiv - lcmiv1);
	  b = (float) (gcmix - lcmiv1 - 1);
	  dx = (x[gcmix] - a* x[lcmiv1]) / (fn*(x[lcmiv] - x[lcmiv1])) - b/fn;
	  --ix;
	  if (dx >= d) {
	    d = dx;
	    ig = ix + 1;
	    ih = iv;
	  }
	}

	/*     If the next point of either the GCM or LCM is from the GCM, */
	/*     calculate the distance here. */

	if (ix < 1)		ix = 1;
	if (iv > icv)	iv = icv;
	if (gcm[ix] == lcm[iv]) break;
      }
    } else { /* icx or icv == 2 */
      d = one / fn;
    }

    if (d < *dip)	goto L_END;

/*     Calculate the DIPs for the current LOW and HIGH. */

    /* The DIP for the convex minorant. */

    dip_l = zero;
    if (ig != icx) {
      icxa = icx - 1;
      for (j = ig; j <= icxa; ++j) {
	temp = one / fn;
	jb = gcm[j + 1];
	je = gcm[j];
	if (je - jb > 1 && x[je] != x[jb]) {
	  a = (float) (je - jb);
	  C = a / (fn * (x[je] - x[jb]));
	  for (jr = jb; jr <= je; ++jr) {
	    b = (float) (jr - jb + 1);
	    t = b / fn - (x[jr] - x[jb]) * C;
	    if (t > temp) { temp = t; }
	  }
	}
	if (dip_l < temp) dip_l = temp;
      }
    }

    /* The DIP for the concave majorant. */

    dip_u = zero;
    if (ih != icv) {
      icva = icv - 1;
      for (k = ih; k <= icva; ++k) {
	temp = one / fn;
	kb = lcm[k];
	ke = lcm[k + 1];
	if (ke - kb > 1 && x[ke] != x[kb]) {
	  a = (float) (ke - kb);
	  C = a / (fn * (x[ke] - x[kb]));
	  for (kr = kb; kr <= ke; ++kr) {
	    b = (float) (kr - kb - 1);
	    t = (x[kr] - x[kb]) * C - b / fn;
	    if (t > temp) temp = t;
	  }
	}
	if (dip_u < temp) dip_u = temp;
      }
    }

    /* Determine the current maximum. */

    dipnew = dip_l;
    if (dip_u > dip_l) dipnew = dip_u;
    if (*dip < dipnew)   *dip = dipnew;
    /*--- The following 'if' is NECESSARY ! ------------------------------
      --- Martin Maechler, Statistics, ETH Zurich, July 30 1994 ---------- */
    if (low == gcm[ig] && high == lcm[ih]) {
      if(*debug) 
	printf("No improvement in  low = %ld  nor  high = %ld --> END\n",
	       low, high);
    } else {
      low  = gcm[ig];
      high = lcm[ih];      goto LOOP_Start; /* Recycle */
    }
/*---------------------------------------------------------------------------*/

L_END:
    *xl = x[low];  *xu = x[high];  *dip = (float)0.5 * *dip;    return 0;
} /* diptst */
