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
				===			     ======

   dip.f -- translated by f2c (version of 22 July 1992	22:54:52).

   Pretty-Edited and extended (debug argument)
   by Martin Maechler <maechler@stat.math.ethz.ch>
      ETH Seminar fuer Statistik
      8092 Zurich	 SWITZERLAND

   $Id: dip.c,v 1.16 2003/10/31 18:17:32 maechler Exp $
*/

#include <R.h>

/* Subroutine */
void diptst(double *x, Sint *n,
	    double *dip, double *xl, double *xu,
	    Sint *ifault,
	    Sint *gcm, Sint *lcm,
	    Sint *mn, Sint *mj, Sint *debug)
{
    /* Local variables */
    int low, high, gcmi1, gcmix,  lcmiv, lcmiv1,
	mnj, mnmnj, mjk, mjmjk,	  ic, icv, icva, icx, icxa,
	ig, ih, iv, ix, j, jb, je, jk, jr, k, kb, ke, kr;
    double dip_l, dip_u, dipnew, d, dx, t, temp, C;
    int N = *n, N1 = N - 1;
    double fN = (double)N;

    /* Parameter adjustments, so I can do "as with index 1" : x[1]..x[N] */
    --mj;    --mn;
    --lcm;   --gcm;
    --x;

/*-------- Function Body ------------------------------ */

    *ifault = 1;    if (N <= 0) return;
    *ifault = 0;

/* Check that X is sorted --- if not, return with  ifault = 2*/

    *ifault = 2;    for (k = 2; k <= N; ++k) if (x[k] < x[k - 1]) return;
    *ifault = 0;

/* Check for all values of X identical, */
/*     and for 1 <= N < 4. */

    if (N < 4 || x[N] == x[1]) {
      *xl = x[1];  *xu = x[N];	 *dip = 0.;	 return;
    }

/* LOW contains the index of the current estimate  of the lower end
   of the modal interval, HIGH contains the index for the upper end.
*/

    low = 1;	high = N; /*-- IDEA:  *xl = x[low];    *xu = x[high]; --*/
    *dip = 1. / N;

    if(*debug) Rprintf("'dip': starting with dip = %5g\n", *dip);

/* Establish the indices   mn[1..N]  over which combination is necessary
   for the convex MINORANT (GCM) fit.
*/
    mn[1] = 1;
    for (j = 2; j <= N; ++j) {
	mn[j] = j - 1;
	while(1) {
	  mnj = mn[j];
	  mnmnj = mn[mnj];
	  if (mnj == 1 ||
	      ( x[j]  - x[mnj]) * (mnj - mnmnj) <
	      (x[mnj] - x[mnmnj]) * (j - mnj)) break;
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
	  if (mjk == N ||
	      ( x[k]  - x[mjk]) * (mjk - mjmjk) <
	      (x[mjk] - x[mjmjk]) * (k - mjk)) break;
	  mj[k] = mjmjk;
	}
    }

/* ----------------------- Start the cycling. ------------------------------- */
LOOP_Start:

    if(*debug) Rprintf( "'dip':LOOP-BEGIN: low, high = %5ld,%5ld\n", low,high);

/* Collect the change points for the GCM from HIGH to LOW. */

    gcm[1] = high;
    for(ic = 1; gcm[ic] > low; ic++)
	gcm[ic+1] = mn[gcm[ic]];
    icx = ic;

/* Collect the change points for the LCM from LOW to HIGH. */

    lcm[1] = low;
    for(ic = 1; lcm[ic] < high; ic++)
	lcm[ic+1] = mj[lcm[ic]];
    icv = ic;

    ig = icx; /* ICX, IX, IG are counters for the convex minorant, */
    ih = icv; /* ICV, IV, IH are counters for the concave majorant. */
    ix = icx - 1;    iv = 2;

/*	Find the largest distance greater than 'DIP' between the GCM and
 *	the LCM from LOW to HIGH. */

    d = 0.;
    if (icx != 2 || icv != 2) {
      do { /* gcm[ix] != lcm[iv]  (after first loop) */
	  gcmix = gcm[ix];
	  lcmiv = lcm[iv];
	  if (gcmix > lcmiv) {
	      /* If the next point of either the GCM or LCM is from the LCM,
	       * calculate the distance here. */
	      gcmi1 = gcm[ix + 1];
	      dx = (lcmiv - gcmi1 + 1) / fN -
		  (x[lcmiv] - x[gcmi1]) * (gcmix-gcmi1)/(N*(x[gcmix] - x[gcmi1]));
	      ++iv;
	      if (dx >= d) {
		  d = dx;
		  ig = ix + 1;
		  ih = iv - 1;
	      }
	  }
	  else {
	      /* If the next point of either the GCM or LCM is from the GCM,
	       * calculate the distance here. */
	      lcmiv1 = lcm[iv - 1];
/* original
	      dx = (x[gcmix] - (lcmiv - lcmiv1)* x[lcmiv1]) /
		  (N*(x[lcmiv] - x[lcmiv1])) - (gcmix - lcmiv1 - 1) / fN;
*/
/* Fix by Yong Lu  {is more symmetric to the above case}:*/
	      dx = (x[gcmix] - x[lcmiv1]) * (lcmiv-lcmiv1) /
		  (N*(x[lcmiv] - x[lcmiv1]))- (gcmix - lcmiv1 - 1) / fN;
	      --ix;
	      if (dx >= d) {
		  d = dx;
		  ig = ix + 1;
		  ih = iv;
	      }
	  }

	  if (ix < 1)	ix = 1;
	  if (iv > icv)	iv = icv;
      } while (gcm[ix] != lcm[iv]);
    }
    else { /* icx or icv == 2 */
      d = 1. / fN;
    }

    if (d < *dip)	goto L_END;

/*     Calculate the DIPs for the current LOW and HIGH. */

    /* The DIP for the convex minorant. */

    dip_l = 0.;
    if (ig != icx) {
      icxa = icx - 1;
      for (j = ig; j <= icxa; ++j) {
	temp = 1. / fN;
	jb = gcm[j + 1];
	je = gcm[j];
	if (je - jb > 1 && x[je] != x[jb]) {
	  C = (je - jb) / (fN * (x[je] - x[jb]));
	  for (jr = jb; jr <= je; ++jr) {
	    t = (jr - jb + 1) / fN - (x[jr] - x[jb]) * C;
	    if (t > temp)
		temp = t;
	  }
	}
	if (dip_l < temp)
	  dip_l = temp;
      }
    }

    /* The DIP for the concave majorant. */

    dip_u = 0.;
    if (ih != icv) {
      icva = icv - 1;
      for (k = ih; k <= icva; ++k) {
	temp = 1. / fN;
	kb = lcm[k];
	ke = lcm[k + 1];
	if (ke - kb > 1 && x[ke] != x[kb]) {
	  C = (ke - kb) / (fN * (x[ke] - x[kb]));
	  for (kr = kb; kr <= ke; ++kr) {
	    t = (x[kr] - x[kb]) * C - (kr - kb - 1) / fN;
	    if (t > temp) temp = t;
	  }
	}
	if (dip_u < temp) dip_u = temp;
      }
    }

    /* Determine the current maximum. */

    dipnew = dip_l;
    if (dip_u > dip_l) dipnew = dip_u;
    if (*dip < dipnew)	 *dip = dipnew;
    /*--- The following 'if' is NECESSARY ! ------------------------------
      --- Martin Maechler, Statistics, ETH Zurich, July 30 1994 ---------- */
    if (low == gcm[ig] && high == lcm[ih]) {
      if(*debug)
	Rprintf("No improvement in  low = %ld  nor  high = %ld --> END\n",
	       low, high);
    } else {
      low  = gcm[ig];
      high = lcm[ih];	   goto LOOP_Start; /* Recycle */
    }
/*---------------------------------------------------------------------------*/

L_END:
    *xl = x[low];  *xu = x[high];
    *dip /= 2;
    return;
} /* diptst */
