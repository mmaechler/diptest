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

   $Id: dip.c,v 1.19 2003/10/31 21:44:08 maechler Exp $
*/

#include <R.h>

/* Subroutine */
void diptst(double *x, Sint *n,
	    double *dip, Sint *lo_hi, Sint *ifault,
	    Sint *gcm, Sint *lcm, Sint *mn, Sint *mj, Sint *debug)
{
#define low  lo_hi[0]
#define high lo_hi[1]

    /* Local variables */
    int gcmi1, gcmix,  lcmiv, lcmiv1, mnj, mnmnj, mjk, mjmjk,
	ic, icv, icva, icx, icxa, ig, ih, iv, ix,  j, jb, je, jr,  k, kb, ke, kr;
    double dip_l, dip_u, dipnew, d, dx, t, temp, C;
    int N = *n;

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

/* LOW contains the index of the current estimate  of the lower end
   of the modal interval, HIGH contains the index for the upper end.
*/
    low = 1;	high = N; /*-- IDEA:  *xl = x[low];    *xu = x[high]; --*/

    if (N < 4 || x[N] == x[1]) {
	*dip = 0.;	 return;
    }

/* M.Maechler -- speedup: it saves many divisions by N when we just work with
 * (N * dip) everywhere but the very end! */
    *dip = 1.;

    if(*debug) Rprintf("dip() in C: N = %d; starting with N*dip = 1.\n", N);

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
    for (k = N - 1; k >= 1; k--) {
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
	      dx = (lcmiv - gcmi1 + 1) -
		  (x[lcmiv] - x[gcmi1]) * (gcmix - gcmi1)/(x[gcmix] - x[gcmi1]);
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
	      dx = (x[gcmix] - x[lcmiv1]) * (lcmiv - lcmiv1) /
		  (x[lcmiv] - x[lcmiv1])- (gcmix - lcmiv1 - 1);
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
      d = 1.;
    }

    if (d < *dip)	goto L_END;

/*     Calculate the DIPs for the current LOW and HIGH. */

    /* The DIP for the convex minorant. */

    dip_l = 0.;
    if (ig != icx) {
      icxa = icx - 1;
      for (j = ig; j <= icxa; ++j) {
	temp = 1.;
	jb = gcm[j + 1];
	je = gcm[j];
	if (je - jb > 1 && x[je] != x[jb]) {
	  C = (je - jb) / (x[je] - x[jb]);
	  for (jr = jb; jr <= je; ++jr) {
	    t = (jr - jb + 1) - (x[jr] - x[jb]) * C;
	    if (temp < t)
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
	temp = 1.;
	kb = lcm[k];
	ke = lcm[k + 1];
	if (ke - kb > 1 && x[ke] != x[kb]) {
	  C = (ke - kb) / (x[ke] - x[kb]);
	  for (kr = kb; kr <= ke; ++kr) {
	    t = (x[kr] - x[kb]) * C - (kr - kb - 1);
	    if (temp < t)
		temp = t;
	  }
	}
	if (dip_u < temp)
	    dip_u = temp;
      }
    }

    /* Determine the current maximum. */
    dipnew = (dip_u > dip_l) ? dip_u : dip_l;
    if (*dip < dipnew)
	*dip = dipnew;

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
    /* do this in the caller :
     *   *xl = x[low];  *xu = x[high];
     * rather return the (low, high) indices -- automagically via lo_hi[]  */
    *dip /= (2*N);
    return;
} /* diptst */
#undef low
#undef high
