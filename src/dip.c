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
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)

   Pretty--Edited by 	Martin Maechler <maechler@stat.math.ethz.ch>
   			Seminar fuer Statistik, ETH 8092 Zurich	 SWITZERLAND

   $Id: dip.c,v 1.4 1994/07/28 15:18:46 maechler Exp $
*/

/*---- this is OLD  K&R C -- can use 'cc' --> no problem with S-plus 3.2
/* Subroutine */ 
int diptst (x, n, dip, xl, xu, ifault, gcm, lcm, mn, mj)
     float *x; long *n;
     float *dip, *xl, *xu;
     long *ifault, *gcm, *lcm, *mn, *mj;

/*--- This would be ANSI : --
/* int diptst (float *x, long *n, float *dip, 
/* 	    float *xl, float *xu, long *ifault, 
/* 	    long *gcm, long *lcm, long *mn, long *mj)
 */
{
    /* Initialized data */

    static float zero = (float)0.;
    static float half = (float).5;
    static float one = (float)1.;

    /* Local variables */
    static long high, igcm, icva, icxa;
    static float temp;
    static long igcm1;
    static float a, b, d;
    static long j, k;
    static float t;
    static long igcmx, mjmjk, lcmiv, mnmnj;
    static float const__;
    static long lcmiv1, ic, jb, kb, na, ig, ih;
    static float dl;
    static long je;
    static float fn;
    static long jk, ke;
    static float du, dx;
    static long jr, kr, iv, ix;
    static float dipnew;
    static long mjk, icx, mnj, icv, low, lcm1;

    long N = *n;

    /* Parameter adjustments, so I can do "as with index 1 : x[1]..x[N] */
    --mj;
    --mn;
    --lcm;
    --gcm;
    --x;

/*-------- Function Body ------------------------------ */

    *ifault = 1;    if (N <= 0) { return 0; }
    *ifault = 0;

/* Check that X is sorted --- if not, return with  ifault = 2*/

    *ifault = 2;    for (k = 2; k <= N; ++k) if (x[k] < x[k - 1]) return 0;
    *ifault = 0;

/* Check for all values of X identical, */
/*     and for 1 <= N < 4. */

    if (N < 4 || x[N] == x[1]) {	
      *xl = x[1];
      *xu = x[N];
      *dip = zero;      return 0;
    }

/*     LOW contains the index of the current estimate  of the lower end
       of the modal interval, HIGH contains the index for the upper end. 
*/

    low = 1;    high = N;
    fn = (float) (N);
    *dip = one / fn;
    /**-- IDEA:  *xl = x[low];    *xu = x[high]; --*/

/*     Establish the indices over which combination is necessary for the 
       convex MINORANT fit.
*/
    mn[1] = 1;
    for (j = 2; j <= N; ++j) {
	mn[j] = j - 1;
L25:
	mnj = mn[j];
	mnmnj = mn[mnj];
	a = (float) (mnj - mnmnj);
	b = (float) (j - mnj);
	if (mnj == 1 || (x[j] - x[mnj]) * a < (x[mnj] - x[mnmnj]) * b) {
	    goto L28;
	}
	mn[j] = mnmnj;
	goto L25;
L28:
	;
    }

/*     Establish the indices over which combination is necessary for the 
       concave MAJORANT fit. 
*/

    mj[N] = N;
    na = N - 1;
    for (jk = 1; jk <= na; ++jk) {
	k = N - jk;
	mj[k] = k + 1;
L32:
	mjk = mj[k];
	mjmjk = mj[mjk];
	a = (float) (mjk - mjmjk);
	b = (float) (k - mjk);
	if (mjk == N || (x[k] - x[mjk]) * a < (x[mjk] - x[mjmjk]) * b) {
	    goto L34;
	}
	mj[k] = mjmjk;
	goto L32;
L34:
	;
    }

/* ----------------------- Start the cycling. ------------------------------- */
/*     Collect the change points for the GCM from HIGH to LOW. */

LOOP_Start:
    ic = 1;
    gcm[1] = high;

    while(1) {
      igcm1 = gcm[ic];
      ++ic;
      gcm[ic] = mn[igcm1];
      if (gcm[ic] <= low) {	exit;    }
    }
    icx = ic;

/*     Collect the change points for the LCM from LOW to HIGH. */

    ic = 1;
    lcm[1] = low;
L44:
    lcm1 = lcm[ic];
    ++ic;
    lcm[ic] = mj[lcm1];
    if (lcm[ic] < high) {	goto L44;    }
    icv = ic;

/*     ICX, IX, IG are counters for the convex minorant, */
/*     ICV, IV, IH are counters for the concave majorant. */

    ig = icx;
    ih = icv;

/*     Find the largest distance greater than 'DIP' between the GCM and */
/*     the LCM from LOW to HIGH. */

    ix = icx - 1;    iv = 2;    d = zero;
    if (icx != 2 || icv != 2) {	goto L50;    }
    d = one / fn;
    goto L65;
L50:
    igcmx = gcm[ix];
    lcmiv = lcm[iv];
    if (igcmx > lcmiv) {	goto L55;    }

/*     If the next point of either the GCM or LCM is from the LCM, */
/*     calculate the distance here. */

    lcmiv1 = lcm[iv - 1];
    a = (float) (lcmiv - lcmiv1);
    b = (float) (igcmx - lcmiv1 - 1);
    dx = (x[igcmx] - x[lcmiv1] * a) / (fn * (x[lcmiv] - x[lcmiv1])) - b / fn;
    --ix;
    if (dx < d) {	goto L60;    }
    d = dx;
    ig = ix + 1;
    ih = iv;
    goto L60;

/*     If the next point of either the GCM or LCM is from the GCM, */
/*     calculate the distance here. */

L55:
    lcmiv = lcm[iv];
    igcm = gcm[ix];
    igcm1 = gcm[ix + 1];
    a = (float) (lcmiv - igcm1 + 1);
    b = (float) (igcm - igcm1);
    dx = a / fn - (x[lcmiv] - x[igcm1]) * b / (fn * (x[igcm] - x[igcm1]));
    ++iv;
    if (dx < d) {	goto L60;    }
    d = dx;
    ig = ix + 1;
    ih = iv - 1;
L60:
    if (ix < 1) {
	ix = 1;
    }
    if (iv > icv) {
	iv = icv;
    }
    if (gcm[ix] != lcm[iv]) {	goto L50;    }
L65:
    if (d < *dip) {	goto L_END;    }

/*     Calculate the DIPs for the current LOW and HIGH. */

/*     The DIP for the convex minorant. */

    dl = zero;
    if (ig != icx) {
      icxa = icx - 1;
      for (j = ig; j <= icxa; ++j) {
	temp = one / fn;
	jb = gcm[j + 1];
	je = gcm[j];
	if (je - jb > 1 && x[je] != x[jb]) {
	  a = (float) (je - jb);
	  const__ = a / (fn * (x[je] - x[jb]));
	  for (jr = jb; jr <= je; ++jr) {
	    b = (float) (jr - jb + 1);
	    t = b / fn - (x[jr] - x[jb]) * const__;
	    if (t > temp) { temp = t; }
	  }
	}

	if (dl < temp) { dl = temp; }
      }
    }

/*     The DIP for the concave majorant. */

    du = zero;
    if (ih != icv) {
      icva = icv - 1;
      for (k = ih; k <= icva; ++k) {
	temp = one / fn;
	kb = lcm[k];
	ke = lcm[k + 1];
	if (ke - kb > 1 && x[ke] != x[kb]) {
	  a = (float) (ke - kb);
	  const__ = a / (fn * (x[ke] - x[kb]));
	  for (kr = kb; kr <= ke; ++kr) {
	    b = (float) (kr - kb - 1);
	    t = (x[kr] - x[kb]) * const__ - b / fn;
	    if (t > temp) { temp = t; }
	  }
	}
	if (du < temp) { du = temp; }
      }
    }

/*     Determine the current maximum. */

    dipnew = dl;
    if (du > dl) {	dipnew = du;    }
    if (*dip < dipnew) { *dip = dipnew;    }
    low = gcm[ig];
    high = lcm[ih];

    goto LOOP_Start; /* Recycle */
/* ---------------------------------------------------------------------------*/

L_END:
    *dip = half * *dip;
    *xl = x[low];
    *xu = x[high];

    return 0;
} /* diptst */
