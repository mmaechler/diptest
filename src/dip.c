/* dip.f -- translated by f2c (version 19940714).
   You must link the resulting object file with the libraries:
	-lF77 -I77 -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int diptst_(x, n, dip, xl, xu, ifault, gcm, lcm, mn, mj)
real *x;
integer *n;
real *dip, *xl, *xu;
integer *ifault, *gcm, *lcm, *mn, *mj;
{
    /* Initialized data */

    static real zero = (float)0.;
    static real half = (float).5;
    static real one = (float)1.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer high, igcm, icva, icxa;
    static real temp;
    static integer igcm1;
    static real a, b, d;
    static integer j, k;
    static real t;
    static integer igcmx, mjmjk, lcmiv, mnmnj;
    static real const_;
    static integer lcmiv1, ic, jb, kb, na, ig, ih;
    static real dl;
    static integer je;
    static real fn;
    static integer jk, ke;
    static real du, dx;
    static integer jr, kr, iv, ix;
    static real dipnew;
    static integer mjk, icx, mnj, icv, low, lcm1;


/*     ALGORITHM AS 217 APPL. STATIST. (1985) VOL.34, NO.3 */

/*     Does the dip calculation for an ordered vector X using the */
/*     greatest convex minorant and the least concave majorant, skipping 
*/
/*     through the data using the change points of these distributions. */
/*     It returns the dip statistic 'DIP' and the modal interval */
/*     (XL, XU). */

    /* Parameter adjustments */
    --mj;
    --mn;
    --lcm;
    --gcm;
    --x;

    /* Function Body */

    *ifault = 1;
    if (*n <= 0) {
	return 0;
    }
    *ifault = 0;

/*     Check if N = 1 */

    if (*n == 1) {
	goto L4;
    }

/*     Check that X is sorted */

    *ifault = 2;
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	if (x[k] < x[k - 1]) {
	    return 0;
	}
/* L3: */
    }

    *ifault = 0;

/*     Check for all values of X identical, */
/*     and for 1 < N < 4. */

    if (x[*n] > x[1] && *n >= 4) {
	goto L5;
    }
L4:
    *xl = x[1];
    *xu = x[*n];
    *dip = zero;
    return 0;

/*     LOW contains the index of the current estimate of the lower end */
/*     of the modal interval, HIGH contains the index for the upper end. 
*/

L5:
    fn = (real) (*n);
    low = 1;
    high = *n;
    *dip = one / fn;
    *xl = x[low];
    *xu = x[high];

/*     Establish the indices over which combination is necessary for the 
*/
/*     convex minorant fit. */

    mn[1] = 1;
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	mn[j] = j - 1;
L25:
	mnj = mn[j];
	mnmnj = mn[mnj];
	a = (real) (mnj - mnmnj);
	b = (real) (j - mnj);
	if (mnj == 1 || (x[j] - x[mnj]) * a < (x[mnj] - x[mnmnj]) * b) {
	    goto L28;
	}
	mn[j] = mnmnj;
	goto L25;
L28:
	;
    }

/*     Establish the indices over which combination is necessary for the 
*/
/*     concave majorant fit. */

    mj[*n] = *n;
    na = *n - 1;
    i__1 = na;
    for (jk = 1; jk <= i__1; ++jk) {
	k = *n - jk;
	mj[k] = k + 1;
L32:
	mjk = mj[k];
	mjmjk = mj[mjk];
	a = (real) (mjk - mjmjk);
	b = (real) (k - mjk);
	if (mjk == *n || (x[k] - x[mjk]) * a < (x[mjk] - x[mjmjk]) * b) {
	    goto L34;
	}
	mj[k] = mjmjk;
	goto L32;
L34:
	;
    }

/*     Start the cycling. */
/*     Collect the change points for the GCM from HIGH to LOW. */

L40:
    ic = 1;
    gcm[1] = high;
L42:
    igcm1 = gcm[ic];
    ++ic;
    gcm[ic] = mn[igcm1];
    if (gcm[ic] > low) {
	goto L42;
    }
    icx = ic;

/*     Collect the change points for the LCM from LOW to HIGH. */

    ic = 1;
    lcm[1] = low;
L44:
    lcm1 = lcm[ic];
    ++ic;
    lcm[ic] = mj[lcm1];
    if (lcm[ic] < high) {
	goto L44;
    }
    icv = ic;

/*     ICX, IX, IG are counters for the convex minorant, */
/*     ICV, IV, IH are counters for the concave majorant. */

    ig = icx;
    ih = icv;

/*     Find the largest distance greater than 'DIP' between the GCM and */
/*     the LCM from LOW to HIGH. */

    ix = icx - 1;
    iv = 2;
    d = zero;
    if (icx != 2 || icv != 2) {
	goto L50;
    }
    d = one / fn;
    goto L65;
L50:
    igcmx = gcm[ix];
    lcmiv = lcm[iv];
    if (igcmx > lcmiv) {
	goto L55;
    }

/*     If the next point of either the GCM or LCM is from the LCM, */
/*     calculate the distance here. */

    lcmiv1 = lcm[iv - 1];
    a = (real) (lcmiv - lcmiv1);
    b = (real) (igcmx - lcmiv1 - 1);
    dx = (x[igcmx] - x[lcmiv1] * a) / (fn * (x[lcmiv] - x[lcmiv1])) - b / fn;
    --ix;
    if (dx < d) {
	goto L60;
    }
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
    a = (real) (lcmiv - igcm1 + 1);
    b = (real) (igcm - igcm1);
    dx = a / fn - (x[lcmiv] - x[igcm1]) * b / (fn * (x[igcm] - x[igcm1]));
    ++iv;
    if (dx < d) {
	goto L60;
    }
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
    if (gcm[ix] != lcm[iv]) {
	goto L50;
    }
L65:
    if (d < *dip) {
	goto L100;
    }

/*     Calculate the DIPs for the current LOW and HIGH. */

/*     The DIP for the convex minorant. */

    dl = zero;
    if (ig == icx) {
	goto L80;
    }
    icxa = icx - 1;
    i__1 = icxa;
    for (j = ig; j <= i__1; ++j) {
	temp = one / fn;
	jb = gcm[j + 1];
	je = gcm[j];
	if (je - jb <= 1) {
	    goto L74;
	}
	if (x[je] == x[jb]) {
	    goto L74;
	}
	a = (real) (je - jb);
	const_ = a / (fn * (x[je] - x[jb]));
	i__2 = je;
	for (jr = jb; jr <= i__2; ++jr) {
	    b = (real) (jr - jb + 1);
	    t = b / fn - (x[jr] - x[jb]) * const_;
	    if (t > temp) {
		temp = t;
	    }
/* L72: */
	}
L74:
	if (dl < temp) {
	    dl = temp;
	}
/* L76: */
    }

/*     The DIP for the concave majorant. */

L80:
    du = zero;
    if (ih == icv) {
	goto L90;
    }
    icva = icv - 1;
    i__1 = icva;
    for (k = ih; k <= i__1; ++k) {
	temp = one / fn;
	kb = lcm[k];
	ke = lcm[k + 1];
	if (ke - kb <= 1) {
	    goto L86;
	}
	if (x[ke] == x[kb]) {
	    goto L86;
	}
	a = (real) (ke - kb);
	const_ = a / (fn * (x[ke] - x[kb]));
	i__2 = ke;
	for (kr = kb; kr <= i__2; ++kr) {
	    b = (real) (kr - kb - 1);
	    t = (x[kr] - x[kb]) * const_ - b / fn;
	    if (t > temp) {
		temp = t;
	    }
/* L84: */
	}
L86:
	if (du < temp) {
	    du = temp;
	}
/* L88: */
    }

/*     Determine the current maximum. */

L90:
    dipnew = dl;
    if (du > dl) {
	dipnew = du;
    }
    if (*dip < dipnew) {
	*dip = dipnew;
    }
    low = gcm[ig];
    high = lcm[ih];

/*     Recycle */

    goto L40;

L100:
    *dip = half * *dip;
    *xl = x[low];
    *xu = x[high];

    return 0;
} /* diptst_ */

