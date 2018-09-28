/* plasma.f -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"
#include <stdio.h>
#include "run_control.h"

/* Common Block Declarations */

struct cag05b_1_ {
    doublereal store1, store2;
};
struct cag05b_2_ {
    doublereal normal, gamma;
};

#define cag05b_1 (*(struct cag05b_1_ *) &cag05b_)
#define cag05b_2 (*(struct cag05b_2_ *) &cag05b_)

struct cag05a_1_ {
    integer ix, iy, iz;
};

#define cag05a_1 (*(struct cag05a_1_ *) &cag05a_)

/* Initialized data */

struct {
    doublereal e_1[2];
    } cag05b_ = { 1., -1. };

struct {
    integer e_1[3];
    } cag05a_ = { 1, 255, 25555 };


/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static integer c__0 = 0;
static integer c__2 = 2;
static doublereal c_b172 = 1.;
static integer c__12 = 12;
static doublereal c_b417 = 1.1424;
static doublereal c_b419 = .5712;
static doublereal c_b429 = 1.5;
static integer c__20000 = 20000;
static doublereal c_b581 = 0.;
static doublereal c_b582 = .11200000000000002;
static doublereal c_b587 = .14;
static doublereal c_b589 = .8;
static doublereal c_b614 = .001;

/* ------------------------------------------------------ */
/* ������� �.�., */
doublereal g05dde_(doublereal *a, doublereal *b)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal half = .5;
    static doublereal d__[41] = { 0.,.674489750196082,1.150349380376008,
	    1.534120544352546,1.862731867421652,2.153874694061456,
	    2.417559016236505,2.66006746861746,2.885634912426757,
	    3.097269078198785,3.297193345691964,3.487104104114431,
	    3.668329285121323,3.841930685501911,4.008772594168585,
	    4.169569323349106,4.324919040826046,4.475328424654204,
	    4.621231001499247,4.763001034267814,4.900964207963193,
	    5.035405969463927,5.166578119728753,5.294704084854598,
	    5.419983174916868,5.54259405780294,5.662697617459439,
	    5.780439324478935,5.895951216739571,6.009353565530745,
	    6.120756285971941,6.230260137989044,6.33795775455379,
	    6.443934526538564,6.548269367831731,6.651035379893011,
	    6.752300431407015,6.852127665896068,6.95057594791675,
	    7.047700256664409,7.14355203435219 };

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer n;
    static doublereal t, u, v, w, x;
    extern doublereal wrapg05cae_(doublereal *);

    u = cag05b_1.store1;
    for (n = 1; n <= 39; ++n) {
	if (u > half) {
	    goto L40;
	}
	u += u;
/* L20: */
    }
    n = 40;
L40:
    t = d__[n - 1];
    u = wrapg05cae_(&x);
L60:
    w = (d__[n] - t) * u;
    v = w * (w * half + t);
L80:
    u = wrapg05cae_(&x);
    if (v <= u) {
	goto L100;
    }
    v = wrapg05cae_(&x);
    if (u > v) {
	goto L80;
    }
    u = (v - u) / (one - u);
    goto L60;
L100:
    u = (u - v) / (one - v);
    if (u > half) {
	goto L120;
    }
    cag05b_1.store1 = u + u;
    ret_val = *a + *b * (w + t);
    return ret_val;
L120:
    cag05b_1.store1 = u + u - one;
    ret_val = *a - *b * (w + t);
    return ret_val;
} /* g05dde_ */

/* ------------------------------------------------------------------ */
doublereal g05cae_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal ai;
    static integer ii;
    static doublereal ax, ay, az;

    cag05a_1.ix = (cag05a_1.ix - cag05a_1.ix / 177 * 177) * 171 - (
	    cag05a_1.ix / 177 << 1);
    cag05a_1.iy = (cag05a_1.iy - cag05a_1.iy / 176 * 176) * 172 - (
	    cag05a_1.iy / 176 << 1);
    cag05a_1.iz = (cag05a_1.iz - cag05a_1.iz / 178 * 178) * 170 - (
	    cag05a_1.iz / 178 << 1);
    if (cag05a_1.ix < 0) {
	cag05a_1.ix += 30269;
    }
    if (cag05a_1.iy < 0) {
	cag05a_1.iy += 30307;
    }
    if (cag05a_1.iz < 0) {
	cag05a_1.iz += 30323;
    }
    ax = (doublereal) cag05a_1.ix;
    ay = (doublereal) cag05a_1.iy;
    az = (doublereal) cag05a_1.iz;
    ai = ax / 30269. + ay / 30307. + az / 30323.;
    ii = (integer) ai;
    ret_val = ai - ii;
    return ret_val;
} /* g05cae_ */

doublereal wrapg05cae_(doublereal *x)
{
    double t  = g05cae_(x);
#ifdef DEBUG_PLASMA
    printf("%25.15e \n",t);
#endif

    return t;
}

doublereal wrapg05dde_(doublereal *a, doublereal *b)
{
    return g05dde_(a,b);
}

double rnd_uniform()
{
    doublereal x;

    return (double)wrapg05cae_(&x);
}

double rnd_gaussian(double a,double b)
{
    return (double)wrapg05dde_(&a,&b);
}

