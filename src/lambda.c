#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "dist.h"
#include "poly.h"

// *******************************************************************************
// integral of f with respect to r from 0 to r
double fr(double r, double w[])
{
  double gamma = w[0]; //changed
  return (1 - exp(- r * r / (2 * gamma * gamma))) / (2 * PI);//changed
}
// derivative of fr with respect to gamma
double dgamma_fr(double r, double w[])
{
  double gamma = w[0]; //changed
  return - exp(- r * r / (2 * gamma * gamma)) * r * r / (2 * PI * gamma * gamma * gamma);//changed
}
// derivative of fr with respect to C
/* double dC_fr(double r, double w[])
{
  double gamma = w[0], C = w[1];//changed
  return (r * r / (2 * C)) * pow(1 + C * r * r, - gamma) - (pow(1 + C * r * r, 1 - gamma) - 1) / (2 * C * C * (1 - gamma));//changed
} */
// derivative of fr with respect to q
/* double dq_fr(double r, double w[])
{
  double gamma = w[0], D = w[1], q = w[2], mag = w[3], sig = D*exp( gamma*mag );
  return pow(1 + r * r / sig, 1 - q) * log(1 + r * r / sig) / (2 * PI);
} */

// *******************************************************************************
// compute lambda in Eq 2
double clambdaj(double *theta,
               int j,
               double *t,
               double *x,
               double *y,
/*                double *m, */
               double *bk)
{
  // extract model parameters
  double mu = theta[0] * theta[0],
    A     = theta[1] * theta[1],
    alpha = theta[2] * theta[2];
/*     gamma = theta[3] * theta[3]; */


  double part1, part2, part3, delta, r2;
   // first term of Eq 1
  double s = mu * bk[j];

  for (register int i = 0; i < j; i++)
    {
	// exp term of K function Eq 2
      /* part1 = exp(alpha * m[i]); */
	  part1 = A; //changed
// delta t in g function Eq 2
      delta = t[j] - t[i];
	  // g function Eq 2
      part2 = alpha * exp(- alpha * delta); //changed
/* // denominator of f function Eq 2
      sig = D * exp(gamma * m[i]); */
	  // numerator the power term of f function Eq 2
/*       r2 = dist2(x[j], y[j], x[i], y[i]); */
	  // f function Eq 2
/*       part3 = exp(- r2 / (2 * gamma * gamma)) / (2 * PI * gamma * gamma); //changed */
// Eq 2, summation part added to s in line 55
      /* s += A * part1 * part2 * part3; */
	  s += part1 * part2; //changed
    }
  return s;
}

// *******************************************************************************
// Compute lambda in eq 2 and the diferentation of lamba with respect
// to each parameter
void clambdajGr(double *theta,
               int j,
               double *t,
               double *x,
               double *y,
               /* double *m, */
               double *bk,
	       double *fv,
	       double *dfv)
{
  // extract model parameters
  double mu = theta[0] * theta[0],
    A     = theta[1] * theta[1],
    alpha = theta[2] * theta[2];
/*     gamma = theta[3] * theta[3]; */

  double part1, part2, part3, part2_alpha, part3_gamma,
    delta, r2, sg1, sg2 = 0, sg3 = 0,
    sg4 = 0, sg5 = 0, sg6 = 0, sg7 = 0;

// first term of Eq 1
  double s = mu * bk[j];
  // derivative of first term with respect to mu
      sg1 = bk[j];

  for (register int i = 0; i < j; i++)
    {
      /* part1 = exp(alpha * m[i]); */
	  part1 = A;

      delta = t[j] - t[i];
	  // function g
	  part2 = alpha * exp(- alpha * delta); //changed

/*       r2 = dist2(x[j], y[j], x[i], y[i]); */
	  // function f
/*       part3 = exp(- r2 / (2 * gamma * gamma)) / (2 * PI * gamma * gamma); //changed */
// Eq 2, summation part added to s in line 105
      /* s    += A * part1 * part2 * part3; */
	  s    += part1 * part2;//changed
	  
	  // derivative of s with respect to A
/*       sg2  += part1 * part2 * part3; */
	  sg2  += part2;//changed
	  
// derivative of part 2 and s with respect to alpha
      part2_alpha = (1 - alpha * delta) * exp(- alpha * delta);
      sg3    += part1 * part2_alpha;//changed
	  
// derivative of part 3 and s with respect to gamma
/*       part3_gamma = exp(- r2 / (2 * gamma * gamma)) * (r2 - 2 * gamma * gamma) / (2 * PI * gamma * gamma * gamma * gamma * gamma);
	  sg4    += part1 * part2 * part3_gamma; */
	  
    }
// Eq 2, summation part added to s in line 55
  *fv      = s;
  // differentiation of eq 2 with respect of each of the 8 parameters
  // Since they is a sqrt on the parameters you need the coeff 2*theta
  dfv[ 0 ] = sg1 * 2 * theta[0];
  dfv[ 1 ] = sg2 * 2 * theta[1];
  dfv[ 2 ] = sg3 * 2 * theta[2];
/*   dfv[ 3 ] = sg4 * 2 * theta[3]; */
}

// *******************************************************************************
//Compute terms inside summation of eauation 6
double cintegj(double *theta,
              int j,
              double *t,
              double *x,
              double *y,
             /*  double *m, */
              int *np,
              double *px,
              double *py,
              double *tstart2,
              double *tlength)
{

  // extract model parameters
  double //mu = theta[0] * theta[0],
    A     = theta[1] * theta[1],
    alpha = theta[2] * theta[2];
/*     gamma = theta[3] * theta[3]; */

  double ttemp, ttemp1, ttemp2, gi, gi1, gi2, w[2], si, sk;//changed w from 4 to 2

// time integral of g
  if (t[j] > *tstart2)
    {
      ttemp = *tlength - t[j];
      gi  = 1 - exp(- alpha * ttemp);//changed
    }
  else
    {
      ttemp1 = *tstart2 - t[j];
      ttemp2 = *tlength - t[j];

      gi1  = - exp(- alpha * ttemp1);//changed 
      gi2  = - exp(- alpha * ttemp2);//changed 
      gi   = gi2 - gi1;
    }

 /*  w[ 0 ] = gamma; */
/*   w[ 1 ] = C; //changed */
  /* w[ 2 ] = q; */
  /* w[ 3 ] = m[j]; */

//double integral of f in equation 6
/*   si = polyinteg(fr, w, np, px, py, x[j], y[j]); */
  
// Calculation of K  
/*   sk = A * exp(alpha * m[j]); */
  sk = A;//changed 
  return sk * gi;//changed
}

// *******************************************************************************
 // find the triple integration of the second term in eq 6
 // and its derivation with respect of each of 8 parameter
void cintegjGr(double *theta,
              int j,
              double *t,
              double *x,
              double *y,
              /* double *m, */
              int *np,
              double *px,
              double *py,
              double *tstart2,
              double *tlength,
              double *fv,
	      double *dfv)
{
  // extract model parameters
  double //mu = theta[0] * theta[0],
    A     = theta[1] * theta[1],
    alpha = theta[2] * theta[2];
/*     gamma = theta[3] * theta[3]; */

  double ttemp, ttemp1, ttemp2, gi, gi1, gi2, gialpha, gialpha1, gialpha2,
    w[2], si, sigamma, sk;//changed w from 4 to 2

  if (t[j] > *tstart2)
    {
      ttemp = *tlength - t[j];

// time integral of g function in eq 6
      gi  = 1 - exp(- alpha * ttemp);//changed
// derivative of gi with respect to alpha	  
      gialpha = ttemp * exp(- alpha * ttemp);//changed
    }
  else
    {
      ttemp1 = *tstart2 - t[j];
      ttemp2 = *tlength - t[j];
// integral of g function in eq 6 the two terms separately
      gi1  = - exp(- alpha * ttemp1);//changed 
      gi2  = - exp(- alpha * ttemp2);//changed 
// derivative of gi with respect to B and niu the two terms separately
      gialpha1 = ttemp1 * exp(- alpha * ttemp1);//changed
      gialpha2 = ttemp2 * exp(- alpha * ttemp2);//changed
// Subtracting the two terms
      gi  = gi2 - gi1;
      gialpha = gialpha2 - gialpha1;
    }

/*   w[0] = gamma; */
/*   w[1] = C; */
/*   w[2] = q;
  w[3] = m[j]; */
// double integral of f in eq 6 
/*   si      = polyinteg(fr, w, np, px, py, x[j], y[j]);
  // derivative of si with respect to gamma
  sigamma = polyinteg(dgamma_fr, w, np, px, py, x[j], y[j]);//changed */

// function k in eq 6
  /* sk = A * exp(alpha * m[j]); */
  sk = A;
  // second term in eq 6
  *fv      = sk * gi;//changed
  // derivation with respect of each parameter of second term in eq 6
  dfv[ 0 ] = 0;
  dfv[ 1 ] = sk * gi / A        * 2 * theta[1];
  dfv[ 2 ] = sk * gialpha       * 2 * theta[2];
/*   dfv[ 3 ] = sk * gi * sigamma        * 2 * theta[3]; */
  return;
}

// *******************************************************************************

SEXP clambdax(SEXP rt,
             SEXP rx,
             SEXP ry,
             SEXP theta,
             SEXP revents)
{
  SEXP dim;

  // extract events
  PROTECT(dim = allocVector(INTSXP, 2));
  dim = getAttrib(revents, R_DimSymbol);
  int N = INTEGER(dim)[0];
  double *events = REAL(revents);
  double t[N], x[N], y[N], m[N];
  for (int i = 0; i < N; i++)
    {
      t[i] = events[i];
      x[i] = events[N + i];
      y[i] = events[2 * N + i];
      m[i] = events[3 * N + i];
    }

  // extract model parameters
  double *tht = REAL(theta);
  double //mu = tht[0] * tht[0],
    A     = tht[1] * tht[1],
    c     = tht[2] * tht[2],
    alpha = tht[3] * tht[3],
    p     = tht[4] * tht[4],
    D     = tht[5] * tht[5],
    q     = tht[6] * tht[6],
    gamma = tht[7] * tht[7];

  // extract arguments
  double tt = *REAL(rt), xx = *REAL(rx), yy = *REAL(ry);

  double part1, part2, part3, delta, sig, r2;

  double s = 0;

  int i = 0;
  while (t[i] < tt && i < N)
    {
      part1 = exp(alpha * m[i]);

      delta = tt - t[i];
      part2 = (p - 1)/c * pow(1 + delta/c, - p);

      sig = D * exp(gamma * m[i]);
      r2 = dist2(xx, yy, x[i], y[i]);
      part3 = (q - 1) / (sig * PI) * pow(1 + r2 / sig, - q);

      s += A * part1 * part2 * part3;
      i++;
    }

  SEXP out;
  PROTECT(out = allocVector(REALSXP, 1));
  double *outP = REAL(out);
  *outP = s;
  UNPROTECT(2);
  return out;
}

// *******************************************************************************

