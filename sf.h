#ifndef INCLUDE__SF_H
#define INCLUDE__SF_H
#include <gsl/gsl_sf.h>
#include <cstdlib>
#include <cmath>
#include <complex>

using namespace std;

namespace SpherHarmonics{
  double legendre(int ell,double ctheta);
  complex<double> Ylm(int ell,int m,double theta,double phi);
  complex <double> Ylm(int ell, int m, double x, double y, double z);

  double norm_legendre(int ell, int em, double ctheta);

  void YlmUpToL(int lmax, double x, double y, double z, complex<double> *ylms);
  void YlmUpToL(int lmax, double ctheta, double phi, complex<double> *ylms);

  double ReYlm(int ell, int m, double theta, double phi);
  double ReYlm(int ell, int m, double x, double y, double z);
  double ImYlm(int ell, int m, double theta, double phi);
  double ImYlm(int ell, int m, double x, double y, double z);
  
};

#endif

