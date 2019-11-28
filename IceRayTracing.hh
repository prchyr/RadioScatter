/*
This is IceRayTracing made for RadioScatter. Author: Uzair Latif 
released under GPL3.
*/
#ifndef IRT_HEAD
#define IRT_HEAD

#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "TAxis.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <sys/time.h>

using namespace std;

namespace IceRayTracing{

  /* Set the value of pi */
  static constexpr double pi=3.14159265359;
  /* Set the value of the asymptotic parameter of the refractive index model */ 
  static constexpr double A_ice=1.78;
  /* Set the value of the speed of light in m/s */ 
  static constexpr double c_light_ms=299792458;
  
  /* Get the value of the B parameter for the refractive index model */
  double GetB(double z);
  
  /* Get the value of the C parameter for the refractive index model */
  double GetC(double z);
  
  /* Get the value of refractive index model for a given depth  */
  double Getnz(double z);

  /* E-feild Power Fresnel coefficient for S-polarised wave which is perpendicular to the plane of propogation/incidence. This function gives you back the reflectance. The transmittance is T=1-R */
  double_t Refl_S(double thetai);

  /* E-feild Power Fresnel coefficient for P-polarised wave which is parallel to the plane of propogation/incidence. This function gives you back the reflectance. The transmittance is T=1-R */
  double Refl_P(double thetai);
   
  /* The temperature and attenuation model have been taken from AraSim which also took it from here http://icecube.wisc.edu/~araproject/radio/ . This is basically Matt Newcomb's icecube directory which has alot of information, plots and codes about South Pole Ice activities. Please read it if you find it interesting. */
  
  /* Temperature model:The model takes in value of depth z in m and returns the value of temperature in Celsius.*/
  double GetIceTemperature(double z);

  /* Ice Attenuation Length model: Takes in value of frequency in Ghz and depth z and returns you the value of attenuation length in m */
  double GetIceAttenuationLength(double z, double frequency);
 
  /* Use GSL minimiser which uses Brent's Method to find root for a given function. This will be used to find roots wherever it is needed in my code.  */
  double FindFunctionRoot(gsl_function F,double x_lo, double x_hi);

  /* Define the function that will be minimized to get the value of the depth of the turning point for a given refracted ray. This function basically requires the value of the L parameter to find the depth of the turning point.  This comes from the part of the fDnfR function where sqrt( n(z) - L ). This imposes the constraint then n(z)=L at the turning point of the ray from which we can find zmax. */
  struct Minnz_params { double a,l; };
  double GetMinnz(double x,void *params);
  
  /* Get the value of the depth of the turning point for the refracted ray */
  double GetZmax(double A, double L);
  
  /* Analytical solution describing ray paths in ice as function of depth */
  struct fDnfR_params { double a, b, c, l; };
  double fDnfR(double x,void *params);
  
  /* Analytical solution describing the ray path in ice as a function of the L parameter */
  struct fDnfR_L_params { double a, b, c, z; };
  double fDnfR_L(double x,void *params);
  
  /* The function used to calculate ray propogation time in ice */
  struct ftimeD_params { double a, b, c, speedc,l; };
  double ftimeD(double x,void *params);
  
  /* This function is minimised to find the launch angle (or the L parameter) for the direct ray */
  struct fDanfRa_params { double a, z0, x1, z1; };
  double fDa(double x,void *params);
  
  /* This function is minimised to find the launch angle (or the L parameter) for the reflected ray */
  double fRa(double x,void *params);

  /* This function is minimised to find the launch angle (or the L parameter) for the refracted ray */
  double fRaa(double x,void *params);

  /* This functions works for the Direct ray and gives you back the launch angle, receive angle and propagation time of the ray together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. */
  double* GetDirectRayPar(double z0, double x1, double z1);

  /* This functions works for the Reflected ray and gives you back the launch angle, receive angle and propagation times (of the whole ray and the two direct rays that make it up) together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. */
  double *GetReflectedRayPar(double z0, double x1 ,double z1);
  
  /* This functions works for the Refracted ray and gives you back the launch angle, receive angle and propagation times (of the whole ray and the two direct rays that make it up) together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. It requires the launch angle of the reflected ray as an input. */
  double *GetRefractedRayPar(double z0, double x1 ,double z1, double LangR,double RangR);
  
  /* This function returns the x and z values for the full Direct ray path in a TGraph and also prints out the ray path in a text file */
  TGraph* GetFullDirectRayPath(double z0, double x1, double z1,double lvalueD);

  /* This function returns the x and z values for the full Reflected ray path in a TGraph and also prints out the ray path in a text file */
  TGraph* GetFullReflectedRayPath(double z0, double x1, double z1,double lvalueR);

  /* This function returns the x and z values for the full Refracted ray path in a TGraph and also prints out the ray path in a text file */
  TGraph* GetFullRefractedRayPath(double z0, double x1, double z1, double zmax, double lvalueRa);

  /* function for plotting and storing all the rays */
  void PlotAndStoreRays(double x0,double z0, double z1, double x1, double zmax, double lvalues[3], double checkzeroes[3]);

  /* This is the main raytracing function. x0 always has to be zero. z0 is the Tx depth in m and z1 is the depth of the Rx in m. Both depths are negative. x1 is the distance between them */
  double *IceRayTracing(double x0, double z0, double x1, double z1);

}
#endif
