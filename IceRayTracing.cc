/*
This is IceRayTracing made for RadioScatter. Author: Uzair Latif 
released under GPL3.
*/

#include "IceRayTracing.hh"

/* Get the value of the B parameter for the refractive index model */
double IceRayTracing::GetB(double z){
  z=fabs(z);
  double B=0;

  B=-0.43;
  return B;
}

/* Get the value of the C parameter for the refractive index model */
double IceRayTracing::GetC(double z){
  z=fabs(z);
  double C=0;
  
  C=0.0132;
  return C;
}

/* Get the value of refractive index model for a given depth  */
double IceRayTracing::Getnz(double z){
  z=fabs(z);
  return IceRayTracing::A_ice+GetB(z)*exp(-GetC(z)*z);
}

/* E-feild Power Fresnel coefficient for S-polarised wave which is perpendicular to the plane of propogation/incidence. This function gives you back the reflectance. The transmittance is T=1-R */
double IceRayTracing::Refl_S(double thetai){
  double Nair=1;
  double Nice=IceRayTracing::Getnz(0); 
  double n1=Nice;
  double n2=Nair;
  
  double sqterm=sqrt(1-pow(1-(n1/n2)*(sin(thetai)),2));
  double num=n1*cos(thetai)-n2*sqterm;
  double den=n1*cos(thetai)+n2*sqterm;
  double RS=(num*num)/(den*den);
  return (RS);
}

/* E-feild Power Fresnel coefficient for P-polarised wave which is parallel to the plane of propogation/incidence. This function gives you back the reflectance. The transmittance is T=1-R */
double IceRayTracing::Refl_P(double thetai){
  double Nair=1;
  double Nice=IceRayTracing::Getnz(0); 
  double n1=Nice;
  double n2=Nair;
  
  double sqterm=sqrt(1-pow(1-(n1/n2)*(sin(thetai)),2));
  double num=n1*sqterm-n2*cos(thetai);
  double den=n1*sqterm+n2*cos(thetai);
  double RP=(num*num)/(den*den);
  return (RP);
}


/* Use GSL minimiser which uses Brent's Method to find root for a given function. This will be used to find roots wherever it is needed in my code.  */
double IceRayTracing::FindFunctionRoot(gsl_function F,double x_lo, double x_hi)
{
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_set_error_handler_off();
  status = gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  // cout<<x_lo<<" "<<x_hi<<" "<<endl;
  // printf ("error: %s\n", gsl_strerror (status));
  
  // printf ("using %s method\n", gsl_root_fsolver_name (s));
  // printf ("%5s [%9s, %9s] %9s %9s\n","iter", "lower", "upper", "root", "err(est)");

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,0, 0.000001);
	
      if (status == GSL_SUCCESS){
	// printf ("Converged:");
	// printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi,r,x_hi - x_lo);
      }
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  return r;
}

/* Define the function that will be minimized to get the value of the depth of the turning point for a given refracted ray. This function basically requires the value of the L parameter to find the depth of the turning point.  This comes from the part of the fDnfR function where sqrt( n(z) - L ). This imposes the constraint then n(z)=L at the turning point of the ray from which we can find zmax. */
double IceRayTracing::GetMinnz(double x,void *params){
  struct IceRayTracing::Minnz_params *p= (struct IceRayTracing::Minnz_params *) params;
  double A = p->a;
  double L = p->l;
  return A+GetB(x)*exp(-GetC(x)*x)-L;
}

/* Get the value of the depth of the turning point for the refracted ray */
double IceRayTracing::GetZmax(double A, double L){
  gsl_function F1;
  struct IceRayTracing::Minnz_params params1= {A,L};
  F1.function = &GetMinnz;
  F1.params = &params1;
  double zmax=FindFunctionRoot(F1,0.0,5000);
  return zmax;
}

/* Analytical solution describing ray paths in ice as function of depth */
double IceRayTracing::fDnfR(double x,void *params){
  
  struct IceRayTracing::fDnfR_params *p= (struct IceRayTracing::fDnfR_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double L = p->l;
  
  return (L/C)*(1.0/sqrt(A*A-L*L))*(C*x-log(A*Getnz(x)-L*L+sqrt(A*A-L*L)*sqrt(pow(Getnz(x),2)-L*L)));;
}

/* Analytical solution describing the ray path in ice as a function of the L parameter */
double IceRayTracing::fDnfR_L(double x,void *params){
  
  struct IceRayTracing::fDnfR_L_params *p= (struct IceRayTracing::fDnfR_L_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double Z = p->z;
  
  double result=(x/C)*(1.0/sqrt(A*A-x*x))*(C*Z-log(A*Getnz(Z)-x*x+sqrt(A*A-x*x)*sqrt(pow(Getnz(Z),2)-x*x)));
  return result;
}

/* The function used to calculate ray propogation time in ice */
double IceRayTracing::ftimeD(double x,void *params){

  struct IceRayTracing::ftimeD_params *p= (struct IceRayTracing::ftimeD_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double Speedc = p->speedc;
  double L = p->l;
  
  return (1.0/(Speedc*C*sqrt(pow(Getnz(x),2)-L*L)))*(pow(Getnz(x),2)-L*L+(C*x-log(A*Getnz(x)-L*L+sqrt(A*A-L*L)*sqrt(pow(Getnz(x),2)-L*L)))*(A*A*sqrt(pow(Getnz(x),2)-L*L))/sqrt(A*A-L*L) +A*sqrt(pow(Getnz(x),2)-L*L)*log(Getnz(x)+sqrt(pow(Getnz(x),2)-L*L)) );
}

/* This function is minimised to find the launch angle (or the L parameter) for the direct ray */
double IceRayTracing::fDa(double x,void *params){
  struct IceRayTracing::fDanfRa_params *p= (struct IceRayTracing::fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;

  struct IceRayTracing::fDnfR_L_params params1a = {A, GetB(z1), GetC(z1), z1};
  struct IceRayTracing::fDnfR_L_params params1b = {A, GetB(z0), GetC(z0), z0};
  
  return fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b) - x1;
}

/* This function is minimised to find the launch angle (or the L parameter) for the reflected ray */
double IceRayTracing::fRa(double x,void *params){
  struct IceRayTracing::fDanfRa_params *p= (struct IceRayTracing::fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;

  struct IceRayTracing::fDnfR_L_params params1a = {A, GetB(z1), -GetC(z1), -z1};
  struct IceRayTracing::fDnfR_L_params params1b = {A, GetB(z0), -GetC(z0), -z0};
  struct IceRayTracing::fDnfR_L_params params1c = {A, GetB(0.0000001), -GetC(0.0000001), 0.0000001};

  return fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b) - 2*( fDnfR_L(x,&params1c) - fDnfR_L(x,&params1b) ) - x1;
}

/* This function is minimised to find the launch angle (or the L parameter) for the refracted ray */
double IceRayTracing::fRaa(double x,void *params){
  struct IceRayTracing::fDanfRa_params *p= (struct IceRayTracing::fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;

  double zmax= GetZmax(A,x)+0.0000001;

  struct IceRayTracing::fDnfR_L_params params1a = {A, GetB(z1), -GetC(z1), -z1};
  struct IceRayTracing::fDnfR_L_params params1b = {A, GetB(z0), -GetC(z0), -z0};
  struct IceRayTracing::fDnfR_L_params params1c = {A, GetB(zmax), -GetC(zmax), zmax};

  return fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b) - 2*( fDnfR_L(x,&params1c) - fDnfR_L(x,&params1b) ) - x1;
}

/* This functions works for the Direct ray and gives you back the launch angle, receive angle and propagation time of the ray together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. */
double* IceRayTracing::GetDirectRayPar(double z0, double x1, double z1){

  double *output=new double[5];
  
  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* First we setup the fDa function that will be minimised to get the launch angle (or the L parameter) for the direct ray. */
  gsl_function F1;
  struct IceRayTracing::fDanfRa_params params1= {IceRayTracing::A_ice, z0, x1, z1};
  F1.function = &fDa;
  F1.params = &params1;

  /* In my raytracing solution given in the function fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0 and z1 and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */ 
  double UpLimnz[]={Getnz(z1),Getnz(z0)};
  double *UpperLimitL=min_element(UpLimnz,UpLimnz+2);

  /* Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function. */
  double lvalueD=FindFunctionRoot(F1,0.0000001,UpperLimitL[0]);
  double LangD=asin(lvalueD/Getnz(z0))*(180.0/IceRayTracing::pi);
  double checkzeroD=fDa(lvalueD,&params1);

  /* Get the propagation time for the direct ray using the ftimeD function after we have gotten the value of the L parameter. */
  struct IceRayTracing::ftimeD_params params2a = {IceRayTracing::A_ice, GetB(z0), -GetC(z0), IceRayTracing::c_light_ms,lvalueD};
  struct IceRayTracing::ftimeD_params params2b = {IceRayTracing::A_ice, GetB(z1), -GetC(z1), IceRayTracing::c_light_ms,lvalueD};

  /* we do the subtraction because we are measuring the time taken between the Tx and Rx positions */
  double timeD=+ftimeD(-z0,&params2a) - ftimeD(-z1,&params2b);

  /* Setup the function that will be used to calculate the angle of reception for all the rays */
  gsl_function F5;
  struct IceRayTracing::fDnfR_params params5a = {IceRayTracing::A_ice, GetB(z1), -GetC(z1), lvalueD};
  double result, abserr;
  F5.function = &fDnfR;

  /* Calculate the recieve angle for direc rays by calculating the derivative of the function at the Rx position */
  F5.params = &params5a;
  gsl_deriv_central (&F5, -z1, 1e-8, &result, &abserr);
  double RangD=atan(result)*(180.0/IceRayTracing::pi);

  /* When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line. */
  if(z1==z0 && isnan(RangD)==true){
    RangD=180-LangD;
  }
  
  /* This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation. */
  if(z1!=z0 && isnan(RangD)==true){
    RangD=90;
  }
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangD;
  output[1]=LangD;
  output[2]=timeD;
  output[3]=lvalueD;
  output[4]=checkzeroD;

  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangD;
    output[1]=180-RangD;
  }
  
  return output;
}

/* This functions works for the Reflected ray and gives you back the launch angle, receive angle and propagation times (of the whole ray and the two direct rays that make it up) together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. */
double *IceRayTracing::GetReflectedRayPar(double z0, double x1 ,double z1){

  double *output=new double[8];

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* First we setup the fRa function that will be minimised to get the launch angle (or the L parameter) for the reflected ray. */
  gsl_function F3;
  struct IceRayTracing::fDanfRa_params params3= {IceRayTracing::A_ice, z0, x1, z1};
  F3.function = &fRa;
  F3.params = &params3;

  /* In my raytracing solution given in the function fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0, z1 and also 0.0000001 m and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */
  double UpLimnz[]={Getnz(z1),Getnz(z0),Getnz(0.0000001)};
  double *UpperLimitL=min_element(UpLimnz,UpLimnz+3);
  
  /* Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fRa function. */
  double lvalueR=FindFunctionRoot(F3,0.0000001,UpperLimitL[0]);
  double LangR=asin(lvalueR/Getnz(z0))*(180.0/IceRayTracing::pi);
  double checkzeroR=fRa(lvalueR,&params3); 

  /* Get the propagation time for the reflected ray using the ftimeD function after we have gotten the value of the L parameter. */
  struct IceRayTracing::ftimeD_params params3a = {IceRayTracing::A_ice, GetB(z0), GetC(z0), IceRayTracing::c_light_ms,lvalueR};
  struct IceRayTracing::ftimeD_params params3b = {IceRayTracing::A_ice, GetB(z1), GetC(z1), IceRayTracing::c_light_ms,lvalueR};
  struct IceRayTracing::ftimeD_params params3c = {IceRayTracing::A_ice, GetB(0.0000001), GetC(0.0000001), IceRayTracing::c_light_ms,lvalueR};

  /* we do the subtraction because we are measuring the time taken between the Tx and Rx positions. In the reflected case we basically have two direct rays 1) from Tx to surface 2) from surface to Rx. */
  double timeR= 2*ftimeD(-0.0000001,&params3c) - ftimeD(z0,&params3a) - ftimeD(z1,&params3b);

  /* Also get the time for the two individual direct rays separately */
  double timeR1=ftimeD(-0.0000001,&params3c) - ftimeD(z0,&params3a);
  double timeR2=ftimeD(-0.0000001,&params3c) - ftimeD(z1,&params3b);

  /* flip the times back if the original positions were flipped */
  if(Flip==true){
    double dumR=timeR2;
    timeR2=timeR1;
    timeR1=dumR;
  }
  timeR1=timeR1;
  timeR2=timeR2;

  /* Setup the function that will be used to calculate the angle of reception for all the rays */
  gsl_function F5;
  struct IceRayTracing::fDnfR_params params5b = {IceRayTracing::A_ice, GetB(z1), GetC(z1), lvalueR};
  double result, abserr;
  F5.function = &fDnfR;
  
  /* Calculate the recieve angle for reflected ray by calculating the derivative of the function at the Rx position */
  F5.params = &params5b;
  gsl_deriv_central (&F5, z1, 1e-8, &result, &abserr);
  double RangR=180-atan(result)*(180.0/IceRayTracing::pi);

  /* When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line. */
  if(z1==z0 && isnan(RangR)==true){
    RangR=180-LangR;
  }

  /* This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation. */
  if(z1!=z0 && isnan(RangR)==true){
    RangR=90;
  }

  /* Calculate the angle of incidence of the reflected ray at the surface ice. This will be used to calculate the Fresnel Coefficients. The angle is calculated by calculating the derivative of the ray path fucnction at the surface*/
  struct IceRayTracing::fDnfR_params paramsIAngB = {A_ice, GetB(0.0000001), GetC(0.0000001), lvalueR};
  F5.function = &fDnfR; 
  F5.params = &paramsIAngB;
  gsl_deriv_central (&F5, -0.0000001, 1e-8, &result, &abserr);
  double IncidenceAngleInIce=atan(result)*(180.0/pi);

  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangR;
  output[1]=LangR;
  output[2]=timeR;
  output[3]=lvalueR;
  output[4]=checkzeroR;
  output[5]=timeR1;
  output[6]=timeR2;
  output[7]=IncidenceAngleInIce;
  
  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangR;
    output[1]=180-RangR;
  } 
  
  return output;
}

/* This functions works for the Refracted ray and gives you back the launch angle, receive angle and propagation times (of the whole ray and the two direct rays that make it up) together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. It requires the launch angle of the reflected ray as an input. */
double *IceRayTracing::GetRefractedRayPar(double z0, double x1 ,double z1, double LangR, double RangR){

  double *output=new double[8];

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }

  dsw=0;
  if(Flip==true){
    dsw=180-LangR;
    LangR=180-RangR;
    RangR=dsw;
  }
  
  /* Set up all the variables that will be used to get the parameters for the refracted ray */
  double lvalueR=sin(LangR*(IceRayTracing::pi/180))*Getnz(z0);
  double lvalueRa=0;
  double LangRa=0;
  double checkzeroRa=1000;

  double timeRa=0;
  double timeRa1=0;
  double timeRa2=0;
  double raytime=0;
  double RangRa=0;
  double zmax=10;

  /* In my raytracing solution given in the function fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0, z1 and also 0.0000001 m and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */
  double UpLimnz[]={Getnz(z1),Getnz(z0)};
  double *UpperLimitL=min_element(UpLimnz,UpLimnz+2);
  
  /* First we setup the fRa function that will be minimised to get the launch angle (or the L parameter) for the refracted ray. */
  gsl_function F4;
  struct IceRayTracing::fDanfRa_params params4= {IceRayTracing::A_ice, z0, x1, z1};
  F4.function = &fRaa;
  F4.params = &params4;

  /* Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fRaa function. The thing to note here is the lower limit of the minimisation function is set to the L value corresponding to the reflected ray launch angle. Since we know the refracted ray always has bigger launch angle the reflected ray this reduces our range and makes the function more efficient at finding the refracted ray launch angle. */
  lvalueRa=FindFunctionRoot(F4,Getnz(z0)*sin((LangR*(IceRayTracing::pi/180.0))),UpperLimitL[0]);
  LangRa=asin(lvalueRa/Getnz(z0))*(180.0/IceRayTracing::pi);
  checkzeroRa=fRaa(lvalueRa,&params4);

  /* If the above strategy did not work then we start decreasing the reflected ray launch angle in steps of 5 degree and increase our range for minimisation to find the launch angle (or the L parameter). Sometimes the refracted and reflected rays come out to be the same in that case also I forced my solution to try harder by changing the minimisation range. */
  double iangstep=5;
  while( (isnan(checkzeroRa)==true || fabs(checkzeroRa)>0.5 || fabs(lvalueRa-lvalueR)<pow(10,-5) || fabs(LangRa-LangR)<pow(10,-1)) && LangR>iangstep && iangstep<90){
    //cout<<"2nd try to get Refracted ray "<<isnan(checkzeroRa)<<" "<<fabs(checkzeroRa)<<endl;
    lvalueRa=FindFunctionRoot(F4,Getnz(z0)*sin(((LangR-iangstep)*(IceRayTracing::pi/180.0))),UpperLimitL[0]);
    LangRa=asin(lvalueRa/Getnz(z0))*(180.0/IceRayTracing::pi);
    checkzeroRa=fRaa(lvalueRa,&params4);
    iangstep=iangstep+5;
  }///end the second attempt    

  /* If we still did not find a refracted ray then set the check zero parameter to 1000 to make sure my code does not output this as a possible solution */
  if(isnan(checkzeroRa)==true){
    checkzeroRa=1000;
  }

  /* If we did find a possible refracted ray then now we need to find the depth at which the ray turns back down without hitting the surface. */
  zmax=GetZmax(IceRayTracing::A_ice,lvalueRa)+0.0000001;

  /* If the turning point depth also came out to be zero then now we are sure that there is no refracted ray */
  if(zmax==0.0000001){
    checkzeroRa=1000;
  }

  /* Set parameters for ftimeD function to get the propagation time for the refracted ray */
  struct IceRayTracing::ftimeD_params params4a = {IceRayTracing::A_ice, GetB(z0), GetC(z0), IceRayTracing::c_light_ms,lvalueRa};
  struct IceRayTracing::ftimeD_params params4b = {IceRayTracing::A_ice, GetB(z1), GetC(z1), IceRayTracing::c_light_ms,lvalueRa};
  struct IceRayTracing::ftimeD_params params4c = {IceRayTracing::A_ice, GetB(zmax), GetC(zmax), IceRayTracing::c_light_ms,lvalueRa};

  /* This if condition checks if the function has not gone crazy and given us a turning point of the ray which is lower than both Tx and Rx and is shallower in depth than both */
  if((z0<-zmax || zmax<-z1)){
    /* we do the subtraction because we are measuring the time taken between the Tx and Rx positions. In the refracted case we basically have two direct rays 1) from Tx to turning point 2) from turning point to Rx. */
    raytime=2*ftimeD(zmax,&params4c) - ftimeD(z0,&params4a) - ftimeD(z1,&params4b);

    /* Also get the time for the two individual direct rays separately */
    timeRa1=ftimeD(zmax,&params4c) - ftimeD(z0,&params4a);
    timeRa2=ftimeD(zmax,&params4c) - ftimeD(z1,&params4b);
    if(Flip==true){
      double dumRa=timeRa2;
      timeRa2=timeRa1;
      timeRa1=dumRa;
    }
    timeRa1=timeRa1;
    timeRa2=timeRa2;
  }
  timeRa=raytime;

  /* Setup the function that will be used to calculate the angle of reception for all the rays */
  gsl_function F5;
  struct IceRayTracing::fDnfR_params params5c = {IceRayTracing::A_ice, GetB(z1), GetC(z1), lvalueRa};
  double result, abserr;
  F5.function = &fDnfR;
    
  /* Calculate the recieve angle for refacted ray by calculating the derivative of the function at the Rx position */
  F5.params = &params5c;
  gsl_deriv_central (&F5, z1, 1e-8, &result, &abserr);
  RangRa=180-atan(result)*(180.0/IceRayTracing::pi);

  /* When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line. */
  if(z1==z0 && isnan(RangRa)==true){
    RangRa=180-LangRa;
  }

  /* This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation. */
  if(z1!=z0 && isnan(RangRa)==true){
    RangRa=90;
  }

  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangRa;
  output[1]=LangRa;
  output[2]=timeRa;
  output[3]=lvalueRa;
  output[4]=checkzeroRa;
  output[5]=timeRa1;
  output[6]=timeRa2;
  output[7]=zmax;

  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangRa;
    output[1]=180-RangRa;
  }
  
  return output;
}

/* This function returns the x and z values for the full Direct ray path in a TGraph and also prints out the ray path in a text file */
TGraph* IceRayTracing::GetFullDirectRayPath(double z0, double x1, double z1,double lvalueD){

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
   
  /* Set the name of the text files */
  ofstream aoutD("DirectRay.txt");
  /* Set the step size for plotting */
  double h=0.1;
  /* Set the total steps required for looping over the whole ray path */
  int dmax=100000;
  /* Set the values to start the rays from */
  double zn=z1;
  double xn=0;

  /* Map out the direct ray path */
  int npnt=0;
  double checknan=0;
  struct IceRayTracing::fDnfR_params params6a;
  struct IceRayTracing::fDnfR_params params6b;
  
  TGraph *gr1=new TGraph();
  for(int i=0;i<dmax;i++){
    params6a = {IceRayTracing::A_ice, GetB(zn), GetC(zn), lvalueD};
    params6b = {IceRayTracing::A_ice, GetB(z0), GetC(z0), lvalueD};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
    checknan=fDnfR(zn,&params6a);
    if(isnan(checknan)==false && Flip==false){
      gr1->SetPoint(npnt,xn,zn);
      aoutD<<npnt<<" "<<xn<<" "<<zn<<endl;;
      npnt++;
    }

    if(isnan(checknan)==false && Flip==true){
      gr1->SetPoint(npnt,x1-xn,zn);
      aoutD<<npnt<<" "<<x1-xn<<" "<<zn<<endl;;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;      
    }  
  }

  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  return gr1;

}

/* This function returns the x and z values for the full Reflected ray path in a TGraph and also prints out the ray path in a text file */
TGraph* IceRayTracing::GetFullReflectedRayPath(double z0, double x1, double z1,double lvalueR){

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* Set the name of the text files */
  ofstream aoutR("ReflectedRay.txt");
  /* Set the step size for plotting. */
  double h=0.1;
  /* Set the total steps required for looping over the whole ray path */
  int dmax=100000;
  /* Set the values to start the rays from */
  double zn=z1;
  double xn=0;
  
  /* Map out the direct ray path */
  int npnt=0;
  double checknan=0;  
  struct IceRayTracing::fDnfR_params params6a;
  struct IceRayTracing::fDnfR_params params6b;
  struct IceRayTracing::fDnfR_params params6c;

  /* Map out the 1st part of the reflected ray */
  TGraph *gr2=new TGraph();
  for(int i=0;i<dmax;i++){
    params6a = {IceRayTracing::A_ice, GetB(zn), -GetC(zn), lvalueR};
    params6b = {IceRayTracing::A_ice, GetB(z0), -GetC(z0), lvalueR};
    params6c = {IceRayTracing::A_ice, GetB(0.0000001), -GetC(0.0000001), lvalueR};
    xn=(fDnfR(-zn,&params6a)-fDnfR(-z0,&params6b)+2*fabs(fDnfR(0.0000001,&params6c)-fDnfR(-z0,&params6b)));
    checknan=fDnfR(-zn,&params6a);
    if(isnan(checknan)==false && zn<=0 && Flip==false){
      gr2->SetPoint(npnt,xn,zn);
      aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(isnan(checknan)==false && zn<=0 && Flip==true){
      gr2->SetPoint(npnt,x1-xn,zn);
      aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
      
    zn=zn+h;
    if(zn>0){
      i=dmax+2;      
    }
  }

  /* Map out the 2nd part of the reflected ray */
  zn=-0.0000001;
  for(int i=0;i<dmax;i++){  
    params6a = {IceRayTracing::A_ice, GetB(zn), GetC(zn), lvalueR};
    params6b = {IceRayTracing::A_ice, GetB(z0), GetC(z0), lvalueR};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
    checknan=fDnfR(zn,&params6a);
    if(isnan(checknan)==false && Flip==false){
      gr2->SetPoint(npnt,xn,zn);
      aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }
      
    if(isnan(checknan)==false && Flip==true){
      gr2->SetPoint(npnt,x1-xn,zn);
      aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;
    }
  }

  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  return gr2;

}

/* This function returns the x and z values for the full Refracted ray path in a TGraph and also prints out the ray path in a text file */
TGraph* IceRayTracing::GetFullRefractedRayPath(double z0, double x1, double z1, double zmax, double lvalueRa){

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* Set the name of the text files */
  ofstream aoutRa("RefractedRay.txt");
  /* Set the step size for plotting. */
  double h=0.1;
  /* Set the total steps required for looping over the whole ray path */
  int dmax=100000;
  /* Set the values to start the rays from */
  double zn=z1;
  double xn=0;

  /* Map out the direct ray path */
  int npnt=0;
  double checknan=0;
  struct IceRayTracing::fDnfR_params params6a;
  struct IceRayTracing::fDnfR_params params6b;
  struct IceRayTracing::fDnfR_params params6c;  

  /* Map out the 1st part of the refracted ray */
  TGraph *gr3=new TGraph();    
  for(int i=0;i<dmax;i++){
    params6a = {IceRayTracing::A_ice, GetB(zn), -GetC(zn), lvalueRa};
    params6b = {IceRayTracing::A_ice, GetB(z0), -GetC(z0), lvalueRa};
    params6c = {IceRayTracing::A_ice, GetB(zmax), -GetC(zmax), lvalueRa};
    xn=(fDnfR(-zn,&params6a)-fDnfR(-z0,&params6b)+2*fabs(fDnfR(zmax,&params6c)-fDnfR(-z0,&params6b)));
    checknan=fDnfR(-zn,&params6a);
    if(isnan(checknan)==false && zn<=0 && Flip==false){
      gr3->SetPoint(npnt,xn,zn);
      aoutRa<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(isnan(checknan)==false && zn<=0 && Flip==true){
      gr3->SetPoint(npnt,x1-xn,zn);
      aoutRa<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
    
    zn=zn+h;
    if(zn>-zmax){
      i=dmax+2;      
    }
  }

  /* Map out the 2nd part of the refracted ray */
  zn=-zmax;
  for(int i=0;i<dmax;i++){  
    params6a = {IceRayTracing::A_ice, GetB(zn), GetC(zn), lvalueRa};
    params6b = {IceRayTracing::A_ice, GetB(z0), GetC(z0), lvalueRa};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
    checknan=fDnfR(zn,&params6a);
    if(isnan(checknan)==false && Flip==false){
      gr3->SetPoint(npnt,xn,zn);
      aoutRa<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(isnan(checknan)==false && Flip==true){
      gr3->SetPoint(npnt,x1-xn,zn);
      aoutRa<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
    
    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;
    }
  }

  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  return gr3;
  
}

/* function for plotting and storing all the rays */
void IceRayTracing::PlotAndStoreRays(double x0,double z0, double z1, double x1, double zmax, double lvalues[3], double checkzeroes[3]){
  
  double lvalueD=lvalues[0];
  double lvalueR=lvalues[1];
  double lvalueRa=lvalues[2];

  double checkzeroD=checkzeroes[0];
  double checkzeroR=checkzeroes[1];
  double checkzeroRa=checkzeroes[2]; 

  TMultiGraph *mg=new TMultiGraph();
  
  TGraph *gr1=GetFullDirectRayPath(z0,x1,z1,lvalueD);
  TGraph *gr2=GetFullReflectedRayPath(z0,x1,z1,lvalueR);
  TGraph *gr3=new TGraph();
  if((fabs(checkzeroR)>0.5 || fabs(checkzeroD)>0.5) && fabs(checkzeroRa)<0.5){
    gr3=GetFullRefractedRayPath(z0,x1,z1,zmax,lvalueRa);
  }

  gr1->SetMarkerColor(kBlue);
  gr2->SetMarkerColor(kBlue);
  gr3->SetMarkerColor(kBlue);
  
  /* Plot the all the possible ray paths on the canvas */
  TGraph *gr4=new TGraph();
  gr4->SetPoint(0,x1,z1);
  gr4->SetMarkerStyle(20);
  gr4->SetMarkerColor(kRed);

  TGraph *gr4b=new TGraph();
  gr4b->SetPoint(0,0,z0);
  gr4b->SetMarkerStyle(20);
  gr4b->SetMarkerColor(kGreen);
    
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerColor(2);

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerColor(2);

  double zlower=z0;
  if(fabs(z0)<fabs(z1)){
    zlower=z1;
  }
  if(fabs(z0)>fabs(z1)){
    zlower=z0;
  }
  TGraph *gr5=new TGraph();
  gr5->SetPoint(0,0,zlower-50);
  gr5->SetPoint(1,x1+50,0);

  if(fabs(checkzeroD)<0.5){
    mg->Add(gr1);
  }
  if(fabs(checkzeroR)<0.5){
    mg->Add(gr2);
  }
  if(fabs(checkzeroRa)<0.5){
    mg->Add(gr3);
  }
  mg->Add(gr4);
  mg->Add(gr4b);
  //mg->Add(gr5);

  TString title="Depth vs Distance, Tx at x=";
  title+=x0;
  title+=" m,z=";
  title+=(int)z0;
  title+=" m, Rx at x=";
  title+=x1;
  title+=" m,z=";
  title+=(int)z1;
  title+=" m; Distance (m);Depth (m)";
  mg->SetTitle(title);
  
  TCanvas *cRay=new TCanvas("cRay","cRay");
  cRay->cd();
  mg->Draw("AP");
  mg->GetXaxis()->SetNdivisions(20);
  cRay->SetGridx();
  cRay->SetGridy();
}

/* This is the main raytracing function. x0 always has to be zero. z0 is the Tx depth in m and z1 is the depth of the Rx in m. Both depths are negative. x1 is the distance between them */
double *IceRayTracing::IceRayTracing(double x0, double z0, double x1, double z1){

  /* define a pointer to give back the output of raytracing */ 
  double *output=new double[12];

  /* Store the ray paths in text files */
  bool PlotRayPaths=false;
  /* calculate the attenuation (not included yet!) */
  bool attcal=false;
  
  double Txcor[2]={x0,z0};/* Tx positions */
  double Rxcor[2]={x1,z1};/* Rx Positions */
  
  /*  ********This part of the code will try to get the Direct ray between Rx and Tx.********** */
  double* GetDirectRay=GetDirectRayPar(z0,x1,z1);
  double RangD=GetDirectRay[0];
  double LangD=GetDirectRay[1];
  double timeD=GetDirectRay[2];
  double lvalueD=GetDirectRay[3];
  double checkzeroD=GetDirectRay[4];
  delete []GetDirectRay;
  
  /* ********This part of the code will try to get the Reflected ray between Rx and Tx.********** */
  double* GetReflectedRay=GetReflectedRayPar(z0,x1,z1);
  double RangR=GetReflectedRay[0];
  double LangR=GetReflectedRay[1];
  double timeR=GetReflectedRay[2];
  double lvalueR=GetReflectedRay[3];
  double checkzeroR=GetReflectedRay[4];
  double timeR1=GetReflectedRay[5];
  double timeR2=GetReflectedRay[6];
  double AngleOfIncidenceInIce=GetReflectedRay[7];
  delete []GetReflectedRay;

  /* ********This part of the code will try to get the Refracted ray between Rx and Tx.********** */
  double RangRa=0;
  double LangRa=0;
  double timeRa=0;
  double lvalueRa=0; 
  double checkzeroRa=0;
  double timeRa1=0;
  double timeRa2=0;
  double zmax=0;
  
  /* This if condition makes sure that we only try to find a refracted ray if we don't get two possible ray paths from the direct and reflected case. This saves us alot of time since we know that between each Tx and Rx position we only expect 2 rays. */
  if(fabs(checkzeroR)>0.5 || fabs(checkzeroD)>0.5){
    double* GetRefractedRay=GetRefractedRayPar(z0,x1,z1,LangR,RangR);
    RangRa=GetRefractedRay[0];
    LangRa=GetRefractedRay[1];
    timeRa=GetRefractedRay[2];
    lvalueRa=GetRefractedRay[3]; 
    checkzeroRa=GetRefractedRay[4];
    timeRa1=GetRefractedRay[5];
    timeRa2=GetRefractedRay[6];
    zmax=GetRefractedRay[7];
    delete []GetRefractedRay;
  }

  /* Fill in the output pointer after calculating all the results */
  output[0]=LangD;
  output[1]=LangR;
  output[2]=LangRa;
  output[3]=timeD;
  output[4]=timeR;
  output[5]=timeRa;
  output[6]=RangD;
  output[7]=RangR;
  output[8]=RangRa;

  /* This part of the code can be used if the user wants to plot the individual ray paths. This part of the code prints out the individual ray paths in text files and also plots them on a canvas */
  if(PlotRayPaths==true){
    double lvalues[3];
    lvalues[0]=lvalueD;
    lvalues[1]=lvalueR;
    lvalues[2]=lvalueRa;

    double checkzeroes[3];
    checkzeroes[0]=checkzeroD;
    checkzeroes[1]=checkzeroR;
    checkzeroes[2]=checkzeroRa;
    
    PlotAndStoreRays(x0,z0,z1,x1,zmax,lvalues,checkzeroes);
  }
  
  /* print out all the output from the code */
  //cout<<0<<" ,x0= "<<x0<<" ,z0= "<<z0<<" ,x1= "<<x1<<" ,z1= "<<z1<<" ,langRa= "<<output[2]<<" ,langR= "<<output[1]<<" ,langD= "<<output[0]<<" ,langD-langR= "<<output[0]-output[1]<<" ,langD-langRa= "<<output[0]-output[2]<<" ,RangRa= "<<output[8]<<" ,RangR= "<<output[7]<<" ,RangD= "<<output[6]<<" ,RangR-RangD= "<<output[7]-output[6]<<" ,RangRa-RangD= "<<output[8]-output[6]<<" ,timeRa= "<<output[5]<<" ,timeR= "<<output[4]<<" ,timeD= "<<output[3]<<" ,timeR-timeD= "<<output[4]-output[3]<<" ,timeRa-timeD= "<<output[5]-output[3]<<" ,lvalueRa "<<lvalueRa<<" ,lvalueR "<<lvalueR<<" "<<" ,lvalueD "<<lvalueD<<" ,checkzeroRa "<<checkzeroRa<<" ,checkzeroR "<<checkzeroR<<" ,checkzeroD "<<checkzeroD<<endl;

  /* fill in the output array part where you fill in the times for the two parts of the reflected or refracted rays */
  if(fabs(checkzeroR)<0.5){
    output[9]=timeR1;
    output[10]=timeR2;
  }
  
  if(fabs(checkzeroRa)<0.5){
    output[9]=timeRa1;
    output[10]=timeRa2;
  }

  output[11]=AngleOfIncidenceInIce;

  /* Set the recieve angle to be zero for a ray which did not give us a possible path between Tx and Rx. I use this as a flag to determine which two rays gave me possible ray paths. */
  if(fabs(checkzeroD)>0.5){
    output[6]=0;
  }
  if(fabs(checkzeroR)>0.5){
    output[7]=0;
  }
  if(fabs(checkzeroRa)>0.5){
    output[8]=0;
  }
  
  return output;
}
