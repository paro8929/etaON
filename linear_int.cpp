#include <math.h>

double linear_int(double x0,double x1,double y0, double y1, double x)
{
  double temp=0;
  double t=x1-x0;
  if (isfinite(1/t))
    temp=((x-x0)*y1+(x1-x)*y0)/t;
  else
    temp=y0;
  return temp;
}


double quadratic_int(double x0,double x1,double x2, double z0, double z1, double z2, double x)
{
  double temp=0;
  double t1=x1-x0;
  double t2=x2-x1;
  double t3=x2-x0;
  double a,b,c;

  if (isfinite(1/t1)&&isfinite(1/t2)&&isfinite(1/t3))
    {
      a=(-x1*z0 + x2*z0 + x0*z1 - x2*z1 - x0*z2 + x1*z2)/(t1*t2*t3);
      b=(x1*x1*z0 - x2*x2*z0 - x0*x0*z1 + x2*x2*z1 + x0*x0*z2 - x1*x1*z2)/(t1*t2*t3);
      c=(-x1*x1*x2*z0 + x1*x2*x2*z0 + x0*x0*x2*z1 - x0*x2*x2*z1 - x0*x0*x1*z2 +  x0*x1*x1*z2)/(t1*t2*t3);   
      temp=a*x*x+b*x+c;
    }
  //else
      //printf("Error quadratic_int: cannot handle degenerate points\n");

  return temp;
}

//inspired by numerical recipes
//x0,x1: grid points in x-direction
//y0,y1: grid points in y-direction
//f0-f3: function value starting at x0,y0, continue counterclockwise
//put differently: f0=f(x0,y0)
//f1=f(x1,y0)
//f2=f(x1,y1)
//f3=f(x0,y1)
double bilinear_int(double x0,double x1,double y0, double y1, double f0, double f1, double f2, double f3, double x, double y)
{
  double temp=0;
  double t=(x-x0)/(x1-x0);
  double u=(y-y0)/(y1-y0);


  if ((isfinite(u)==1)&&(isfinite(t)==1))
    {
      temp=(1-t)*(1-u)*f0+t*(1-u)*f1+t*u*f2+(1-t)*u*f3;
    }
  else
    {
      if (isfinite(u)==1)
	temp=linear_int(y0,y1,f0,f2,y);	
      else if (isfinite(t)==1)
	temp=linear_int(x0,x1,f0,f2,x);
      else
	  temp=f0;
    }
  
  

  //printf("t=%f u=%f f0=%f f1=%f f2=%f f3=%f temp=%f aha %i\n",t,u,f0,f1,f2,f3,temp,isfinite(t));
  return temp;
}


double trilinear_int(double x0,double x1,double y0, double y1, double z0, double z1,double f000, double f001, double f010, double f011, double f100, double f101,double f110, double f111, double x, double y,double z)
{
  double temp=0;
  double t=(x-x0)/(x1-x0);
  double u=(y-y0)/(y1-y0);
  double v=(z-z0)/(z1-z0);

  int mysum=isfinite(u)+isfinite(t)+isfinite(v);

  if (mysum==3)
    {
      temp=(1-t)*(1-u)*(1-v)*f000;
      temp+=(1-t)*(1-u)*v*f001;
      temp+=(1-t)*u*(1-v)*f010;
      temp+=(1-t)*u*v*f011;
      temp+=t*(1-u)*(1-v)*f100;
      temp+=t*(1-u)*v*f101;
      temp+=t*u*(1-v)*f110;
      temp+=t*u*v*f111;      
    }
  if (mysum==2)
    {
      if ((isfinite(u))&&(isfinite(v)))
	temp=bilinear_int(y0,y1,z0,z1,f000,f010,f011,f001,y,z);
      if((isfinite(u))&&(isfinite(t)))
	temp=bilinear_int(x0,x1,y0,y1,f000,f100,f110,f010,x,y);
      if ((isfinite(t))&&(isfinite(v)))
	temp=bilinear_int(x0,x1,z0,z1,f000,f100,f101,f001,x,z);
    }
  if (mysum==1)
    {
      if (isfinite(t))
	temp=linear_int(x0,x1,f000,f100,x);
      if (isfinite(u))
	temp=linear_int(y0,y1,f010,f010,y);
      if (isfinite(v))
	temp=linear_int(z0,z1,f001,f001,z);
    }
  if (mysum==0)
    temp=f000;

  //printf("t=%f u=%f f0=%f f1=%f f2=%f f3=%f temp=%f aha %i\n",t,u,f0,f1,f2,f3,temp,isfinite(t));
  return temp;
}
