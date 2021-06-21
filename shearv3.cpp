#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <Eigen/Dense>
#include <linear_int.h>
#include <rootfinding.cpp>

using namespace Eigen;
using namespace std;



int D=2; //number of space dimensions
//int LL=100; //lambda parameter
int quad=4; //quadrature order

int ABORT=0;
int singlerun=0;
double lambda=1;
int vlevel=0;

//tech specs for integratos
int MAX=100000; //max number of integrator steps
int RM=40;//reduced number of integrator steps for phi-integration
int RM2=400;
int RM3=400;
double EPS=0.0002;//integrator spacing
double safe=0.01; //standard: use 0.025
double E2=10.; 


//for interpolation of Pi
int NN=80;
double **dataR1,**dataR2,**dataRe;
int SF3=400;
double **m3func;


#include<mathfuncs.cpp>
#include<pifuncs.cpp>
#include<quadrature.cpp>








double massfuncD2(double *args)
{
  double m=args[0];
  double lambda=args[5];
  //printf("this is massfunc for m=%f lambda=%f\n",m,lambda);
  double res=(m+2*log(1-exp(-m)))/M_PI+m*m/lambda;
  return res;
}

double massfuncD3(double *args)
{
  double m=args[0];
  double lambda=args[5];
  //printf("this is massfunc for m=%f lambda=%f\n",m,lambda);
  double res=m-8*eval_lin_int(m3func,m,SF3)/(4*M_PI*M_PI/lambda+1+2*log(1./m));
  return res;
}


double massfunc(double *args)
{
  double res;
  if (D==2)
    res=massfuncD2(args);
  else if (D==3)
    res=massfuncD3(args);
  else printf("ERROR D=%i not defined for mass function\n",D);

  return res;
}



double rhoD(double p0,double p,double m,double lambda)
{
  //printf("rhoD %.16g %.16g\n",p0,p);
  
  //double im=ImPi0(p0,p,m)+ImPiT(p0,p,m);
  //double re=RePi0(p0,p,m)+RePiT(p0,p,m);

  double im=ImPi0(p0,p,m)+ImInt(p0,p,m,dataR1,dataR2,NN);
  
  double re=RePi0(p0,p,m)+Reint(p0,p,dataRe,NN);
  re+=0.5/lambda;

  if (D==3)
    re+=log(1./m)/(4*M_PI*M_PI);

  if (isnan(im))
    {
      printf("im nan: %f where %f\n",im,ImInt(p0,p,m,dataR1,dataR2,NN));
    }
  if (isnan(re))
    {
      printf("re(p0=%f,p=%f,m=%f,lambda=%f) nan: %f where %f check %f\n",p0,p,m,lambda,re,Reint(p0,p,dataRe,NN),RePiT(p0,p,m));
    }
  
  //printf("%f %f\n",im,re);
  
  double res=-im/(im*im+re*re);
  
  
  
  
  return res;
}

double WrhoD(double p0,double p,double m,double lambda)
{
  //printf("rhoD %.16g %.16g\n",p0,p);
  
  //double im=ImPi0(p0,p,m)+ImPiT(p0,p,m);
  //double re=RePi0(p0,p,m)+RePiT(p0,p,m);

  double im=ImPi0(p0,p,m)+ImInt(p0,p,m,dataR1,dataR2,NN);
  
  
  
  double res=-im;
  
  
  
  
  return res;
}

double DR2(double p0,double p,double m,double lambda)
{
  

  double im=ImPi0(p0,p,m)+ImInt(p0,p,m,dataR1,dataR2,NN);
  
  double re=RePi0(p0,p,m)+Reint(p0,p,dataRe,NN);
  re+=0.5/lambda;

  if (D==3)
    re+=log(1./m)/(4*M_PI*M_PI);

  
  
  double res=1./(im*im+re*re);
  
  
  
  
  return res;
}

double intrhoDminus(double p,double k,double m,double lambda)
{
  double res=0;
  double meps=M_PI/RM;
  double epek=sqrt(p*p+m*m)-sqrt(k*k+m*m);
  double add=0;
  for (int i=0;i<2*RM;i++)
    {
      double phi=i*meps;
      if (D==2)
	add=rhoD(epek,sqrt(p*p+k*k-2*p*k*cos(phi)),m,lambda);
      else if (D==3)
	add=0.5*rhoD(epek,sqrt(p*p+k*k-2*p*k*cos(phi)),m,lambda)*fabs(sin(phi));
      //printf("phi=%f add=%f\n",phi,add);
      res+=add;      
    }
  res/=2*RM;
  return res;
}

double WintrhoDminus(double p,double k,double m,double lambda)
{
  double res=0;
  double meps=M_PI/RM;
  double epek=sqrt(p*p+m*m)-sqrt(k*k+m*m);
  double add=0;
  for (int i=0;i<2*RM;i++)
    {
      double phi=i*meps;
      if (D==2)
	add=WrhoD(epek,sqrt(p*p+k*k-2*p*k*cos(phi)),m,lambda);
      else if (D==3)
	add=0.5*WrhoD(epek,sqrt(p*p+k*k-2*p*k*cos(phi)),m,lambda)*fabs(sin(phi));
      //printf("phi=%f add=%f\n",phi,add);
      res+=add;      
    }
  res/=2*RM;
  return res;
}

double WintrhoDplus(double p,double k,double m,double lambda)
{
  //printf("D2 p=%f k=%f m=%f l=%f\n",p,k,m,lambda);
  double res=0;
  double meps=M_PI/RM;
  double epek=sqrt(p*p+m*m)+sqrt(k*k+m*m);
  double add=0;
  for (int i=0;i<2*RM;i++)
    {
      double phi=i*meps;
      if (D==2)
	add=WrhoD(epek,sqrt(p*p+k*k-2*p*k*cos(phi)),m,lambda);
      else if (D==3)
	add=0.5*WrhoD(epek,sqrt(p*p+k*k-2*p*k*cos(phi)),m,lambda)*fabs(sin(phi));
      //printf("phi=%f add=%f res=%f\n",phi,add,res);
      res+=add;
      if (isnan(res))
	{
	  printf("i=%i p=%f k=%f m=%f lambda=%f epek=%f sq=%f\n",i,p,k,m,lambda,epek,p*p+k*k-2*p*k*cos(phi));
	  exit(1);
	}
    }
  res/=2*RM;
  return res;
}

double Deltabar(double x0,double x, double y,double m)
{
  double res=0;
  double disc=4*x*x*y*y-(m*m-x0*x0+x*x+y*y)*(m*m-x0*x0+x*x+y*y);

  if (disc>0)
    {
      if (D==3)
	res=0.25/(x*y)*sign(x0);
    }
  return res;
}

double xyDeltabar(double x0,double x, double y,double m)
{
  double res=0;
  double m8=(m*m-x0*x0+x*x+y*y);
  double disc=4*x*x*y*y-m8*m8;

  if (disc>0)
    {
      if (D==3)
	res=0.25/(x*y)*sign(x0)*P2(m8/(2*x*y));
      else if (D==2)
	res=sign(x0)*f(m8/(2*x*y))/sqrt(disc);
    }
  // if (isnan(res))
  // {
  //   printf("disc=%f 2xy=%f f(cos)=%f\n",disc,2*x*y,f(m8/(2*x*y)));
  // }
  
  return res;
}

double G4minus(double q, double k,double m,double lambda)
{
  double res=0;
  double meps=10./RM;
  double q0=sqrt(q*q+m*m);
  double k0=sqrt(k*k+m*m);
  for (int i=1;i<RM;i++)
    for (int j=0;j<2*RM;j++)
      {
	double r=i*meps;
	double r0=(j-RM)*meps;
	res+=r*Deltabar(r0-q0,r,q,m)*Deltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
      }
  res*=meps*meps*8/(4*M_PI*M_PI);
  return res;
}

double xyG4minus(double q, double k,double m,double lambda)
{
  double res=0;
  double res2=0;
  double res3=0;
  double meps=E2/RM3;
  double q0=sqrt(q*q+m*m);
  double k0=sqrt(k*k+m*m);

  for (int i=1;i<RM3;i++)
    {
      double r=i*meps;
      double r0kp=sqrt(m*m+(r+k)*(r+k));
      double r0km=sqrt(m*m+(r-k)*(r-k));
      double r0qp=sqrt(m*m+(r+q)*(r+q));
      double r0qm=sqrt(m*m+(r-q)*(r-q));
      //positive part:
      double PPr0min=max(k0+r0km,q0+r0qm);
      double PPr0max=min(k0+r0kp,q0+r0qp);
      //negative part:
      double MMr0min=max(k0-r0kp,q0-r0qp);
      double MMr0max=min(k0-r0km,q0-r0qm);

      /*
      for (int j=0;j<2*RM3;j++)
	{
	  
	  double r0=(j-RM3)*meps;
	  
	  if (D==3)
	    res+=r*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	  else if (D==2)
	    {
	      double addme=r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	      res+=addme;
	      if ((r0>MMr0min)&&(r0<MMr0max))
		res3+=addme;
	      if ((r0>PPr0min)&&(r0<PPr0max))
		res3+=addme;
	      
	      
	      if (fabs(addme)>0)
		if ((r0<MMr0min)||(r0>MMr0max))
		{
		  if ((r0<PPr0min)||(r0>PPr0max))
		    {
		      printf("added for r0=%f my bounds (%f,%f)\n",r0,PPr0min,PPr0max);
		      printf("my NN bounds (%f,%f)\n",MMr0min,MMr0max);
		      printf("all roots km=%f kp=%f qm=%f qp=%f\n",r0km,r0kp,r0qm,r0qp);
		      printf("all other r=%f q=%f k=%f m=%f\n",r,q,k,m);
		    }
		}
		  
	    }
	    }*/

      //first part
      double eps1=(MMr0max-MMr0min)/RM3*0.01;
      for (int j=1;j<RM3;j++)
	{
	  double r0=MMr0min+j*eps1;
	  if (D==3)
	    res2+=eps1*r*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	  else if (D==2)
	    res2+=eps1*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);	  
	}
      //second part
      eps1=(MMr0max-MMr0min)/RM3*0.98;
      for (int j=0;j<RM3;j++)
	{
	  double r0=MMr0min+(MMr0max-MMr0min)*0.01+j*eps1;
	  if (D==3)
	    res2+=eps1*r*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	  else if (D==2)
	    res2+=eps1*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);	  
	}
      //third part
      eps1=(MMr0max-MMr0min)/RM3*0.01;
      for (int j=0;j<RM3;j++)
	{
	  double r0=MMr0min+(MMr0max-MMr0min)*0.99+j*eps1;
	  if (D==3)
	    res2+=eps1*r*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	  else if (D==2)
	    res2+=eps1*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);	  
	}


      
      eps1=(PPr0max-PPr0min)/RM3*0.01;
      for (int j=1;j<RM3;j++)
	{
	  double r0=PPr0min+j*eps1;
	  if (D==3)
	    res2+=eps1*r*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	  else if (D==2)
	    res2+=eps1*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	}
      eps1=(PPr0max-PPr0min)/RM3*0.98;
      for (int j=0;j<RM3;j++)
	{
	  double r0=PPr0min+(PPr0max-PPr0min)*0.01+j*eps1;
	  if (D==3)
	    res2+=eps1*r*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	  else if (D==2)
	    res2+=eps1*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	}
      eps1=(PPr0max-PPr0min)/RM3*0.01;
      for (int j=0;j<RM3;j++)
	{
	  double r0=PPr0min+(PPr0max-PPr0min)*0.99+j*eps1;
	  if (D==3)
	    res2+=eps1*r*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	  else if (D==2)
	    res2+=eps1*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	}
      
      
      
      
      }
  //printf("got res=%f vs res2=%f res3=%f\n",res*meps,res2,res3*meps);

  res2*=meps*8/(4*M_PI*M_PI);
  
  //res*=meps*meps*8/(4*M_PI*M_PI);

  
  
  return res2;
}

double xyG4plus(double q, double k,double m,double lambda)
{
  double res=0;
  double meps=E2/RM3;
  double q0=-sqrt(q*q+m*m);
  double k0=sqrt(k*k+m*m);


  
  for (int i=1;i<RM3;i++)
    {
      double r=i*meps;
      double r0kp=sqrt(m*m+(r+k)*(r+k));
      double r0km=sqrt(m*m+(r-k)*(r-k));
      double r0qp=sqrt(m*m+(r+q)*(r+q));
      double r0qm=sqrt(m*m+(r-q)*(r-q));
      //positive part:
      //double PPr0min=max(k0+r0km,q0+r0qm);
      //double PPr0max=min(k0+r0kp,q0+r0qp);
      //negative part:
      //double MMr0min=max(k0-r0kp,q0-r0qp);
      //double MMr0max=min(k0-r0km,q0-r0qm);
      //mixed part
      double MIr0min=max(k0-r0kp,q0+r0qm);
      double MIr0max=min(k0-r0km,q0+r0qp);

      /*
      for (int j=0;j<2*RM3;j++)
	{
	  
	  double r0=(j-RM3)*meps;
	  if (D==3)
	    res+=r*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	  else if (D==2)
	    {
	      double addme=r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	      res+=addme;
	      if (fabs(addme)>0)
		{
		  if ((r0<MIr0min)||(r0>MIr0max))
		    {
		      printf("added for r0=%f my bounds (%f,%f)\n",r0,MIr0min,MIr0max);
		      printf("my NN bounds (%f,%f)\n",MMr0min,MMr0max);
		      printf("all roots km=%f kp=%f qm=%f qp=%f\n",r0km,r0kp,r0qm,r0qp);
		      printf("all other r=%f q=%f k=%f m=%f\n",r,q,k,m);
		    }
		}
		  
	    }
	}
      res*=meps*meps*8/(4*M_PI*M_PI);

      */

      //first part
      double eps1=(MIr0max-MIr0min)/RM3*0.01;
      for (int j=1;j<RM3;j++)
	{
	  double r0=MIr0min+j*eps1;
	  if (D==3)
	    res+=eps1*r*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	  else if (D==2)
	    res+=eps1*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);	  
	}
      //second part
      eps1=(MIr0max-MIr0min)/RM3*0.98;
      for (int j=0;j<RM3;j++)
	{
	  double r0=MIr0min+(MIr0max-MIr0min)*0.01+j*eps1;
	  if (D==3)
	    res+=eps1*r*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	  else if (D==2)
	    res+=eps1*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);	  
	}
      //third part
      eps1=(MIr0max-MIr0min)/RM3*0.98;
      for (int j=0;j<RM3;j++)
	{
	  double r0=MIr0min+(MIr0max-MIr0min)*0.99+j*eps1;
	  if (D==3)
	    res+=eps1*r*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);
	  else if (D==2)
	    res+=eps1*r*xyDeltabar(r0-q0,r,q,m)*xyDeltabar(r0-k0,r,k,m)*(nB(r0-k0)-nB(r0-q0))*DR2(r0,r,m,lambda);	  
	}
      
    }
  res*=meps*8/(4*M_PI*M_PI);
      
  return res;
}



double intrhoDplus(double p,double k,double m,double lambda)
{
  //printf("D2 p=%f k=%f m=%f l=%f\n",p,k,m,lambda);
  double res=0;
  double meps=M_PI/RM;
  double epek=sqrt(p*p+m*m)+sqrt(k*k+m*m);
  double add=0;
  for (int i=0;i<2*RM;i++)
    {
      double phi=i*meps;
      if (D==2)
	add=rhoD(epek,sqrt(p*p+k*k-2*p*k*cos(phi)),m,lambda);
      else if (D==3)
	add=0.5*rhoD(epek,sqrt(p*p+k*k-2*p*k*cos(phi)),m,lambda)*fabs(sin(phi));
      //printf("phi=%f add=%f res=%f\n",phi,add,res);
      res+=add;
      if (isnan(res))
	{
	  printf("i=%i p=%f k=%f m=%f lambda=%f epek=%f sq=%f\n",i,p,k,m,lambda,epek,p*p+k*k-2*p*k*cos(phi));
	  exit(1);
	}
    }
  res/=2*RM;
  return res;
}



//angular average with extra angular weight from kx ky
double xyintrhoDminus(double p,double k,double m,double lambda)
{
  double res=0;
  double meps=M_PI/RM;
  double epek=sqrt(p*p+m*m)-sqrt(k*k+m*m);
  double add=0;
  for (int i=0;i<2*RM;i++)
    {
      double phi=i*meps;
      if (D==2)
	add=f(cos(phi))*rhoD(epek,sqrt(p*p+k*k-2*p*k*cos(phi)),m,lambda);
      else if (D==3)
	add=0.5*P2(cos(phi))*rhoD(epek,sqrt(p*p+k*k-2*p*k*cos(phi)),m,lambda)*fabs(sin(phi));
      //printf("phi=%f add=%f\n",phi,add);
      res+=add;      
    }
  res/=2*RM;
  return res;
}

//angular average with extra angular weight from kx ky
double xyintrhoDplus(double p,double k,double m,double lambda)
{
  //printf("D2 p=%f k=%f m=%f l=%f\n",p,k,m,lambda);
  double res=0;
  double meps=M_PI/RM;
  double epek=sqrt(p*p+m*m)+sqrt(k*k+m*m);
  double add=0;
  for (int i=0;i<2*RM;i++)
    {
      double phi=i*meps;
      if (D==2)
	add=f(cos(phi))*rhoD(epek,sqrt(p*p+k*k-2*p*k*cos(phi)),m,lambda);
      else if (D==3)
	add=0.5*P2(cos(phi))*rhoD(epek,sqrt(p*p+k*k-2*p*k*cos(phi)),m,lambda)*fabs(sin(phi));
      //printf("phi=%f add=%f res=%f\n",phi,add,res);
      res+=add;
      if (isnan(res))
	{
	  printf("i=%i p=%f k=%f m=%f lambda=%f epek=%f sq=%f\n",i,p,k,m,lambda,epek,p*p+k*k-2*p*k*cos(phi));
	  exit(1);
	}
    }
  res/=2*RM;
  return res;
}



double ImSigma(double p,double m,double lambda)
{
  //printf("\t INFO: called for Im Sigma at p=%.16g\n",p);
  double res=0;
  double ep=sqrt(p*p+m*m);
  double add=0;
  double add2;


  //#pragma omp parallel for private(add, add2) reduction(+: res)
  for (int i=0;i<RM2;i++)
    {
      double meps=10*p/RM2;
      double ek=m+i*meps;
      double k=sqrt(ek*ek-m*m);
      if (fabs(k-p)>safe)
	add=-(intrhoDminus(p,k,m,lambda)*(coth(0.5*ek)+coth(0.5*(ep-ek))));
      else
	{
	  double proxy=intrhoDminus(p,p+safe,m,lambda)/(ep-sqrt(m*m+(p+safe)*(p+safe)));
	  //printf("actual %f proxy %f\n",intrhoDminus(p,k,m,lambda)/(ep-ek),proxy);
	  add=-proxy*(ep-ek)*(coth(0.5*ek)+coth(0.5*(ep-ek)));
	}

      
      //second contribution      
      add2=-(intrhoDplus(p,k,m,lambda)*(coth(0.5*ek)-coth(0.5*(ep+ek))));
      //printf("k=%.16g add=%f add2=%f where %f\n",k,add,add2,intrhoDplus(p,k,m,lambda));
      add+=add2;


      if (D==3)
	add*=k;
      
      if (i==0)
	add*=0.5;
      if (i==RM2-1)
	add*=0.5;
      res+=add*meps;
      if (isnan(res))
	{
	  if (fabs(k-p)>safe)
	    add=-(intrhoDminus(p,k,m,lambda)*(coth(0.5*ek)+coth(0.5*(ep-ek))));
	  else
	    {
	      double proxy=intrhoDminus(p,p+safe,m,lambda)/(ep-sqrt(m*m+(p+safe)*(p+safe)));
	      //printf("actual %f proxy %f\n",intrhoDminus(p,k,m,lambda)/(ep-ek),proxy);
	      add=-proxy*(ep-ek)*(coth(0.5*ek)+coth(0.5*(ep-ek)));
	    }
	  
	  
	  //second contribution      
	  add2=-(intrhoDplus(p,k,m,lambda)*(coth(0.5*ek)-coth(0.5*(ep+ek))));
	  printf("i=%i add=%f add2=%f\n",i,add,add2);
	  printf("problem here %f\n",intrhoDminus(p,k,m,lambda));
	  exit(1);
	}
    }

  //printf("\t\treturning %f\n",res);
  
  return res/(2*M_PI);
}

double WeakImSigma(double p,double m,double lambda)
{
  //printf("\t INFO: called for Im Sigma at p=%.16g\n",p);
  double res=0;
  double ep=sqrt(p*p+m*m);
  double add=0;
  double add2;

  //m=0.01;

  //#pragma omp parallel for private(add, add2) reduction(+: res)
  for (int i=0;i<RM2;i++)
    {
      double meps=10*p/RM2;
      double ek=m+i*meps;
      double k=sqrt(ek*ek-m*m);
      if (fabs(k-p)<safe)
	{
	  double proxy=WintrhoDminus(p,p+safe,m,lambda)/(ep-sqrt(m*m+(p+safe)*(p+safe)));
	  //printf("actual %f proxy %f\n",intrhoDminus(p,k,m,lambda)/(ep-ek),proxy);
	  add=-proxy*(ep-ek)*(coth(0.5*ek)+coth(0.5*(ep-ek)));
	  
	}
      else 
	{
	  add=-(WintrhoDminus(p,k,m,lambda)*(coth(0.5*ek)+coth(0.5*(ep-ek))));
	  //add=-WrhoD(sqrt(p*p+m*m)-sqrt(k*k+m*m),p,m,lambda)*(coth(0.5*ek)+coth(0.5*(ep-ek)));
	}
      

      
      //second contribution      
      add2=-(WintrhoDplus(p,k,m,lambda)*(coth(0.5*ek)-coth(0.5*(ep+ek))));
      //add2=-WintrhoDplus(p,k,m,lambda)*(coth(0.5*ek)-1);
      //add2=-WrhoD(sqrt(p*p+m*m)+sqrt(k*k+m*m),p,m,lambda)*(coth(0.5*ek)-1);
      //printf("k=%.16g add=%f add2=%f where %f\n",k,add,add2,intrhoDplus(p,k,m,lambda));
      add+=add2;


      if (D==3)
	add*=k;
      
      if (i==0)
	add*=0.5;
      if (i==RM2-1)
	add*=0.5;
      res+=add*meps;
      
    }

  //printf("\t\treturning %f\n",res);
  
  return 4*res/(2*M_PI);
}


void VertexKernel(double **polys,double p,double m,double lambda,double* res)
{
  //printf("\t INFO: called for Im Pi at p=%.16g\n",p);

  double ep=sqrt(p*p+m*m);
  double add=0;
  double add2;



  // #pragma omp parallel for private(add, add2) reduction(+: res)
  for (int i=0;i<RM2;i++)
    {
      double meps=10*p/RM2;
      double eq=m+i*meps;
      double q=sqrt(eq*eq-m*m);
      if (fabs(q-p)>safe)
	add=-(intrhoDminus(p,q,m,lambda)*(coth(0.5*eq)+coth(0.5*(ep-eq))));
      else
	{
	  double proxy=intrhoDminus(p,p+safe,m,lambda)/(ep-sqrt(m*m+(p+safe)*(p+safe)));
	  //printf("actual %f proxy %f\n",intrhoDminus(p,q,m,lambda)/(ep-eq),proxy);
	  add=-proxy*(ep-eq)*(coth(0.5*eq)+coth(0.5*(ep-eq)));
	}

      
      //second contribution      
      add2=(intrhoDplus(p,q,m,lambda)*(coth(0.5*eq)-coth(0.5*(ep+eq))));
      //printf("q=%.16g add=%f add2=%f where %f\n",q,add,add2,intrhoDplus(p,q,m,lambda));
      add+=add2;

      if (D==3)
	add*=q;
      
      if (i==0)
	add*=0.5;
      if (i==RM2-1)
	add*=0.5;

     

      for (int i=0;i<quad;i++)
	res[i]+=add*meps*poly_eval(polys,i,eq)/(2*M_PI);
      
    }

  //printf("\t\treturning %f\n",res);
  
  //return res/(4*M_PI);
}

void xyVertexKernel(double **polys,double p,double m,double lambda,double* res)
{
  //printf("\t INFO: called for Im Pi at p=%.16g\n",p);

  double ep=sqrt(p*p+m*m);
  double add=0;
  double add2;



#pragma omp parallel for private(add, add2) 
  for (int i=0;i<RM2;i++)
    {
      double meps=10*p/RM2;
      double eq=m+i*meps;
      double q=sqrt(eq*eq-m*m);
      if (fabs(q-p)>safe)
	{
	  add=xyintrhoDminus(p,q,m,lambda);
	  add+=xyG4minus(p,q,m,lambda);
	}
      else
	{
	  double proxy=xyintrhoDminus(p,p+safe,m,lambda)/(ep-sqrt(m*m+(p+safe)*(p+safe)));
	  //printf("actual %f proxy %f\n",intrhoDminus(p,q,m,lambda)/(ep-eq),proxy);
	  add=proxy*(ep-eq);

	  proxy=xyG4minus(p,p+safe,m,lambda)/(ep-sqrt(m*m+(p+safe)*(p+safe)));
	  add+=proxy*(ep-eq);
	}

      
      //second contribution      
      add2=xyintrhoDplus(p,q,m,lambda);
      add2-=xyG4plus(p,q,m,lambda);

      add*=(coth(0.5*eq)+coth(0.5*(ep-eq)));
      add2*=coth(0.5*eq)-coth(0.5*(ep+eq));
      //printf("q=%.16g add=%f add2=%f where %f\n",q,add,add2,intrhoDplus(p,q,m,lambda));
      //add-=add2;

      if (D==3)
	{
	  add*=q;
	  add2*=q;
	}
      
      if (i==0)
	add*=0.5;
      if (i==RM2-1)
	add*=0.5;

      //for debugging !!! REMOVE!!!!
      //add*=3;
#pragma omp critical
      {
	for (int i=0;i<quad;i++)
	  res[i]-=(add-add2)*meps*poly_eval(polys,i,eq)/(2*M_PI);
      }
      
    }

  //printf("\t\treturning %f\n",res);
  
  //return res/(4*M_PI);
}



void printint(double p,double m,double lambda)
{
  double ep=sqrt(p*p+m*m);
  double meps=0.025;
  for (int i=0;i<30;i++)
    {
      double k=i*meps+0.001;
      double ek=sqrt(k*k+m*m);
      double myint=intrhoDminus(p,k,m,lambda)*coth(0.5*(ep-ek));
      printf("k=%f myint=%f\n",k,myint);
    }
}

void etaos(double *resmat, int verbose)
{

//prepare quadrature

  double stencils[quad];
  double weights[quad];
  double *polys[quad];
  for (int i=0;i<quad;i++)
    polys[i]=new double[quad];

  
  
  double lambda=resmat[0];

  double Rag[6]={1,50,0,1,0,double(lambda)};
  if (D==3)
    {
      Rag[2]=0.002;
      Rag[3]=1.9;
    }
  double m=convergeroot2(&massfunc,Rag);

  //debugging!!!! remove for real physics run!
  //m=1.18492;
  
  if (verbose)
    printf("\t INFO: mass parameter for lambda=%f is m=%f\n",lambda,m);
  resmat[1]=m;

  
  
  if (verbose)
    printf("\t INFO: generating %i-order quadrature for m=%f\n",quad,m);

  
  
  
  getquadrature(stencils,weights,quad,m,verbose);

  getpolys(stencils,weights,polys,quad,m,verbose);

  //  for (int i=0;i<quad;i++)
  //printf("Poly[%i]=%f+%f x+ %f x^2+%f x^3\n",i,polys[i][0],polys[i][1],polys[i][2],polys[i][3]);
  
  if (verbose)
    printf("\t\t DONE.\n");

  if (verbose)
    printf("\t INFO: generating %ix%i interpolation tables for Pi(m=%f)\n",NN,NN,m);
 
  interpolate_allPi(NN,m,dataR1,dataR2,dataRe);
  if (verbose)
    printf("\t\t DONE.\n");
  
 
  
  double isi[quad];
  double vp[quad][quad];
  
  for (int i=0;i<quad;i++)
    {
      double p=sqrt(stencils[i]*stencils[i]-m*m);
      isi[i]=ImSigma(p,m,lambda);
      
      double res[quad];
      for (int j=0;j<quad;j++)
	res[j]=0;
      xyVertexKernel(polys,p,m,lambda,res);
      for (int j=0;j<quad;j++)
	vp[i][j]=res[j];
      //printf("isi =%f vertex %f, %f %f\n",isi[i],res[0],res[1],res[2]);
    }
  
  VectorXd aj=VectorXd::Random(quad);
  for (int j=0;j<quad;j++)
    {
      aj(j)=0;
      for (int i=0;i<quad;i++)
	aj(j)+=weights[i]*poly_eval(polys,j,stencils[i])/isi[i];
      //printf("a[%i]=%f\n",j,aj(j));
    }

  VectorXd ajk2=VectorXd::Random(quad);
  for (int j=0;j<quad;j++)
    {
      ajk2(j)=0;
      for (int i=0;i<quad;i++)
	{
	  double k2=stencils[i]*stencils[i]-m*m;
	  ajk2(j)+=weights[i]*poly_eval(polys,j,stencils[i])/isi[i]*k2;
	}
      //printf("a[%i]=%f\n",j,aj(j));
    }
  
  MatrixXd A = MatrixXd::Random(quad,quad);

  for (int i=0;i<quad;i++)
    {
      double sfunc=0;
      for (int j=0;j<quad;j++)
	sfunc+=aj[j]*poly_eval(polys,j,stencils[i]);
      
      //printf("x=%f isi=%f sfunc=%f\n",stencils[i],isi[i],1/sfunc);
    }
  
  for (int a=0;a<quad;a++)
    for (int b=0;b<quad;b++)
      {
	A(a,b)=0;
	if (a==b)
	  A(a,b)=1;
	for (int i=0;i<quad;i++)
	  {
	    double sfunc=0;
	    for (int j=0;j<quad;j++)
	      sfunc+=aj[j]*poly_eval(polys,j,stencils[i]);
	    A(a,b)-=weights[i]*vp[i][a]*sfunc*poly_eval(polys,b,stencils[i]);
	  }
	//printf("%i %i vertex=%f\n",a,b,Vertex[a][b]);
      }
  
  VectorXd out =  A.colPivHouseholderQr().solve(ajk2);

  //printf("final coefficients\n");
  //cout << out;
  
  //  for (int a=0;a<quad;a++)
  //  printf("1 loop %f vertex %f poly %f %f\n",ajk2(a),out(a),polys[a][0],polys[a][1]);

  double co[6]={0,0,0,0,0,0};
 
  //for (int a=0;a<quad;a++)
  //  for (int b=0;b<quad;b++)
  //	co[b]+=out(a)*polys[a][b];
  
  //printf("{%f,%f,%f,%f}\n",co[0],co[1],co[2],co[3]);
  
  double num1=0;
  double num2=0;
  double num3=0;
  double den=0;
  for (int i=0;i<quad;i++)
    {
      double k2=stencils[i]*stencils[i]-m*m;
      for (int j=0;j<quad;j++)
	{	      
	  num1+=weights[i]*k2*k2*aj(j)*poly_eval(polys,j,stencils[i]);
	  num3+=weights[i]*k2*ajk2(j)*poly_eval(polys,j,stencils[i]);
	  num2+=weights[i]*k2*out(j)*poly_eval(polys,j,stencils[i]);	      
	}
      den+=weights[i]*stencils[i]*k2;
    }


  
  double s=den/(4*M_PI);
  if (D==3)
    s=den/(6*M_PI*M_PI);
  double eos1=num1/(4*den);
  double eos2=num2/(4*den);
  double eos3=num3/(4*den);

  if (D==3)
    {
      eos1*=0.8;
      eos2*=0.8;
      eos3*=0.8;
    }

  resmat[2]=s;
  resmat[3]=eos1;
  resmat[4]=eos2;

  //  if (!singlerun)
  //printf("%f\t%f\t%f\t%f=?%f\t%f\n",lambda,m,s,eos1,eos3,eos2);


  
  for (int i=0;i<quad;i++)
    delete [] polys[i];
  
  
}


void weaketaos(double *resmat, int verbose)
{

  //assumes all the calculations have already been done
double stencils[quad];
  double weights[quad];
  double *polys[quad];
  for (int i=0;i<quad;i++)
    polys[i]=new double[quad];

  double lambda=resmat[0];
  double m=resmat[1];
  printf("\t INFO: mass parameter for lambda=%f is m=%f\n",lambda,m);
  getquadrature(stencils,weights,quad,m,verbose);
  getpolys(stencils,weights,polys,quad,m,verbose);

  
  
  double isi[quad];
  double vp[quad][quad];
  
  for (int i=0;i<quad;i++)
    {
      double p=sqrt(stencils[i]*stencils[i]-m*m);
      isi[i]=WeakImSigma(p,m,lambda);
      
    }
  
  VectorXd aj=VectorXd::Random(quad);
  for (int j=0;j<quad;j++)
    {
      aj(j)=0;
      for (int i=0;i<quad;i++)
	aj(j)+=weights[i]*poly_eval(polys,j,stencils[i])/isi[i];
      //printf("a[%i]=%f\n",j,aj(j));
    }

  

  for (int i=0;i<quad;i++)
    {
      double sfunc=0;
      for (int j=0;j<quad;j++)
	sfunc+=aj[j]*poly_eval(polys,j,stencils[i]);
      
      //printf("x=%f isi=%f sfunc=%f\n",stencils[i],isi[i],1/sfunc);
    }
  
  
  
  //printf("{%f,%f,%f,%f}\n",co[0],co[1],co[2],co[3]);
  
  double num1=0;
 
  double den=0;
  for (int i=0;i<quad;i++)
    {
      double k2=stencils[i]*stencils[i]-m*m;
      for (int j=0;j<quad;j++)
	{	      
	  num1+=weights[i]*k2*k2*aj(j)*poly_eval(polys,j,stencils[i]);	     
	}
      den+=weights[i]*stencils[i]*k2;
    }


  
  double s=den/(4*M_PI);
  if (D==3)
    s=den/(6*M_PI*M_PI);
  double eos1=num1/(4*den);

  if (D==3)
    {
      eos1*=0.8;
    
    }

  printf("eos1=%f\n",eos1);


  
  for (int i=0;i<quad;i++)
    delete [] polys[i];
  
  
}


void read_params(char *paramsfile)
{
  
  ifstream parameters;
  //parameters.open(paramsfile, ios::in);
  parameters.open(paramsfile);

  
  if(parameters.is_open())
    {
      string line;
      for (int i=0;i<8;i++)
	getline(parameters,line);
      
      int dummyint;
      char dummychar[255];
      double dummydouble;

      parameters >> dummychar;      
      parameters >> D;
      //D=dummyint;
      //printf("got %s %i D=%i\n",dummychar,dummyint,D);
      for (int i=0;i<4;i++)
	parameters >> dummychar;
      
      parameters >> dummychar;      
      parameters >> NN;

      for (int i=0;i<3;i++)
	parameters >> dummychar;
      
      parameters >> dummychar;      
      parameters >> quad;

      
      //NN=dummyint;
      //printf("%s %i NN= %i\n",dummychar,dummyint,NN);
      //getline(parameters,line);

      for (int i=0;i<6;i++)
	parameters >> dummychar;

      
      
      parameters >> dummychar;

      //printf("we have %s %i and need %i\n",dummychar,dummychar[0],(int)'Y');
      
      if (dummychar[0]=='Y')
	singlerun=1;
      else
	singlerun=0;

      //printf("We have d=%s\n",dummychar);
      
      parameters >> lambda;

      for (int i=0;i<4;i++)
	parameters >> dummychar;

      parameters >> vlevel;

      for (int i=0;i<3;i++)
	parameters >> dummychar;

      parameters >> RM;

      for (int i=0;i<3;i++)
	parameters >> dummychar;

      parameters >> RM2;

      for (int i=0;i<3;i++)
	parameters >> dummychar;

      parameters >> RM3;
      
      if (singlerun)
	printf("INFO:\t Doing a single run with lambda=%f\n",lambda);
      else
	printf("INFO:\t Doing multiple runs as a function of lambda\n");
      printf("INFO:\t found D=%i NN=%i quad=%i verbose=%i RM=%i RM2=%i RM3=%i\n",D,NN,quad,vlevel,RM,RM2,RM3);
      

      
    }
  else
    {
      printf("\nERROR: params file %s could not be opened\n",paramsfile);
      ABORT=1;
    }
  parameters.close();
  

  

}


int main(int argc, char *argv[])
{
  
  
  char infile[255];
  sprintf(infile,"params.txt");

  
  printf("INFO:\t This is the shear viscosity solver for O(N) model, v1.0\n");
  printf("INFO:\t Using parameter file %s for input parameters\n",infile);

  read_params(infile);

  
  m3func=new double*[SF3];
  for (int i=0;i<SF3;i++)
    m3func[i]=new double[2];
  
  
  
  //prepare interpolation workspace
  dataR1=new double*[NN];
  dataR2=new double*[NN];
  dataRe=new double*[NN];
  for (int i=0;i<NN;i++)
    {
      dataR1[i]=new double[NN];
      dataR2[i]=new double[NN];
      dataRe[i]=new double[NN];
    }


  if (D==3)
    {
      //read special functions table for interpolation
      ifstream mfile;
      char mname[255];
      sprintf(mname,"data-tables/D3sumk1.dat");
      
      mfile.open(mname);
      if (mfile.is_open())
	{
	  for (int i=0;i<SF3;i++)
	    {
	      mfile >> m3func[i][0];
	      mfile >> m3func[i][1];
	    }
	}
      else
	{
	  printf("ERROR: could not read file %s\n",mname);
	}
    

      //printf("a check %f \n",eval_lin_int(m3func,0.5,SF3));
    }
  
  

  //old_stuff();
  
  
  FILE *outfile;
  
  double resmat[5]={0,0,0,0,0};

  if (!singlerun)
    {
      printf("INFO: generating table of eta/s values as a function of lambda\n");
      string name="results/etaos-D"+to_string(D)+"-Q"+to_string(quad)+"-NN"+to_string(NN)+"-RM"+to_string(RM)+"-RM2-"+to_string(RM2)+"RM3-"+to_string(RM3)+".dat";
      outfile = fopen(name.c_str(),"w");
      fprintf(outfile,"#lambda/T\t mass/T\t\t s/T^3 \t\t eta/s(R2)\t eta/s (full)\n");
      printf("#lambda/T\t mass/T\t\t s/T^3 \t\t eta/s(R2)\t eta/s (full)\n");

      int lowr=-12;
      int highr=34;
      double step=0.5;
      
      if (D==3)
	{
	  step=1;
	  lowr=-6;
	  highr=4;

	  //step=0.5;
	  //lowr=1;
	  //highr=17;
	  
	}
      
      for (int run=lowr;run<highr;run++)
	{
	  lambda=pow(2,run*step);
	  //lambda=run*0.5/24;
	  //printf("lambda=%f\n",lambda);
	  resmat[0]=lambda;
	  etaos(resmat,vlevel);
	  
	  for (int u=0;u<5;u++)
	    {
	      fprintf(outfile,"%f\t",resmat[u]);
	      printf("%f\t",resmat[u]);
	    }
	  fprintf(outfile,"\n");
	  printf("\n");
	  fflush(outfile);
	  
	}
      fclose(outfile);    
    }
  else
    {
      
      printf("INFO: generating eta/s for a single value of lambda=%f\n",lambda);
      /*
      double p0=4;
      double p=2;
      double m=0.720076;
      printf("basic check %f %f\n",ImPi0(p0,p,m),ImPiT(p0,p,m));
      printf("basic check %f %f\n",RePi0(p0,p,m),RePiT(p0,p,m));
      exit(1);
      */

      resmat[0]=lambda;
      
      etaos(resmat,vlevel);

      //weaketaos(resmat,vlevel);
      
      //OK, let's see

      
      
      
      /*
      double p=2;
      double k=2.5;
      double m=resmat[1];
      printf("rhominus %f rhoplus %f\n",intrhoDminus(p,k,m,lambda),intrhoDplus(p,k,m,lambda));
     

      printf("xyrhominus %f xyrhoplus %f\n",xyintrhoDminus(p,k,m,lambda),xyintrhoDplus(p,k,m,lambda));
      printf("xyG4minus %f xyG4plus %f\n",xyG4minus(p,k,m,lambda),xyG4plus(p,k,m,lambda));
      printf("check %f xyG4plus %f\n",xyG4minus(k,p,m,lambda),xyG4plus(k,p,m,lambda));

      
      printf("samples im=%f %f re=%f %f %f\n",ImPi0(p0,p,m),ImInt(p0,p,m,dataR1,dataR2,NN),RePi0(p0,p,m),Reint(p0,p,dataRe,NN),0.5/lambda);
      */
      /*
      fstream outf;
      outf.open("test.dat",ios::out);
      if (outf.is_open())
	{
	  for (int i=0;i<400;i++)
	    {
	      double w=i*0.05;
	      double lambdam=1/(1/lambda-log(resmat[1])*0.5/M_PI/M_PI);
	      outf << w << "\t";
	      outf << -8*rhoD(w,2,resmat[1],lambda)/(24*lambdam)/(24*lambdam) << "\t";
	      outf << lambdam << endl;
	    }
	  
      
	  outf.close();
	}
      else
	printf("ERROR couldn open outfile\n");

      */
      //printf("result for thermal width: isi(p=5)=%f\n",ImSigma(5,resmat[1],lambda));
      
      //printf("check rhoDplus %f\n",intrhoDplus(5,4,resmat[1],lambda));
      //printf("check rhoDminus %f\n",intrhoDminus(5,4,resmat[1],lambda));
      /*
      double res=0;
      double res2=0;
      for (int i=0;i<2*RM;i++)
	{
	  double meps=M_PI/RM;
	  double p=5;
	  double k=4;
	  double m=resmat[1];
	  double phi=i*meps;
	  double y=rhoD(sqrt(p*p+m*m)+sqrt(k*k+m*m),sqrt(p*p+k*k-2*p*k*cos(phi)),m,lambda);
	  res+=y*fabs(sin(phi))/(4*RM);
	  res2+=sin(phi)*sin(phi)/(2*RM);
	  printf("phi=%f y=%f res=%f res2=%f\n",phi,y,res,res2);
	  }

      double res=0;
      double p=5;
      double m=resmat[1];
      double ep=sqrt(p*p+m*m);
      double add=0;
      double add2;
      double safe=0.025;
      for (int i=0;i<RM2;i++)
	{
	  double meps=10*p/RM2;
	  double ek=m+i*meps;
	  double k=sqrt(ek*ek-m*m);
	  if (fabs(k-p)>safe)
	    add=-2*(intrhoDminus(p,k,m,lambda)*(coth(0.5*ek)+coth(0.5*(ep-ek))));
	  else
	    {
	      double proxy=intrhoDminus(p,p+safe,m,lambda)/(ep-sqrt(m*m+(p+safe)*(p+safe)));
	      //printf("actual %f proxy %f\n",intrhoDminus(p,k,m,lambda)/(ep-ek),proxy);
	      add=-2*proxy*(ep-ek)*(coth(0.5*ek)+coth(0.5*(ep-ek)));
	    }
	  //second contribution      
	  add2=-2*(intrhoDplus(p,k,m,lambda)*(coth(0.5*ek)-coth(0.5*(ep+ek))));
	  printf("k=%.16g add=%f add2=%f where %f\n",k,add,add2,intrhoDplus(p,k,m,lambda));
	  add+=add2;
	  res+=add*k*meps;
	}
	printf("get %f\n",res/(2*M_PI));*/
      
      printf("Result for lambda=%f\t m=%f\t s=%f eos1=%f eos2=%f\n",lambda,resmat[1],resmat[2],resmat[3],resmat[4]);
      //printf("Im Sigma/l^2 %f vs %f\n",ImSigma(5,resmat[1],lambda)/(lambda*lambda),WeakImSigma(5,resmat[1],lambda));
    }
  
  
  

  for (int i=0;i<NN;i++)
    {
      delete [] dataR1[i];
      delete [] dataR2[i];
      delete [] dataRe[i];
    }
  for (int i=0;i<SF3;i++)
    delete [] m3func[i];
}
