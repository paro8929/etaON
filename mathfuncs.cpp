int sgn(double x)
{
  int out=0;
  if (x>0) out=1;
  else if (x<0) out=-1;

  return out;
}

double min(double a,double b)
{
  double res=0;
  if (a<b)
    res=a;
  else
    res=b;
  return res;
}

double coth(double x)
{
  double res=cosh(x)/sinh(x);
  return res;
}

double nB(double x)
{
  double res=1/(exp(x)-1);
  return res;
}

double nBprime(double x)
{
  double res=1/(2-exp(x)-exp(-x));
  return res;
}

double polylog(int n,double x)
{
  double res=0;
  for (int i=1;i<100;i++)
    {
      double c=pow(x,i)/pow(i,n);
      //if (isnan(c)&&(c>1.e-10))
      res+=c;
    //printf("i=%i res=%f x=%f\n",i,res,x);
    }
  return res;
}

double eval_lin_int(double **data, double x,int size)
{
  double res=0;
  int i=0;
  if ((x>=data[0][0])&&(x<data[size-1][0]))
    {
      while (data[i][0]<x)
	i++;
      i--;
      res=linear_int(data[i][0],data[i+1][0],data[i][1],data[i+1][1],x);
      //printf("here with %f found %i %f\n",x,i,res);
    }
  else
    printf("WARNING: %f out of bounds [%f,%f]\n",x,data[0][0],data[size-1][0]);
  return res;
}


//int dE nB^\prime(E)*sqrt(E^2-m^2) E^n
double D3ints(int n, double m,long steps,double eps)
{
  double res=0;
  for (int i=1;i<steps;i++)
    {
      double y=m+i*eps;
      res-=nBprime(y)*sqrt(y*y-m*m)*pow(y,n);
    }
  return res*eps;
}

//adapt integrator steps to get to answer
double D3ints(int n, double m)
{

  double eps=EPS*10;
  long steps=MAX/10;
  int cont=1;
  double re=1;
  double res1=D3ints(n,m,steps,eps);

  //relative error aim
  double AIM=1.e-6;

  
  while(cont)
    {

      steps*=2;
      double res2=D3ints(n,m,steps,eps);
      re=fabs(res1-res2)/res2;
      //printf("s=%li iteration r1=%f r2=%f re=%f %i\n",steps,res1,res2,re,(re>1.e-2));
      
      res1=res2;
      if ((re<AIM)||(steps>1.e6))
	cont=0;
      
    }

  if (re>AIM)
    printf("\t D3 integrator warning: could not reach relative error limit for %li steps; err=%.2g>%.2g=AIM\n",steps,re,AIM);
  
  cont=1;

   while(cont)
    {

      steps*=2;
      eps/=2.;
      double res2=D3ints(n,m,steps,eps);
      re=fabs(res1-res2)/res2;
      //printf("s=%li iteration r1=%f r2=%f re=%f %i\n",steps,res1,res2,re,(re>1.e-2));
      
      res1=res2;
      if ((re<AIM)||(steps>1.e6))
	cont=0;
      
    }

   if (re>AIM)
    printf("\t D3 integrator warning: could not reach relative error limit for %li steps; %.2g>%.2g\n",steps,re,AIM);
   
  
  return res1;
}


double P2(double x)
{
  double res=3*x*x-1;
  return 0.5*res;
}

double f(double x)
{
  double res=2*x*x-1;
  return res;
}
