double RePi0(double p0, double p, double m)
{
  double res=0;
  double pm2=p0*p0-p*p;
  if (D==2)
    {
      if (pm2<0)
	res=atan(sqrt(-pm2)/(2*m))/(2*M_PI*sqrt(-pm2));
      else if (pm2>0)
	res=log(fabs((2*m+sqrt(pm2))/(2*m-sqrt(pm2))))/(4*M_PI*sqrt(pm2));
      else
	res=0.25/(m*M_PI);
    }
  else if (D==3)
    {
      
      double arg=1-4*m*m/pm2;
      res=1;
      if ((pm2>0)&&(arg<0))
	{
	  double b=sqrt(-arg);
	  res-=b*atan(1./b);
	  //printf("h1, res=%f\n",res);
	}
      if (pm2<0)
	{
	  double b=sqrt(arg);
	  res+=0.5*b*log(fabs((1-b)/(1+b)));
	  //printf("h2, res=%f\n",res);
	}
      if (pm2-4*m*m>0)
	{
	  double b=sqrt(arg);
	  res+=0.5*b*log(fabs((1-b)/(1+b)));
	  //printf("h3, res=%f\n",res);
	}
      
      res/=4*M_PI*M_PI;
	  
    }
  else
    {
      printf("D=%i not coded up\n",D);
      exit(1);
    }
  return res;
}

double ImPi0(double p0, double p, double m)
{
  double res=0;
  double pm2=p0*p0-p*p;
  if (D==2)
    {
      if (pm2>4*m*m)
	{
	  res=sgn(p0)/(4*sqrt(pm2));
	}
    }
  else if (D==3)
    {
      if (pm2>4*m*m)
	{
	  res=sgn(p0)*sqrt(1-4*m*m/pm2);
	  res/=8*M_PI;
	}
    }
  else
    {
      printf("D=%i not coded up\n",D);
      exit(1);
    }
  return res;
}



double RePiT(double p0,double p,double m)
{
    
  double res=0;
  double pm2=p0*p0-p*p;

  //printf("here with %f %f\n",p0,p);
  if (p0<0)
    {
      p0*=-1;
    }
  if (D==2)
    {
      if (pm2>0)
	{
	  if (p>0)
	    {
	      
	      double lowc=m;
	      double ek;
	      double c1,c2;
	      ek=lowc;
	      c1=0.5*EPS*nB(ek)/sqrt((pm2+2*p0*ek)*(pm2+2*p0*ek)-4*p*p*(ek*ek-m*m));
	      for (int i=1;i<MAX;i++)
		{
		  ek=i*EPS+lowc;	      
		  c1+=EPS*nB(ek)/sqrt((pm2+2*p0*ek)*(pm2+2*p0*ek)-4*p*p*(ek*ek-m*m));  
		}
	      res+=c1;
	      //printf("c1=%f\n",c1);
	      
	      //printf("c1=%F res=%f\n",c1,res);
	      if (pm2<4*m*m)
		{
		  ek=lowc;
		  c2=0.5*EPS*sgn(pm2-2*p0*ek)*nB(ek)/sqrt((pm2-2*p0*ek)*(pm2-2*p0*ek)-4*p*p*(ek*ek-m*m));
		  for (int i=1;i<MAX;i++)
		    {
		      ek=i*EPS+lowc;	      
		      c2+=EPS*sgn(pm2-2*p0*ek)*nB(ek)/sqrt((pm2-2*p0*ek)*(pm2-2*p0*ek)-4*p*p*(ek*ek-m*m));  
		    }
		  res+=c2;
		  
		}
	      else
		{
		  //printf("\t why are we here %f\n",pm2);
		  //printf("so far%f\n",res);
		  lowc=0.5*p0-0.5*p*sqrt(1-4*m*m/pm2);
		  double highc=0.5*p0+0.5*p*sqrt(1-4*m*m/pm2);
		  double meps=(lowc-m)/(MAX-1);
		  if (meps>0)
		    {
		      ek=m;
		      c2=nB(ek)/sqrt(pm2)*asinh(sqrt(lowc-ek)/sqrt(highc-lowc));
		      c2+=0.5*meps*nBprime(ek)/sqrt(pm2)*asinh(sqrt(lowc-ek)/sqrt(highc-lowc));
		      //printf("boundary %f where %.16g %f ss=%f\n",c2,lowc,highc,lowc-ek);
		      //printf("why %.16g %.16g\n",0.5*p0,0.5*p*sqrt(1-4*m*m/pm2));
		      //printf("p0=%f p=%f meps=%f\n",p0,p,meps);
		      
		      for (int i=1;i<MAX-1;i++)
			{
			  ek=i*meps+m;	      
			  c2+=meps*nBprime(ek)/sqrt(pm2)*asinh(sqrt(lowc-ek)/sqrt(highc-lowc));
			  //printf("adding %f\n",meps*nBprime(ek)/sqrt(pm2)*asinh(sqrt(lowc-ek)/sqrt(highc-lowc)));
			}		  
		      c2*=sgn(pm2-2*p0*ek);
		      res+=c2;
		      
		      
		      
		      //printf("c2=%f res=%f\n",c2,res);
		    }
		  //printf("c2= %f meps=%f lowc=%f, highc=%f\n",c2,meps,lowc,highc);
		  
		  double c3=0;
		  ek=highc;
		  c3=-0.5*nB(ek)/sqrt(pm2)*log(pm2*(2*ek-p0));
		  //printf("boundary c3  %f\n",c3);
		  
		  c3-=0.5*0.5*EPS*nBprime(ek)/sqrt(pm2)*log(pm2*(2*ek-p0));
		  for (int i=1;i<MAX;i++)
		    {
		      ek=i*EPS+highc;	      
		      c3-=EPS*0.5*nBprime(ek)/sqrt(pm2)*log(pm2*(2*ek-p0)+sqrt(pm2)*sqrt(4*m*m*p*p+pm2*(4*ek*ek-4*ek*p0+pm2)));
		      //printf("ek=%f hu %f\n",ek,0.5*nBprime(ek)/sqrt(pm2)*log(pm2*(2*ek-p0)+sqrt(pm2)*sqrt(4*m*m*p*p+pm2*(4*ek*ek-4*ek*p0+pm2))));
		    }
		  c3*=sgn(pm2-2*p0*ek);
		  res+=c3;
		  //printf("c3=%f res=%f\n",c3,res);
		  //printf("c3= %f\n",c3);
		  //printf("myres=%f\n",res);
		}
	    }
	  else //p=0
	    {
	      double lowc=m;
	      double ek;
	      ek=lowc;

	      if (p0<2*m)
		{
		  for (int i=0;i<MAX;i++)
		    {
		      ek=i*EPS+lowc;
		      double r=p0*p0-4*ek*ek;
		      res+=2*EPS*nB(ek)*r/(r*r+4*p0*p0*1.e-10);
		    }
		}
	      else
		{
		  double meps=(0.5*p0-lowc)/MAX;
		  for (int i=0;i<MAX;i++)
		    {
		      ek=lowc+i*meps;
		      double r=p0*p0-4*ek*ek;
		      res+=2*meps*nB(ek)*r/(r*r+4*p0*p0*1.e-10);		      
		    }
		  for (int i=0;i<MAX;i++)
		    {
		      ek=0.5*p0+i*meps;
		      double r=p0*p0-4*ek*ek;
		      res+=2*meps*nB(ek)*r/(r*r+4*p0*p0*1.e-10);		      
		    }
		  
		  for (int i=0;i<MAX;i++)
		    {
		      ek=i*EPS+p0-m;
		      double r=p0*p0-4*ek*ek;
		      res+=2*EPS*nB(ek)*r/(r*r+4*p0*p0*1.e-10);
		    }
		  
		}
	    }
	}
      else if (pm2<0)
	{
	  double lowc=-p0*0.5-0.5*p*sqrt(1-4*m*m/pm2);
	  double highc=-p0*0.5+0.5*p*sqrt(1-4*m*m/pm2);
	  double meps=(highc-m)/(MAX-1);
	  double ek;
	  double c1,c2;

	  //printf("lowc=%f highc=%f\n",lowc,highc);
	  if (meps>0)
	    {
	      ek=m;
	      c1=-nB(ek)*asin(sqrt((ek-lowc)/(highc-lowc)));
	      c1-=meps*0.5*nBprime(ek)*asin(sqrt((ek-lowc)/(highc-lowc)));
	      for (int i=1;i<MAX-1;i++)
		{
		  ek=i*meps+m;
		  c1-=meps*nBprime(ek)*asin(sqrt((ek-lowc)/(highc-lowc))); 
		}
	      ek=highc;
	      c1-=0.5*meps*nBprime(ek)*M_PI*0.5;
	      c1+=nB(ek)*M_PI*0.5;
	      c1*=sgn(pm2+2*p0*m)/sqrt(-pm2);
	      //printf("c1=%f\n",c1);
	      res+=c1;
	    }

	  ek=m;
	  lowc=p0*0.5-0.5*p*sqrt(1-4*m*m/pm2);
	  highc=p0*0.5+0.5*p*sqrt(1-4*m*m/pm2);
	  meps=(highc-m)/(MAX-1);
	  c2=-nB(ek)*asin(sqrt((ek-lowc)/(highc-lowc)));
	  c2-=meps*0.5*nBprime(ek)*asin(sqrt((ek-lowc)/(highc-lowc)));
	  for (int i=1;i<MAX-1;i++)
	    {
	      ek=i*meps+m;
	      c2-=meps*nBprime(ek)*asin(sqrt((ek-lowc)/(highc-lowc))); 
	    }
	  ek=highc;
	  c2-=0.5*meps*nBprime(ek)*M_PI*0.5;
	  c2+=nB(ek)*M_PI*0.5;
	  c2*=sgn(pm2-2*p0*m)/sqrt(-pm2);
	  //printf("c2=%f\n",c2);
	  res+=c2;
	  //if (isnan(res))
	  // printf("trouble here %f  c1=%f c2=%f\n",res,c1,c2);
	}
    }
  else if (D==3)
    {
      
      if (pm2<0)
	{
	  double b2=1-4*m*m/pm2;
	  double low=m;
	  double k=0.5*(p0*sqrt(b2)-p);
	  double ek=low;
	  double high=sqrt(k*k+m*m);
	  double meps=(high-low)/MAX;
	  k=sqrt(ek*ek-m*m);
	  res=meps*0.5*nB(ek)*log(fabs(((k+0.5*p)*(k+0.5*p)-0.25*p0*p0*b2)/((k-0.5*p)*(k-0.5*p)-0.25*p0*p0*b2)));
	  

	  
	  for (int i=1;i<MAX;i++)
	    {
	      ek=low+i*meps;
	      k=sqrt(ek*ek-m*m);
	      res+=meps*nB(ek)*log(fabs(((k+0.5*p)*(k+0.5*p)-0.25*p0*p0*b2)/((k-0.5*p)*(k-0.5*p)-0.25*p0*p0*b2)));
	      
	    }
	  
	  


	  if (p0>0)
	    {
	      low=high;
	      k=0.5*(p0*sqrt(b2)+p);
	      high=sqrt(k*k+m*m);
	      meps=(high-low)/MAX;
	      
	      for (int i=1;i<MAX;i++)
		{
		  ek=low+i*meps;
		  k=sqrt(ek*ek-m*m);
		  res+=meps*nB(ek)*log(fabs(((k+0.5*p)*(k+0.5*p)-0.25*p0*p0*b2)/((k-0.5*p)*(k-0.5*p)-0.25*p0*p0*b2)));
		}
	    }

	  low=high;
	  for (int i=1;i<MAX;i++)
	    {
	      ek=low+i*EPS;
	      double k=sqrt(ek*ek-m*m);
	      res+=EPS*nB(ek)*log(fabs(((k+0.5*p)*(k+0.5*p)-0.25*p0*p0*b2)/((k-0.5*p)*(k-0.5*p)-0.25*p0*p0*b2)));
	    }
	  
	  res/=-4*M_PI*p;
	  

	  
	}
      else if (pm2<4*m*m)
	{
	  if (pm2>0)
	    {
	      double b2=1-4*m*m/pm2;
	      double ek=m;
	      double k=sqrt(ek*ek-m*m);
	      if (p>0)
		{
		  res=0.5*nB(ek)*log(fabs(((k+0.5*p)*(k+0.5*p)-0.25*p0*p0*b2)/((k-0.5*p)*(k-0.5*p)-0.25*p0*p0*b2)));
		  for (int i=1;i<MAX;i++)
		    {
		      ek=m+i*EPS;
		      k=sqrt(ek*ek-m*m);
		      res+=nB(ek)*log(fabs(((k+0.5*p)*(k+0.5*p)-0.25*p0*p0*b2)/((k-0.5*p)*(k-0.5*p)-0.25*p0*p0*b2)));
		    }
		  
		  res*=-EPS;
		  res/=4*M_PI*p;
		}
	      else
		{
		  res=0.5*nB(ek)*k/(ek*ek-p0*p0*0.25);
		  for (int i=1;i<MAX;i++)
		    {
		      ek=m+i*EPS;
		      k=sqrt(ek*ek-m*m);
		      res+=nB(ek)*k/(ek*ek-p0*p0*0.25);
		    }
		  
		  res*=-EPS;
		  res/=2*M_PI;
		}
	      //else it's just zero
	    }
	}

      if (pm2>4*m*m)
	{
	  //2 poles
	  if (p>0)
	    {
	      double b2=1-4*m*m/pm2;
	      double low=m;
	      double k=0.5*(p-p0*sqrt(b2));
	      double high=sqrt(k*k+m*m);
	      double ek=low;
	      k=sqrt(ek*ek-m*m);
	      double meps=(high-low)/MAX;
	      res=meps*0.5*nB(ek)*log(fabs(((k+0.5*p)*(k+0.5*p)-0.25*p0*p0*b2)/((k-0.5*p)*(k-0.5*p)-0.25*p0*p0*b2)));
	      for (int i=1;i<MAX;i++)
		{
		  ek=low+i*meps;
		  k=sqrt(ek*ek-m*m);
		  res+=meps*nB(ek)*log(fabs(((k+0.5*p)*(k+0.5*p)-0.25*p0*p0*b2)/((k-0.5*p)*(k-0.5*p)-0.25*p0*p0*b2)));
		}

	      if (p0>0)
		{
		  low=high;
		  k=0.5*(p0*sqrt(b2)+p);
		  high=sqrt(k*k+m*m);
		  meps=(high-low)/MAX;
		  
		  for (int i=1;i<MAX;i++)
		    {
		      ek=low+i*meps;
		      double k=sqrt(ek*ek-m*m);
		      res+=meps*nB(ek)*log(fabs(((k+0.5*p)*(k+0.5*p)-0.25*p0*p0*b2)/((k-0.5*p)*(k-0.5*p)-0.25*p0*p0*b2)));
		    }
		}

	     
	      
	      
	      low=high;
	      for (int i=1;i<MAX;i++)
		{
		  ek=low+i*EPS;
		  double k=sqrt(ek*ek-m*m);
		  res+=EPS*nB(ek)*log(fabs(((k+0.5*p)*(k+0.5*p)-0.25*p0*p0*b2)/((k-0.5*p)*(k-0.5*p)-0.25*p0*p0*b2)));
		  
		}
	       	       
	      res/=-4*M_PI*p;
	    }
	  else //p=0
	    {
	      double low=m;
	      double high=p0/2;
	      double ek=low;
	      double k=sqrt(ek*ek-m*m);
	      double meps=(high-low)/MAX;

	      res=-0.5*meps*nB(ek)*k/(ek*ek-p0*p0*0.25);
	      for (int i=1;i<MAX;i++)
		{
		  ek=m+i*meps;
		  k=sqrt(ek*ek-m*m);
		  res-=meps*nB(ek)*k/(ek*ek-p0*p0*0.25);
		}
	      
	     
	      low=high;
	      ek=low;
	      k=sqrt(ek*ek-m*m);
	      res-=0.5*EPS*nB(ek)*k/(ek*ek-p0*p0*0.25);
	      for (int i=1;i<MAX;i++)
		{
		  ek=m+i*EPS;
		  k=sqrt(ek*ek-m*m);
		  res-=EPS*nB(ek)*k/(ek*ek-p0*p0*0.25);
		}
	      
	     
	      res/=2*M_PI;
	      
	    }
	}     
    }
  else
    printf("D=%i not coded up\n",D);
  
  //printf("resigmat p0=%f p=%f pm2=%f returning %f\n",p0,p,pm2,res);
  
  res/=-M_PI;
  
  return res;
}





double ImPiT(double p0, double p, double m)
{
  double res=0;
  int sig=1;
  if (p0<0)
    {
      sig=-1;
      p0*=-1;
    }
  double pm2=p0*p0-p*p;
  if (D==2)
    {
      //printf("coming in with %f %f %f\n",p0,p,m);
      if (pm2<0)
	{
	  double lowc=-p0*0.5+0.5*p*sqrt(1-4*m*m/pm2);
	  double ek;
	  double c1,c2;
	  //never happens
	  //if (lowc<m)
	  //{
	  //  lowc=m;
	  //  ek=lowc;
	  //     c1=-nB(ek)*0.5/sqrt(-pm2)*log((-pm2)*(2*ek+p0)+sqrt(-pm2)*sqrt(-pm2*(pm2+4*ek*ek+4*ek*p0)-4*m*m*p*p));
	  //    c1-=0.5*EPS*nBprime(ek)*0.5/sqrt(-pm2)*log((-pm2)*(2*ek+p0)+sqrt(-pm2)*sqrt(-pm2*(pm2+4*ek*ek+4*ek*p0)-4*m*m*p*p));
	  // }
	  //else
	  // {
	  ek=lowc;
	  c1=-nB(ek)*0.5/sqrt(-pm2)*log((-pm2)*(2*ek+p0));
	  c1-=0.5*EPS*nBprime(ek)*0.5/sqrt(-pm2)*log((-pm2)*(2*ek+p0));
	  // }
	  //printf("first %f at ek=%f\n",c1,ek);
	  
	  
	  for (int i=1;i<MAX;i++)
	    {
	      ek=i*EPS+lowc;	      
	      c1-=EPS*nBprime(ek)*0.5/sqrt(-pm2)*log((-pm2)*(2*ek+p0)+sqrt(-pm2)*sqrt(-pm2*(pm2+4*ek*ek+4*ek*p0)-4*m*m*p*p));
	      //printf("i=%i ek=%f %f\n",i,ek,c1);
	    }
	  //printf("ek%f c1=%f\n",ek,c1);
	  res+=c1;
	  
	  //second term:
	  
	  lowc=p0*0.5+0.5*p*sqrt(1-4*m*m/pm2);
	  ek=lowc;
	  c2=-nB(ek)*0.5/sqrt(-pm2)*log((-pm2)*(2*ek-p0));
	  c2-=0.5*EPS*nBprime(ek)*0.5/sqrt(-pm2)*log((-pm2)*(2*ek-p0));
	  for (int i=1;i<MAX;i++)
	    {
	      ek=i*EPS+lowc;	      
	      c2-=EPS*nBprime(ek)*0.5/sqrt(-pm2)*log((-pm2)*(2*ek-p0)+sqrt(-pm2)*sqrt(-pm2*(pm2+4*ek*ek-4*ek*p0)-4*m*m*p*p));
	      //printf("i=%i ek=%f %f\n",i,ek,c1);
		}
	  //printf("c2=%f\n",c2);
	  res-=c2;
	}
      else if (pm2>4*m*m)
	{
	  if (p>0)
	    {
	      double highc=p0*0.5+0.5*p*sqrt(1-4*m*m/pm2);
	      double lowc=p0*0.5-0.5*p*sqrt(1-4*m*m/pm2);
	      double meps=(highc-lowc)/(MAX-1);
	      double ek;
	      double c2=0;
	      ek=lowc;	  
	      for (int i=1;i<MAX-1;i++)
		{
		  ek=i*meps+lowc;	      
		  c2-=meps*nBprime(ek)*asin(sqrt((ek-lowc)/(highc-lowc)));
		  //printf("i=%i ek=%f %f\n",i,ek,c1);
		}
	      ek=highc;
	      c2-=0.5*meps*nBprime(ek)*M_PI*0.5;
	      c2+=nB(ek)*M_PI*0.5;
	      //printf("just boundary %f\n",
	      //printf("just boundary %f how %f\n",nB(ek)*M_PI*0.5,c2);
	      res+=c2/sqrt(pm2);
	    }
	  else
	    {
	      //p=0:
	      res=-M_PI*0.5/p0*nB(0.5*p0);
	    }
	}
    }
  else if (D==3)
    {
      if (pm2<0)
	{
	  double b=sqrt(1-4*m*m/pm2);
	  double pplus=0.5*(p0+p*b);
	  double pminus=0.5*(p0-p*b);
	  res=log((1-exp(-pplus))/(1-exp(pminus)));
	  res/=4*p;
	}
      else if (pm2>4*m*m)
	{
	  if (p>0)
	    {
	      double b=sqrt(1-4*m*m/pm2);
	      double pplus=0.5*(p0+p*b);
	      double pminus=0.5*(p0-p*b);
	      res=log((1-exp(-pplus))/(1-exp(-pminus)));
	      res/=4*p;
	    }
	  else
	    {
	      double b=sqrt(1-4*m*m/pm2);
	      res=0.25*b/nB(0.5*p0);
	    }
	}
    }
  else
    printf("D=%i not coded up\n",D);

  res*=sig/M_PI;
  return res;
}



double TanF(int i,int bigN)
{
  double res=tan(0.5*i*M_PI/(bigN));
  return res;
}

double SC(int i,int bigN)
{
  double res=i*1./(bigN-1);
  return res;
}

double IC(int i,int bigN)
{
  double res=1/(1-i*1./bigN);
  return res;
}

double R1int(double p0,double p,double **data,int bigN)
{
  double res=0;
  int p0flag=0;
  if (p0<0)
    {
      p0*=-1;
      p0flag=1;
    }

  
  if (p>0)
    {  //get p bounds
      
      double x=atan(p)*2*bigN/M_PI;
      int lowp=(int)floor(x);
      int highp=(int)ceil(x);
      double y=p0/p*(bigN-1);
      int lowp0=(int)floor(y);
      int highp0=(int)ceil(y);
      //printf("this is R1int p0=%f p=%f\n",p0,p);
      //printf("lowp=%i highp=%i\n",lowp,highp);
      //printf("lowp0=%i highp0=%i\n",lowp0,highp0);

      double x0=TanF(lowp,bigN);
      double x1=TanF(highp,bigN);
      double y0=p*SC(lowp0,bigN);
      double y1=p*SC(highp0,bigN);

      //printf("%f %f %f %f\n",x0,x1,y0,y1);
      
      double f00=data[lowp][lowp0];
      double f01=data[lowp][highp0];
      double f10=data[highp][lowp0];
      double f11=data[highp][highp0];
     
      
      //printf("%f %f %f %f\n",f00,f01,f10,f11);
     
      res=bilinear_int(x0,x1,y0,y1,f00,f10,f11,f01,p,p0);
      if (p0flag)
	res*=-1;
      
    }
  else
    {
      //nothing, return zero
    }
  return res;
}

double R2int(double p0,double p,double m,double **data,int bigN)
{
  double res=0;
  int p0flag=0;
  if (p0<0)
    {
      p0*=-1;
      p0flag=1;
    }
  

  if (p>0)    
    {
      double maxp0=sqrt(p*p+4*m*m)*IC(bigN-1,bigN);
      if (p0<maxp0)
	{
	  //get p bounds
	  double x=atan(p)*2*bigN/M_PI;
	  int lowp=(int)floor(x);
	  int highp=(int)ceil(x);
	  double y=(1-sqrt(p*p+4*m*m)/p0)*bigN;
	  int lowp0=(int)floor(y);
	  int highp0=(int)ceil(y);
	  //printf("lowp=%i highp=%i\n",lowp,highp);
	  //printf("lowp0=%i highp0=%i\n",lowp0,highp0);
	  double f00=data[lowp][lowp0];
	  double f01=data[lowp][highp0];
	  double f10=data[highp][lowp0];
	  double f11=data[highp][highp0];
	  double x0=TanF(lowp,bigN);
	  double x1=TanF(highp,bigN);
	  double y0=sqrt(p*p+4*m*m)*IC(lowp0,bigN);
	  double y1=sqrt(p*p+4*m*m)*IC(highp0,bigN);
	  
	  //printf("%f %f %f %f\n",f00,f01,f10,f11);
	  //printf("%f %f %f %f\n",x0,x1,y0,y1);
	  res=bilinear_int(x0,x1,y0,y1,f00,f10,f11,f01,p,p0);

	  if (isnan(res))
	    {
	      printf("R2int debug\n");
	      printf("%f %f %f %f\n",f00,f01,f10,f11);
	      printf("%f %f %f %f\n",x0,x1,y0,y1);
	    }
	  
	  if (p0flag)
	    res*=-1;	  
	}
      else
	{
	  res=0;
	}
    }
  else
    {
      //nothing, return zero
    }
  return res;
}

double ImInt(double p0,double p,double m,double **d1,double **d2,int bigN)
{
  double res=0;
  double p2=p*p-p0*p0;
  //printf("this is imint p2=%f\n",p2);

  double maxp=TanF(bigN-1,bigN);

  if (p<maxp)
    {
      if (p2>0)
	{
	  res=R1int(p0,p,d1,NN);
	}
      else if (p2<-4*m*m)
	{
	  res=R2int(p0,p,m,d2,NN);
	}
      //else nothing to do
    }
  else
    {
      res=0; //above interpolation limit
    }

  if (isnan(res))
    {
      printf("ImInt returning nan; p2=%f\n",p2);
    }

  return res;
}

double Reint(double p0,double p,double **data,int bigN)
{
  double res=0;
  if (p0<0)
    p0*=-1;


  double maxp=TanF(bigN-1,bigN);
  if ((p<maxp)&&(p0<maxp))
    {
      double x=atan(p)*2*bigN/M_PI;
      int lowp=(int)floor(x);
      int highp=(int)ceil(x);
      double y=atan(p0)*2*bigN/M_PI;
      int lowp0=(int)floor(y);
      int highp0=(int)ceil(y);
      
      double f00=data[lowp][lowp0];
      double f01=data[lowp][highp0];
      double f10=data[highp][lowp0];
      double f11=data[highp][highp0];
      double x0=TanF(lowp,bigN);
      double x1=TanF(highp,bigN);
      double y0=TanF(lowp0,bigN);
      double y1=TanF(highp0,bigN);
      
      //printf("%f %f %f %f\n",f00,f01,f10,f11);
      //printf("%f %f %f %f\n",x0,x1,y0,y1);
      res=bilinear_int(x0,x1,y0,y1,f00,f10,f11,f01,p,p0);
    }
  else
    {
      //out of bounds
      res=0;
    }

  return res;
}

void interpolate_ImPi(int NN,double m,double** R1,double **R2)
{
  //region 1: p^2-p0^2<0; p0=x*p
  for (int i=0;i<NN;i++)
    for (int j=0;j<NN;j++)
      {
	//just region 1 values
	double p=TanF(i,NN);
	double p0=p*SC(j,NN);
	double p2=-p0*p0+p*p;
	R1[i][j]=ImPiT(p0,p,m);
	//printf("%i %i p0=%f p=%f p2=%f impi=%.16g\n",i,j,p0,p,p2,dataR1[i][j]);
      }

  //region2 : p^2-p0^2<-4m^2; p0=sqrt(p^2+4m^2)/(1-x)
  for (int i=0;i<NN;i++)
    for (int j=0;j<NN;j++)
      {
	//just region 1 values
	double p=TanF(i,NN);
	double p0=sqrt(p*p+4*m*m)*IC(j,NN);
	double p2=-p0*p0+p*p;
	R2[i][j]=ImPiT(p0,p,m);
	//printf("%i %i p0=%f p=%f p2=%f impi=%.16g\n",i,j,p0,p,p2,dataR2[i][j]);
      }
}

void interpolate_RePi(int NN,double m,double** data)
{
  for (int i=0;i<NN;i++)
    for (int j=0;j<NN;j++)
      {
	double p=TanF(i,NN);
	double p0=TanF(j,NN);
	double p2=-p0*p0+p*p;
	data[i][j]=RePiT(p0,p,m);
	//printf("%i %i p0=%f p=%f p2=%f repi=%.16g\n",i,j,p0,p,p2,data[i][j]);
      }
}

void interpolate_allPi(int NN,double m,double** R1,double **R2,double **Re)
{
  






#pragma omp parallel for
  for (int i=0;i<NN;i++)
    for (int j=0;j<NN;j++)
      {
	//just region 1 values
	double p=TanF(i,NN);
	double p0=p*SC(j,NN);
	//double p2=-p0*p0+p*p;
	R1[i][j]=ImPiT(p0,p,m);
	//printf("%i %i p0=%f p=%f p2=%f impi=%.16g\n",i,j,p0,p,p2,dataR1[i][j]);
	
	//region 2
	p0=sqrt(p*p+4*m*m)*IC(j,NN);
	R2[i][j]=ImPiT(p0,p,m);
	
	//real part
	p0=TanF(j,NN);
	Re[i][j]=RePiT(p0,p,m);
	if (isnan(Re[i][j]))
	  printf("nan for p0=%f p=%f\n",p0,p);
      }
  
}

