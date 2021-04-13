double est_root(double xa,double xb, double ya, double yb)
{
  double res=(xb*ya-xa*yb)/(ya-yb);
  return res;
}

double convergeroot1(double (*fun)(double *),double *args)
{
  double xa,xb,ya,yb,res;
  //initial gues:
  xa=args[0];
  int IT=round(args[1]);
  double strength=args[2];
  double brek=args[3];
  int verbose=round(args[4]);
  double mult;
  for (int i=0;i<IT;i++)
    {
      args[0]=xa;
      ya=(*fun)(args);      
      mult=-ya*strength;
      xb=xa+mult;
      args[0]=xb;
      yb=(*fun)(args);
      res=est_root(xa,xb,ya,yb);
      if (verbose)
	printf("i=%i res=%f y(%f)=%f y(%f)=%f\n",i,res,xa,ya,xb,yb);
      xa=xb+(res-xb)*brek;
      
    }
  return res;
}

int sign(double x)
{
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

double convergeroot2(double (*fun)(double *),double *args)
{
  
  int IT=round(args[1]);
  double low=args[2];
  double high=args[3];  
  int verbose=round(args[4]);
  if (low>high)
    {
      printf("converge root 2 WARNING: low>high!\n");
      //ABORT=1;
    }

  
  args[0]=low;
  double ylow=(*fun)(args);
  
  args[0]=high;
  double yhigh=(*fun)(args);

  if (sign(ylow)==sign(yhigh))
    {
      printf("converge root 2 WARNING: root not bracketed!\n");
      printf("low=%f: ylow%f high=%f yhigh %f\n",low,ylow,high,yhigh);
      //ABORT=1;
    }
 
  double middle;
  for (int i=0;i<IT;i++)
    {
      
      
      middle=0.5*(low+high);
      args[0]=middle;
      double ym=(*fun)(args);
      
      if (sign(ym)==sign(ylow))
	{
	  low=middle;
	  ylow=ym;
	}
      else if (sign(ym)==sign(yhigh))
	{
	  high=middle;
	  yhigh=ym;
	}
	  
      if (verbose)
	printf("i=%i %f %f %f %f\n",i,low,ylow,high,yhigh);
      
    }
  return middle;
}


double testf(double *args)
{
  double x=args[0];
  return (x-0.144)*(x-0.144);
}
