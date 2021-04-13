void getquadrature(double *x,double *w,int size,double m,int verbose)
{
  double ints[2*size];
  //printf("my %f\n",polylog(2,exp(-1.)));

  if (D==2)
    {
      for (int i=0;i<2*size;i++)
	{
	  ints[i]=pow(m,i)*nB(m);
	  if (i>0)
	    {
	      ints[i]-=i*pow(m,i-1)*log(1-exp(-m));
	      long long fac=i*(i-1);
	      for (int j=i-2;j>=0;j--)
		{
		  ints[i]+=fac*pow(m,j)*polylog(i-j,exp(-m));
		  //printf("i=%i fac=%llu, pow=%i, poly=%i ints=%f\n",i,fac,j,i-j,ints[i]);
		  fac*=j;	      
		}
	    }
	  //printf("i=%i int[%i]=%f\n",i,i,ints[i]);
	}
    }
  else if (D==3)
    {
      for (int i=0;i<2*size;i++)
	ints[i]=D3ints(i,m);
      
    }
  else
    printf("D=%i not coded  up!\n",D);
  
  switch(size)
    {
    case 1:
      {
	x[0]=ints[1]/ints[0];
	w[0]=ints[0];
	break;
      }
    case 2:
      {
	double a20=(ints[2]*ints[2]-ints[1]*ints[3])/(ints[1]*ints[1]-ints[0]*ints[2]);
	double a21=(-ints[1]*ints[2]+ints[0]*ints[3])/(ints[1]*ints[1]-ints[0]*ints[2]);
	
	x[0]=-0.5*a21-sqrt(a21*a21*0.25-a20);
	x[1]=-0.5*a21+sqrt(a21*a21*0.25-a20);
	w[0]=(ints[0]*x[1]-ints[1])/(x[1]-x[0]);
	w[1]=(ints[1]-ints[0]*x[0])/(x[1]-x[0]);
	break;
      }
    case 3:
      {
	double det=ints[2]*ints[2]*ints[2]-2*ints[1]*ints[2]*ints[3]+ints[0]*ints[3]*ints[3]+ints[1]*ints[1]*ints[4]-ints[0]*ints[2]*ints[4];
	double a32=-(ints[3]*(ints[2]*ints[2]-ints[1]*ints[3])+ints[4]*(-ints[1]*ints[2]+ints[0]*ints[3])+ints[5]*(ints[1]*ints[1]-ints[0]*ints[2]))/det;
	double a31=-(ints[3]*(-ints[2]*ints[3]+ints[1]*ints[4])+ints[4]*(ints[2]*ints[2]-ints[0]*ints[4])+ints[5]*(-ints[1]*ints[2]+ints[0]*ints[3]))/det;
	double a30=-(ints[3]*(ints[3]*ints[3]-ints[2]*ints[4])+ints[4]*(-ints[2]*ints[3]+ints[1]*ints[4])+ints[5]*(ints[2]*ints[2]-ints[1]*ints[3]))/det;

	double p=a31-a32*a32/3.;
	double q=a30-a31*a32/3.+2*a32*a32*a32/27.;
	x[0]=-a32/3.+2*sqrt(-p/3.)*cos(acos(3*q*0.5/p*sqrt(-3/p))/3.-2*M_PI/3.*0);
	x[1]=-a32/3.+2*sqrt(-p/3.)*cos(acos(3*q*0.5/p*sqrt(-3/p))/3.-2*M_PI/3.*1);
	x[2]=-a32/3.+2*sqrt(-p/3.)*cos(acos(3*q*0.5/p*sqrt(-3/p))/3.-2*M_PI/3.*2);
	w[0]=(ints[2]+ints[0]*x[1]*x[2]-ints[1]*(x[1]+x[2]))/(x[0]-x[1])/(x[0]-x[2]);
	w[1]=-(ints[2]+ints[0]*x[0]*x[2]-ints[1]*(x[0]+x[2]))/(x[0]-x[1])/(x[1]-x[2]);
	w[2]=-(ints[2]+ints[0]*x[0]*x[1]-ints[1]*(x[0]+x[1]))/(x[0]-x[2])/(x[2]-x[1]);
	//printf("a32 =%f a31=%f a30=%f\n",a32,a31,a30);
	//printf("get %f\n",ints[3]+a32*ints[2]+a31*ints[1]+a30*ints[0]);
	//for (int i=0;i<size;i++)
	//	  printf("x[%i]=%f\n",i,x[i]);
			      
	break;
      }
    default:
      {
	MatrixXd A = MatrixXd::Random(size,size);
	VectorXd b = VectorXd::Random(size);

	
	for (int i=0;i<size;i++)
	  {
	    b(i)=ints[size+i];
	    for (int j=0;j<size;j++)
	      A(i,j)=ints[size-j-1+i];
	  }
	//cout << A << endl;
	
	VectorXd out =  A.colPivHouseholderQr().solve(b);
	
	
	//cout << "The coefficients are " << out << endl;

	//double err=(A*out-b).norm()/b.norm();

	//cout << "Relative error is \n" << err << endl;

	//out =  A.completeOrthogonalDecomposition().solve(b);
	out =  A.partialPivLu().solve(b);
	//cout << "Second choice " << out << endl;
	double err=(A*out-b).norm()/b.norm();
	if (verbose)
	  printf("\t Quadrature: relative error is %.2g\n",err);

	
	//double relative_error = (A*x - b).norm() / b.norm(); // norm() is L2 norm
	//cout << "The relative error is:\n" << relative_error << endl;

	//create companion matrix:
	for (int i=0;i<size;i++)
	  for (int j=0;j<size;j++)
	    {	      
	      A(i,j)=0;
	      if (j==size-1)
		A(i,j)=out(size-1-i);
	      if (j==i-1)
		A(i,j)=1;
	    }
	//cout << "The companion matrix is " << A << endl;

	EigenSolver<MatrixXd> es(A);

	out = es.eigenvalues().real();

	//cout << "The eigenvalues are " << out << endl;
	for (int i=0;i<size;i++)
	  x[i]=out(i);

	for (int i=0;i<size;i++)
	  for (int j=0;j<size;j++)
	    {
	      A(i,j)=pow(x[j],i);
	    }

	for (int i=0;i<size;i++)
	  b(i)=ints[i];

	VectorXd weights =  A.colPivHouseholderQr().solve(b);
	
	//cout << "The solution is " << weights << endl;
	for (int i=0;i<size;i++)
	  w[i]=weights(i);

	if (verbose>1)
	  {
	    cout << "Stencils are " << out <<endl;
	    cout << "Weights are " << weights <<endl;
	  }
	  
	
	break;
      }
    }

  
  
  
  if (verbose)
    {
      double check;

      for (int j=0;j<2*size-1;j++)
	{
	  check=0;
	  for (int i=0;i<size;i++)
	    check+=w[i]*pow(x[i],j);
	  
	  printf("quadrature check n=%i %.16g, diff=%.16g rdiff=%f\n",j,check,ints[j]-check,(ints[j]-check)/ints[j]);
	}
      /*
      check=0;
      for (int i=0;i<size;i++)
	check+=w[i]*x[i];
      
      printf("quadrature check  %f=?%f\n",check,ints[1]);
      
      check=0;
      for (int i=0;i<size;i++)
	check+=w[i]*x[i]*x[i];
      
      printf("quadrature check  %f=?%f\n",check,ints[2]);
      
      check=0;
      for (int i=0;i<size;i++)
	check+=w[i]*x[i]*x[i]*x[i];
      
	printf("quadrature check  %f=?%f\n",check,ints[3]);*/
    }
}



double poly_eval(double ** polys,int which,double x)
{
  double poly=0;
  for (int j=0;j<quad;j++)
    poly+=polys[which][j]*pow(x,j);
  return poly;
}

void getpolys(double *x,double *w,double** polys,int size,double m,int verbose)
{
  double ints[2*size];

  if (D==2)
    {
      for (int i=0;i<2*size;i++)
	{
	  ints[i]=pow(m,i)*nB(m);
	  if (i>0)
	    {
	      ints[i]-=i*pow(m,i-1)*log(1-exp(-m));
	      long long  fac=i*(i-1);
	      for (int j=i-2;j>=0;j--)
		{
		  ints[i]+=fac*pow(m,j)*polylog(i-j,exp(-m));
		  //printf("i=%i fac=%i, pow=%i, poly=%i\n",i,fac,j,i-j);
		  fac*=j;	      
		}
	    }
	}
    }
  else if (D==3)
    {
      for (int i=0;i<2*size;i++)
	ints[i]=D3ints(i,m);
      
    }
  else printf("D=%i not coded up\n",D);
    

  for (int i=0;i<size;i++)    
    for (int j=0;j<size;j++)
      polys[i][j]=0;
  
  

  //lowest order poly:
  polys[0][0]=1;
  
  double norm=0;
  for (int i=0;i<size;i++)
    {
      for (int j=0;j<size;j++)
	{
	  norm+=w[i]*polys[0][j]*pow(x[i],j);
	  //printf("i=%i j=%i pol=%f add=%f\n",i,j,polys[0][j],polys[0][j]*pow(x[i],size-j-1));
	}
    }
  
  for (int j=0;j<size;j++)
    polys[0][j]/=sqrt(ints[0]);
  
  for (int od=1;od<size;od++)
    {
      if (verbose)
	printf("\t\t INFO: Generating orthogonal polynomial P_%i\n",od);
      MatrixXd A = MatrixXd::Random(od,od);
      VectorXd b = VectorXd::Random(od);
      for (int i=0;i<od;i++)
	{
	  b(i)=ints[od+i];
	  for (int j=0;j<od;j++)
	    A(i,j)=ints[od-j-1+i];
	}
      VectorXd out =  A.partialPivLu().solve(b);
      //cout << "The coefficients are " << out << endl;
      polys[od][od]=1;
      for (int i=0;i<od;i++)
	polys[od][i]=-out(od-1-i);
    }

  for (int od=1;od<size;od++)
    {
      norm=0;
      for (int i=0;i<size;i++)
	{
	  double pol=0;
	  for (int j=0;j<size;j++)
	    pol+=polys[od][j]*pow(x[i],j);
	  norm+=w[i]*pol*pol;
	  //printf("i=%i j=%i pol=%f add=%f\n",i,j,polys[1][j],polys[1][j]*pow(x[i],j));
	  
	}
      for (int j=0;j<size;j++)
	polys[od][j]/=sqrt(norm);
    }
  
  if (verbose>1)
    {
      for (int a=0;a<size;a++)
	for (int b=0;b<size;b++)
	  {
	    norm=0;
	    for (int i=0;i<size;i++)
	      {
		double pol1=0;
		double pol2=0;
		for (int j=0;j<size;j++)
		  {
		    pol1+=polys[a][j]*pow(x[i],j);
		    pol2+=polys[b][j]*pow(x[i],j);
		  }
		norm+=w[i]*pol1*pol2;
		//norm+=w[i]*pol1*poly_eval(polys,b,size,x[i]);
		//printf("i=%i j=%i pol=%f add=%f\n",i,j,polys[1][j],polys[1][j]*pow(x[i],j));	 
	      }
	    printf("{P[%i]P[%i]}= %f\n",a,b,norm);
	  }
    }
}


