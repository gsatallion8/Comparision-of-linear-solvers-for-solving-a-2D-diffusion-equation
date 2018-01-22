#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void TDMA(double*,double*,double*,double*,double*,int);

main(int argc, char*argv[])
{
	int i,j,N,count,k1,k2;
	double L,S,alpha,**T,dx,sp,T1,T2,tol,**temp,error,*tempx,*tempy,*a,*b,*c,*d,ferror,A;
	FILE*f_out;
	double PI=3.14159265;

	f_out=fopen("lineTDMA_IG.txt","w");

	sscanf(argv[1],"%d",&N);
	sscanf(argv[2],"%lf",&L);
	sscanf(argv[3],"%lf",&S);
	sscanf(argv[4],"%lf",&alpha);
	sscanf(argv[5],"%lf",&T1);
	sscanf(argv[6],"%lf",&T2);
	sscanf(argv[7],"%lf",&tol);
	sscanf(argv[8],"%lf",&A);
	sscanf(argv[9],"%d",&k1);
	sscanf(argv[10],"%d",&k2);

	dx=L/N;
	sp=S*dx*dx/alpha;
	T=(double**)malloc(N*sizeof(double*));
	temp=(double**)malloc(N*sizeof(double*));
	tempx=(double*)malloc(N*sizeof(double));
	tempy=(double*)malloc(N*sizeof(double));
	a=(double*)malloc(N*sizeof(double));
	b=(double*)malloc(N*sizeof(double));
	c=(double*)malloc(N*sizeof(double));
	d=(double*)malloc(N*sizeof(double));
	for(i=0;i<N;i++)
	{
		T[i]=(double*)malloc(N*sizeof(double));
		temp[i]=(double*)malloc(N*sizeof(double));
		for (j = 0; j < N; j++)
		{
			T[i][j]=T1+(T2-T1)*(float)(j+1)/(N+1)+A*sin(k1*PI*(j+1)/(N+1))*cos(k2*PI*(i+0.5)/(N+0.5));
		}
	}
	ferror=10*tol;
	count=0;
	while(ferror>tol)
	{
		printf("Iteration no: %d\tError: %.14f\tferror: %.14f\n",count++,error,ferror);


		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				temp[i][j]=T[i][j];
			}
		}

		//Horizontal Sweep:
		printf("Forward horizontal sweep\n");
		//i=0
		for(j=0;j<N;j++)
		{
			a[j]=-4;
			c[j]=1;
			b[j]=1;
			d[j]=-2*T[1][j]-sp;
		}
		d[0]-=T1;
		d[N-1]-=T2;
		TDMA(a,b,c,d,tempy,N);
		for(j=0;j<N;j++)
		{
			T[0][j]=tempy[j];
		}

		//i=1 to N-2
		for(i=1;i<N-1;i++)
		{
			for(j=0;j<N;j++)
			{
				a[j]=-4;
				c[j]=1;
				b[j]=1;
				d[j]=-T[i+1][j]-T[i-1][j]-sp;
			}
			d[0]-=T1;
			d[N-1]-=T2;
			TDMA(a,b,c,d,tempy,N);
			for(j=0;j<N;j++)
			{
				T[i][j]=tempy[j];
			}
		}

		//i=N-1
		for(j=0;j<N;j++)
		{
			a[j]=-4;
			c[j]=1;
			b[j]=1;
			d[j]=-2*T[N-2][j]-sp;
		}
		d[0]-=T1;
		d[N-1]-=T2;
		TDMA(a,b,c,d,tempy,N);
		for(j=0;j<N;j++)
		{
			T[N-1][j]=tempy[j];
		}

		//Horizontal sweep:
		printf("Backward horizontal sweep\n");
		//i=N-1
		for(j=0;j<N;j++)
		{
			a[j]=-4;
			c[j]=1;
			b[j]=1;
			d[j]=-2*T[N-2][j]-sp;
		}
		d[0]-=T1;
		d[N-1]-=T2;
		TDMA(a,b,c,d,tempy,N);
		for(j=0;j<N;j++)
		{
			T[N-1][j]=tempy[j];
		}

		//i=N-2 to 1
		for(i=N-2;i>0;i--)
		{
			for(j=0;j<N;j++)
			{
				a[j]=-4;
				c[j]=1;
				b[j]=1;
				d[j]=-T[i+1][j]-T[i-1][j]-sp;
			}
			d[0]-=T1;
			d[N-1]-=T2;
			TDMA(a,b,c,d,tempy,N);
			for(j=0;j<N;j++)
			{
				T[i][j]=tempy[j];
			}
		}

		//i=0
		for(j=0;j<N;j++)
		{
			a[j]=-4;
			c[j]=1;
			b[j]=1;
			d[j]=-2*T[1][j]-sp;
		}
		d[0]-=T1;
		d[N-1]-=T2;
		TDMA(a,b,c,d,tempy,N);
		for(j=0;j<N;j++)
		{
			T[0][j]=tempy[j];
		}

		//Vertical sweep:
		printf("Forward vertical sweep\n");
		//j=0
		for(i=0;i<N;i++)
		{
			a[i]=-4;
			c[i]=1;
			b[i]=1;
			d[i]=-T1-T[i][1]-sp;
		}
		b[0]+=1;
		c[N-1]+=1;
		TDMA(a,b,c,d,tempx,N);
		for(i=0;i<N;i++)
		{
			T[i][0]=tempx[i];
		}

		//j=1 to N-2
		for (j=1;j<N-1;j++)
		{
			for(i=0;i<N;i++)
			{
				a[i]=-4;
				c[i]=1;
				b[i]=1;
				d[i]=-T[i][j+1]-T[i][j-1]-sp;
			}
			b[0]+=1;
			c[N-1]+=1;
			TDMA(a,b,c,d,tempx,N);
			for(i=0;i<N;i++)
			{
				T[i][j]=tempx[i];
			}	
		}

		//j=N-1
		for(i=0;i<N;i++)
		{
			a[i]=-4;
			c[i]=1;
			b[i]=1;
			d[i]=-T2-T[i][N-2]-sp;
		}
		b[0]+=1;
		c[N-1]+=1;
		TDMA(a,b,c,d,tempx,N);
		for(i=0;i<N;i++)
		{
			T[i][N-1]=tempx[i];
		}

		//Vertical sweep:
		printf("Backward vertical sweep\n");
		//j=N-1
		for(i=0;i<N;i++)
		{
			a[i]=-4;
			c[i]=1;
			b[i]=1;
			d[i]=-T2-T[i][N-2]-sp;
		}
		b[0]+=1;
		c[N-1]+=1;
		TDMA(a,b,c,d,tempx,N);
		for(i=0;i<N;i++)
		{
			T[i][N-1]=tempx[i];
		}

		//j=N-2 to 1
		for (j=N-2;j>0;j--)
		{
			for(i=0;i<N;i++)
			{
				a[i]=-4;
				c[i]=1;
				b[i]=1;
				d[i]=-T[i][j+1]-T[i][j-1]-sp;
			}
			b[0]+=1;
			c[N-1]+=1;
			TDMA(a,b,c,d,tempx,N);
			for(i=0;i<N;i++)
			{
				T[i][j]=tempx[i];
			}	
		}

		//j=0
		for(i=0;i<N;i++)
		{
			a[i]=-4;
			c[i]=1;
			b[i]=1;
			d[i]=-T1-T[i][1]-sp;
		}
		b[0]+=1;
		c[N-1]+=1;
		TDMA(a,b,c,d,tempx,N);
		for(i=0;i<N;i++)
		{
			T[i][0]=tempx[i];
		}

		//Error computation:
		error=0;
		ferror=0;
		for (i = 0; i < N; i++)
		{
			for(j=0;j<N;j++)
			{
				error+=(temp[i][j]-T[i][j])*(temp[i][j]-T[i][j]);
				ferror+=(T1+(T2-T1)*(float)(j+1)/(N+1)-T[i][j])*(T1+(T2-T1)*(float)(j+1)/(N+1)-T[i][j]);
			}
		}
		ferror/=N;
		fprintf(f_out,"%d\t%.14f\n",count,ferror);
	}
	fclose(f_out);
}

void TDMA(double*a,double*b,double*c,double*d,double*x,int n)
{
	int i;
    double*b_,*d_;
    b_=(double*)malloc(n*sizeof(double));
    d_=(double*)malloc(n*sizeof(double));
    b_[0]=b[0]/a[0];
    d_[0]=d[0]/a[0];
    for(i=1;i<n;i++)
    {
    	b_[i]=b[i]/(a[i]-c[i]*b_[i-1]);
    	d_[i]=(d[i]-c[i]*d_[i-1])/(a[i]-c[i]*b_[i-1]);
    }
    x[n-1]=d_[n-1];
    for(i=n-2;i>=0;i--)
    {
        x[i]=d_[i]-b_[i]*x[i+1];
    }
    return;
}