#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void transform(double**,double**,int);

main(int argc, char*argv[])
{
	int i,j,N,count,k1,k2;
	double L,S,alpha,**T,**b,**r,**p,dx,sp,T1,T2,tol,temp,error,step,A;
	FILE*f_out;
	double PI=3.14159265;

	f_out=fopen("GD_IG.txt","w");

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
	b=(double**)malloc(N*sizeof(double*));
	r=(double**)malloc(N*sizeof(double*));
	p=(double**)malloc(N*sizeof(double*));
	for(i=0;i<N;i++)
	{
		T[i]=(double*)malloc(N*sizeof(double));
		b[i]=(double*)malloc(N*sizeof(double));
		r[i]=(double*)malloc(N*sizeof(double));
		p[i]=(double*)malloc(N*sizeof(double));
		for (j = 0; j < N; j++)
		{
			//T[i][j]=T1;
			T[i][j]=T1+(T2-T1)*(float)(j+1)/(N+1)+A*sin(k1*PI*(j+1)/(N+1))*cos(k2*PI*(i+0.5)/(N+0.5));
			b[i][j]=-sp;
		}
		b[i][0]-=T1;
		b[i][N-1]-=T2;
	}
	
	count=0;
	error=10*tol;
	while(error>tol)
	{
		step=0;
		transform(T,r,N);
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				r[i][j]=b[i][j]-r[i][j];
				step+=r[i][j]*r[i][j];
				error+=r[i][j]*r[i][j];
			}
		}
		//printf("%.14f\t",step);
		transform(r,p,N);
		temp=0;
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				temp+=r[i][j]*p[i][j];
			}
		}
		step/=temp;
		//printf("%.14f\n",temp);

		error=0;
		//step=-0.1;
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				T[i][j]+=step*r[i][j];
				error+=r[i][j]*r[i][j];
				//error+=r[i][j]*r[i][j];
			}
		}
		error/=N;

		printf("Iteration no: %d\tError: %.14f\tStepsize: %.14f\n",++count,error,step);
		fprintf(f_out,"%d\t%.14f\n",count,error);
	}
	
	fclose(f_out);
}

void transform(double**x,double**y,int N)
{
	int i,j;

	y[0][0]=2*x[1][0]+x[0][1]-4*x[0][0];
	for(j=1;j<N-1;j++)
	{
		y[0][j]=2*x[1][j]+x[0][j+1]+x[0][j-1]-4*x[0][j];
	}
	y[0][N-1]=2*x[1][N-1]+x[0][N-2]-4*x[0][N-1];

	for(i=1;i<N-1;i++)
	{
		y[i][0]=x[i+1][0]+x[i-1][0]+x[i][1]-4*x[i][0];
		for(j=1;j<N-1;j++)
		{
			y[i][j]=x[i+1][j]+x[i-1][j]+x[i][j+1]+x[i][j-1]-4*x[i][j];
		}
		y[i][N-1]=x[i+1][N-1]+x[i-1][N-1]+x[i][N-2]-4*x[i][N-1];
	}

	y[N-1][0]=2*x[N-2][0]+x[N-1][1]-4*x[N-1][0];
	for(j=1;j<N-1;j++)
	{
		y[N-1][j]=2*x[N-2][j]+x[N-1][j+1]+x[N-1][j-1]-4*x[N-1][j];
	}
	y[N-1][N-1]=2*x[N-2][N-1]+x[N-1][N-2]-4*x[N-1][N-1];
}