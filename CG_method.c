#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void transform(double**,double**,int);

main(int argc, char*argv[])
{
	int i,j,N,count;
	double L,S,alpha,**T,**b,**r,**p,**ptemp,dx,sp,T1,T2,tol,temp,error,step,num,den;

	sscanf(argv[1],"%d",&N);
	sscanf(argv[2],"%lf",&L);
	sscanf(argv[3],"%lf",&S);
	sscanf(argv[4],"%lf",&alpha);
	sscanf(argv[5],"%lf",&T1);
	sscanf(argv[6],"%lf",&T2);
	sscanf(argv[7],"%lf",&tol);

	dx=L/N;
	sp=S*dx*dx/alpha;
	T=(double**)malloc(N*sizeof(double*));
	b=(double**)malloc(N*sizeof(double*));
	r=(double**)malloc(N*sizeof(double*));
	p=(double**)malloc(N*sizeof(double*));
	ptemp=(double**)malloc(N*sizeof(double*));
	for(i=0;i<N;i++)
	{
		T[i]=(double*)malloc(N*sizeof(double));
		b[i]=(double*)malloc(N*sizeof(double));
		r[i]=(double*)malloc(N*sizeof(double));
		p[i]=(double*)malloc(N*sizeof(double));
		ptemp[i]=(double*)malloc(N*sizeof(double));
		for (j = 0; j < N; j++)
		{
			T[i][j]=T2;
			b[i][j]=-sp;
		}
		b[i][0]-=T1;
		b[i][N-1]-=T2;
	}

	transform(T,r,N);
	
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			r[i][j]=b[i][j]-r[i][j];
		}
	}
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			p[i][j]=r[i][j];
		}
	}

	count=0;
	error=10*tol;
	while(error>tol)
	{
		step=0;
		temp=0;
		den=0;
		transform(p,ptemp,N);
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				step+=r[i][j]*r[i][j];
				den+=p[i][j]*ptemp[i][j];
				temp+=p[i][j]*ptemp[i][j];
			}
		}
		step/=temp;
		error=0;
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				T[i][j]+=step*p[i][j];
				r[i][j]-=step*ptemp[i][j];
				error+=r[i][j]*r[i][j];
			}
		}
		num=0;
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				num+=r[i][j]*r[i][j];
			}
		}
		num/=den;
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				p[i][j]=r[i][j]+num*p[i][j];
			}
		}
		printf("Iteration no: %d\tError: %.14f\tStepsize: %.14f\n",++count,error,step);
	}

	FILE*f_out;

	f_out=fopen("T_CG.txt","w");
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			fprintf(f_out, "%lf\t",T[i][j]);
		}
		fprintf(f_out, "\n");
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