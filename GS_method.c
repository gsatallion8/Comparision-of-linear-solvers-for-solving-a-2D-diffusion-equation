#include <stdio.h>
#include <stdlib.h>
#include <string.h>

main(int argc, char*argv[])
{
	int i,j,N,count;
	double L,S,alpha,**T,dx,b,T1,T2,tol,temp,error;

	sscanf(argv[1],"%d",&N);
	sscanf(argv[2],"%lf",&L);
	sscanf(argv[3],"%lf",&S);
	sscanf(argv[4],"%lf",&alpha);
	sscanf(argv[5],"%lf",&T1);
	sscanf(argv[6],"%lf",&T2);
	sscanf(argv[7],"%lf",&tol);

	dx=L/N;
	b=S*dx*dx/alpha;
	T=(double**)malloc(N*sizeof(double*));
	for(i=0;i<N;i++)
	{
		T[i]=(double*)malloc(N*sizeof(double));
		for (j = 0; j < N; j++)
		{
			T[i][j]=T2;
		}
	}

	count=0;
	error=10*tol;
	while(error>tol)
	{
		printf("Iteration no: %d\t error: %lf\n",++count,error);
		error=0;
		temp=(T1+b+2*T[1][0]+T[0][1])/4;
		error+=(T[0][0]-temp)*(T[0][0]-temp);
		T[0][0]=temp;

		for(j=1;j<N-1;j++)
		{
			temp=(b+2*T[1][j]+T[0][j+1]+T[0][j-1])/4;
			error+=(T[0][j]-temp)*(T[0][j]-temp);
			T[0][j]=temp;
		}

		temp=(T2+b+2*T[1][N-1]+T[0][N-2])/4;
		error+=(T[0][N-1]-temp)*(T[0][N-1]-temp);
		T[0][N-1]=temp;

		for(i=1;i<N-1;i++)
		{
			temp=(T1+b+T[i+1][0]+T[i-1][0]+T[i][1])/4;
			error+=(T[i][0]-temp)*(T[i][0]-temp);
			//printf("error: %.14f\t %d\t %d\n",(T[i][0]-temp)*(T[i][0]-temp),i,0);
			T[i][0]=temp;

			for(j=1;j<N-1;j++)
			{
				temp=(b+T[i+1][j]+T[i-1][j]+T[i][j+1]+T[i][j-1])/4;
				error+=(T[i][j]-temp)*(T[i][j]-temp);
				//printf("error: %.14f\t %d\t %d\n",(T[i][j]-temp)*(T[i][j]-temp),i,j);
				T[i][j]=temp;
			}

			temp=(T2+b+T[i+1][N-1]+T[i-1][N-1]+T[i][N-2])/4;
			error+=(T[i][N-1]-temp)*(T[i][N-1]-temp);
			//printf("error: %.14f\t %d\t %d\n",(T[i][N-1]-temp)*(T[i][N-1]-temp),i,N-1);
			T[i][N-1]=temp;
		}

		temp=(T1+b+2*T[N-2][0]+T[N-1][1])/4;
		error+=(T[N-1][0]-temp)*(T[N-1][0]-temp);
		//printf("error: %.14f\t %d\t %d\n",(T[N-1][0]-temp)*(T[N-1][0]-temp),N-1,0);
		T[N-1][0]=temp;

		for(j=1;j<N-1;j++)
		{
			temp=(b+2*T[N-2][j]+T[N-1][j+1]+T[N-1][j-1])/4;
			error+=(T[N-1][j]-temp)*(T[N-1][j]-temp);
			//printf("error: %.14f\t %d\t %d\n",(T[N-1][j]-temp)*(T[N-1][j]-temp),N-1,j);
			T[N-1][j]=temp;
		}

		temp=(T2+b+2*T[N-2][N-1]+T[N-1][N-2])/4;
		error+=(T[N-1][N-1]-temp)*(T[N-1][N-1]-temp);
		//printf("error: %.14f\t %d\t %d\n",(T[N-1][N-1]-temp)*(T[N-1][N-1]-temp),N-1,N-1);
		T[N-1][N-1]=temp;
	}

	FILE*f_out;

	f_out=fopen("T_GS.txt","w");
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