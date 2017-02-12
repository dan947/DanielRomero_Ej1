#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int L = 5, l = 2, d = 1, V0 = 100, m, N;
double h = 0.01953125;

int transformer(int i, int j);
double *init(int x0, int x1, int y0, int y1, double *array);
double *recorrido(int m,int a, int b, int x0, int x1, int y0, int y1, double *V);

int main(int argc, char** argv)
{
  
        MPI_Init(NULL,NULL);
	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int d=256/size;
        int  x0, x1, y0, y1, i=1, j=1, n=0,a,b;
	double average;
	
	m = L/h;
	N = 2*m*m;

	x0 = m/2 - l/(h*2) - 1;
	x1 = m/2 + l/(h*2) - 1;
	y0 = m/2 - d/(h*2) - 1;
	y1 = m/2 + d/(h*2) - 1;
	
	double *V = malloc(m*m*sizeof(double));
	double *V_new = malloc(m*m*sizeof(double));
	double *out_array1 = malloc(m*sizeof(double));
	double *out_array2 = malloc(m*sizeof(double));
	double *in_array1 = malloc(m*sizeof(double));
	double *in_array2 = malloc(m*sizeof(double));
	V = init(x0, x1, y0, y1, V);

	while (n < N)
	{	
	  if(rank==0)
	    {
	      a=1;
	      b=d;
	      V_new = recorrido(m, a, b, x0, x1, y0, y1, V);
	      for(i=1;i < m-1; i++)
		  {
		    out_array1[i] = V_new[transformer(i,b-1)];
		    for(j=a+1;j < b-1; j++)
			  {
			    V[transformer(i,j)] = V_new[transformer(i,j)] ;
			  }
		  }
	      MPI_Recv(&in_array1, m,MPI_DOUBLE, 1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      MPI_Send(&out_array1,m,MPI_DOUBLE, 1, 0, MPI,COMM_WORLD);

	      for(k=1;k<m-1;k++)
		{
		  V[transformer(k,b)]=in_array1[k];
		}
	      
	    }else if(rank==size-1)
	    {
	      a=(size-1)*d-1;
	      b=m-1;
	      V_new = recorrido(m, a, b, x0, x1, y0, y1, V);
	      for(i=1;i < m-1; i++)
		  {
		    out_array1[i] = V_new[transformer(i,a+1)];
		    for(j=a+1;j < b-1; j++)
			  {
			    V[transformer(i,j)] = V_new[transformer(i,j)] ;
			  }
		  }
	      MPI_Recv(&in_array1, m,MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      MPI_Send(&out_array1,m,MPI_DOUBLE, rank-1, 0, MPI,COMM_WORLD);

	      for(k=1;k<m-1;k++)
		{
		  V[transformer(k,a)]=in_array1[k];
		}

	    }else
	    {
	      a=rank*d -1;
	      b=rank*d;
	      V_new = recorrido(m, a, b, x0, x1, y0, y1, V);
	      for(i=1;i < m-1; i++)
		  {
		    out_array1[i]=V_new[transformer(i,a+1)];
		    out_array2[i]=V_new[transformer(i,b-1)]
		    for(j=a+1;j < b-1; j++)
			  {
			    V[transformer(i,j)] = V_new[transformer(i,j)] ;
			  }
		  }
	      MPI_Recv(&in_array1, m,MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      MPI_Recv(&in_array2, m,MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      
	      MPI_Send(&out_array1, m,MPI_DOUBLE, rank-1, 0, MPI,COMM_WORLD);
	      MPI_Send(&out_array2, m,MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);

	      for(k=1;k<m-1;k++)
		{
		  V[transformer(k,a)]=in_array1[k];
		  V[transformer(k,b)]=in_array2[k];
		}
	      
	    }
	}

	for(i=0;i < size; i++)
	{
	  if(rank==i)
	    {
	      if(rank==0)
		{
		  a=0;
		  b=d;
		}
	      else if(rank==size-1)
		{
		  a=(size-1)*d -1;
		  b=m;
		}
	      else
		{
		  a=rank*d -1;
		  b=(rank+1)*d;
		}
	      
	      for(j=a;j < b; j++)
		{
		  for(k=0;k<m;k++)
		    {
			printf("%f\n", V[transformer(j,k)]);
		    }
		}
	    }
		
	}
return 0;
}




int transformer(int i, int j)
{	
	
	return i*m + j;
}

double *init(int x0, int x1, int y0, int y1, double *array)
{	
	int a;
	for(a = x0; a <= x1; a++)
	{
		array[transformer(y0, a)] = V0/2;
		array[transformer(y1, a)] = -V0/2;
	}
	return array;
}

double *recorrido(int m,int a, int b, int x0, int x1, int y0, int y1, double *V)
{
  int up, down, left ,right,i,j;
  double V_new = malloc(m*m*sizeof(double));
  for(i=1;i <m-1 ; i++)
    {
      for(j=a+1;j < b-1; j++)
	{
	  up = transformer(i-1, j);
	  down = transformer(i+1, j);
	  left = transformer(i, j-1);
	  right = transformer(i, j+1);

	  if (!(j >= x0 && j <= x1 && i == y0) && !(j >= x0 && j <= x1 && i == y1))
	    {
	      average = (V[up] + V[down] + V[left] + V[right])/4;
	      V_new[transformer(i,j)] = average;
	    }
	}
    }

  return V_new;
}
