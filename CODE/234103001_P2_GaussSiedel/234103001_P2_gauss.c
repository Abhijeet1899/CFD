#include <stdio.h>
#include <stdlib.h>
#include<math.h>


int main()
{
	int m = 11, n = 11; //Number of grid points 
	int i,j;
	//define the grid size
	double dx, dy;
	dx = 0.1;
	dy = 0.1;
	
	
	//Define two matrices to store the psi values
	double temp_old[m][n], temp_new[m][n];
	
	//Define the S, W, P, E, N points
	double as, aw, ap, ae, an;
	as = (1.0/pow(dy,2.0));
	aw = (1.0/pow(dx,2.0));
	ap = -2.0 * ((1.0/pow(dx,2.0)) + (1.0/pow(dy,2.0)));
	ae = (1.0/pow(dx,2.0));
	an = (1.0/pow(dy,2.0));
	
	printf("dx=%lf  dy=%lf ||| %lf  %lf  %lf  %lf  %lf\n",dx,dy,as,aw,ap,ae,an);
	for(i = 0; i < m; i++)
	{
		for(j = 0; j < n; j++)
		{
			if( i == (m-1) )
			{
				temp_new[i][j] = 0.0; //Top boundary
			}
			else if ( j == 0)
			{
				temp_new[i][j] = 1.0; //left boundary
			}
			else if ( i == 0)
			{
				temp_new[i][j] = 1.0; //bottom boundary
			}
			else if ( j == (n-1) )
			{
				temp_new[i][j] = 1.0;//right boundary
			}
			else
			{
				 temp_new[i][j] = 0.0; //interior points
			}
			
		}
	}
	
		for(i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++)
		{
			printf(" %lf ", temp_new[i][j]);
		}
		printf("\n");
	}
	
		
	
	//Jacobi Iterative Method
	int iteration = 1;
	double error = 1.0;
	
	FILE *file1, *fp;
	file1= fopen("P1_iteratiom_vs_error_PGS.txt","w");//write the errors in a file
    fp=fopen("log_iteration_vs_log_error.plt","w");// iteration vs error plot
        fprintf(fp, "Variables= \"X\",\"Y\"\n");
        fprintf(fp, "zone T= \"BLOCK1\", I=10, J=10\n\n");
	
	do
	{
	
		for(i = 0; i < m; i++)
		{
			for(j = 0; j < n; j++)
			{
				temp_old[i][j] = temp_new[i][j]; //replace old psi with new psi values
			}
		}
		
		for(i = 1; i < (n-1); i++)
		{
			for(j = 1; j< (m-1); j++)
			{
				temp_new[i][j]=((1.0/ap) * (-as*temp_new[i][j-1]-aw*temp_new[i-1][j]-ae*temp_old[i+1][j]-an*temp_old[i][j+1]));
				//temp_new[i][j]=0.25*(temp_old[i-1][j]+temp_old[i+1][j]+temp_old[i][j-1]+temp_old[i][j+1]);
			}
		}
		
	
	
		error = 0.0;
		printf("%lf \n",error);
		for(i = 0; i < m; i++)
		{
			for(j = 0; j < n; j++)
			{
				error = error + pow(temp_new[i][j]-temp_old[i][j],2.0);
				//printf("\n%lf\n",error);
			}
		}
		
		error = sqrt(error/(m*n));
		
		printf ("Iteration %d\t", iteration);
		printf("Error %lf\n", error);
		fprintf(file1, "%d\t%lf\n", iteration, error);
		fprintf(fp, "%lf\t%lf\n", log10(iteration), log10(error));
		iteration++;
		
	}while(error > 1e-8 );
	
	printf("\n");
	
		for(i = 0; i < m; i++ )
	{
		for( j = 0; j < n; j++)
		{
			printf(" %lf ", temp_new[i][j]);
		}
		printf("\n");
	}
	
	 FILE *file2;
        file2=fopen("Temprature.plt","w");
        fprintf(file2, "Variables= \"X\",\"Y\",\"PHI\"\n");
        fprintf(file2, "zone T= \"BLOCK1\", I=11, J=11, F=Point\n\n");
        
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {

                fprintf(file2, "%lf \t %lf \t %lf \n", (i)*0.1, (j)*0.1, temp_new[j][i]);
            }
            
        }
        
	
	
	return 0;
}
