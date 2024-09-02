#include <stdio.h>
#include<stdlib.h>
#include<math.h>
//#include<conio.h>
int main()
{
    int i,j,m, n, iteration = 1;
    double del_x, del_y, error=1.0, errori,errorj;
    double l, h, b; // L=length along x, H=length along y;
    m=31; // no. of column of matrix designated by i , value of i decides the column no. or the verticle line in the grid 
    n=21; //no. of row of the matrix designated by j, value of j decides the row no. or the horizontal line in the grid    del_x = 0.2;
    del_y = 0.2;

    b = 1;

    ////////////////////////tdma/////////////////

    double ai, bi, ci, di[m-2], pi[m-2], qi[m-2];
    double aj, bj, cj, dj[n-2], pj[n-2], qj[n-2];
    double psi_old[m][n], psi_new[m][n],psi;
    
    ai=-2.0*(1+b*b);
    bi=1.0;
    ci=1.0;
    
    aj=-2.0*(1+b*b);
    bj=1.0;
    cj=1.0;


    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            psi_new[i][j] = 0.0; // Interier points
        }
    }

    for (i = 6; i < m; i++)
    {
        psi_new[i][0] = 100.0; // bottom boundary after the inlet to the m
    }

    FILE *file1,*fp;
    file1=fopen("P1_iteratiom_vs_error_ADI.txt", "w");
     fp=fopen("log_iteration_vs_log_error.plt","w");
        fprintf(fp, "Variables= \"X\",\"Y\"\n");
        fprintf(fp, "zone T= \"BLOCK1\", I=10, J=10\n\n");
    do
    {
    	
    	  for(j=0;j<n;j++)
        {
            for(i=0;i<m;i++)
            {

                psi_old[i][j] = psi_new[i][j];
            }
            
        }
    error = 0.0;
    errori=0.0;
    errorj=0.0;
    
    
        for( j = 1; j < n-1; j ++){
        di[0] = -(b*b)*(psi_new[0][j+1] + psi_new[0][j-1]);
        pi[0] = -bi/ai;
        qi[0] = di[0]/ai;

        for( i=1; i < m-1; i++){
            di[i] = -(b*b)*(psi_new[i][j+1] + psi_new[i][j-1]);;
            pi[i] = -(bi/(ai + ci*pi[i-1]));
            qi[i] = (di[i] - ci*qi[i-1])/(ai + ci*pi[i-1]);
        }
////////////// back substitution
        for( i=m-2; i>0; i--){
            psi = psi_new[i][j];
                psi_new[i][j] = pi[i]*psi_new[i+1][j]+qi[i];
              errori+=pow( (psi_new[i][j] - psi),2.0 );
        }

    }
     for (j = 0; j < n; j++){
        psi_new[m-1][j] = psi_new[m-2][j];
    }

    
     for( i = 1; i < m-1; i ++){
        dj[0] = -(b*b)*(psi_new[i+1][0] + psi_new[i-1][0]);
        pj[0] = -bj/aj;
        qj[0] = dj[0]/aj;

        for( j=1; j < n-1; j++){
            dj[j] = -(b*b)*(psi_new[i+1][j] + psi_new[i-1][j]);;
            pj[j] = -(bj/(aj + cj*pj[j-1]));
            qj[j] = (dj[j] - cj*qj[j-1])/(aj + cj*pj[j-1]);
        }

        for( j=n-2; j>0; j--){
            psi = psi_new[i][j];
                psi_new[i][j] = pj[j]*psi_new[i][j+1]+qj[j];
              errorj+=pow( (psi_new[i][j] - psi),2.0 );
        }

    }

    for (j = 0; j < n; j++){
        psi_new[m-1][j] = psi_new[m-2][j];
    }

    errori = sqrt( errori/((m-2)*(n-2)) );
    errorj = sqrt(errorj/((m-2)*(n-2)));
    
    for(i=1;i<m-1;i++)
	{
		for(j=1;j<n-1;j++)
		{
			error+=pow( (psi_new[i][j] - psi_old[i][j]),2.0 );
		}
	}
     
     
    printf("Iteration = %d \t Errori= %lf \t Errorj= %lf \t Error= %lf  \n", iteration, errori,errorj, error);
    fprintf(file1,"Iteration = %d \t Error= %lf \n", iteration, error);
    fprintf(fp, "%lf\t %.10f\n",log10(iteration),log10(error));
    iteration++;

    } while (error > 1e-6);

    fclose(fp);
    
    printf("\nlast iteration\n");
    for(j=n-1;j>=0;j--)
        {
            for(i=0;i<m;i++)
            {

                printf(" %lf",psi_new[i][j]); // printing final values
            }
            printf("\n");
        }
        
         FILE *file2;
        file2=fopen("stream_ADI.plt","w");
        fprintf(file2, "Variables= \"X\",\"Y\",\"PHI\"\n");
        fprintf(file2, "zone T= \"BLOCK1\", I=31, J=21, F=Point\n\n");
        
        for(j=0;j<n;j++)
        {
            for(i=0;i<m;i++)
            {

                fprintf(file2, "%lf \t %lf \t %lf \n", (i)*0.2, (j)*0.2, psi_new[i][j]);
            }
            
        }
        

    
   // getch();
    

    return 0;
}

