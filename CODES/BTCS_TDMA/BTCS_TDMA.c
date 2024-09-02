#include <stdio.h>
#include <stdlib.h>
#include<math.h>


int main()
{
	int m = 100; //Number of grid points 
	int i,j;
	//define the grid size
	
	double t, dt= 0.01;
	
	//Define two matrices to store the psi values
	double u_old[m], u_new[m];
	
	double temp, gamma,error=1;
	
	gamma =   1;
	double ai, bi, ci, d[m-2], p[m-2], q[m-2];
	
	ai=3; bi=-1; ci=-1;
	
	 FILE  *fp1, *file;
    fp1=fopen("log_iteration_vs_log_error.plt","w");
    file=fopen("error1_vs_time1.dat","w");
        fprintf(fp1, "Variables= \"X\",\"Y\"\n");
        fprintf(fp1, "zone T= \"BLOCK1\", I=10, J=1\n\n");
	
	
	
	for(i = 0; i < m; i++)
	{
			if( i == (m-1) )
			{
				u_new[i] = 1.0; //Top boundary
			}
		
			else if ( i == 0)
			{
				u_new[i] = 0.0; //bottom boundary
			}
			
			else
			{
				 u_new[i] = 0.0; //interior points
			}
			
		
	}
	
		for(i = 0; i < m; i++ )
	{
		
		
			printf(" %lf ", u_new[i]);
		
	
	}
	printf("\n");
		
	
	//Jacobi Iterative Method
	int iteration = 0;

	char name[50];

	
	while(error>1e-6)
	{
		if((iteration+100)%100==0)
		{
			sprintf(name,"velocity_%.2f.dat",t);
			FILE *fp;
			fp=fopen(name,"w");
			for(i=0;i<m;i++)
			{
				fprintf(fp, "%lf\t%lf\n",u_new[i],i*.01);
			}
			fclose(fp);
	
		}
		
			for(i = 0; i < m; i++)
		{
			
			
				u_old[i]= u_new[i]; //replace old psi with new psi values
				
			
		}
		
		
		
		d[0]=u_new[0];
		p[0]=1/3;
		q[0]=d[0]/3;
		//p[m-1]=0;
		//q[m-1]=u_new[m-1];
		
		for(i=1;i<m-1;i++)
		{
			d[i]=u_new[i];
			p[i]=1/(3-p[i-1]);
			q[i]=(d[i]+q[i-1])/(3-p[i-1]);
			
			
		}
		error=0.0;
		
		for(i = m-2 ; i > 0; i--)
		{
			
			
				//temp_new[i][j]=((1.0/ap) * (-as*temp_new[i][j-1]-aw*temp_new[i-1][j]-ae*temp_old[i+1][j]-an*temp_old[i][j+1]));
				//temp_new[i][j]=0.25*(temp_old[i-1][j]+temp_old[i+1][j]+temp_old[i][j-1]+temp_old[i][j+1]);
				//u_new[i]= u_old[i]+gamma*(u_old[i+1]-2.0*u_old[i]+u_old[i-1]);
				temp=u_new[i];
				u_new[i]=p[i]*u_new[i+1]+q[i];
					error=error+ pow(u_new[i]-temp,2);
				
							
		}
					error=sqrt(error/(m));
		iteration++;
		fprintf(fp1,"%lf\t %lf\n",t,error);
		fprintf(file,"%lf\t %lf\n",t, error);
		t=t+dt;
		
		printf("time=%lf\t error=%lf\n",t, error);
		
	
	
	
		
		iteration++;
	
	}
	

	
		for(i = 0; i < m; i++ )
	{
		
		printf("%lf\n",u_new[i]);
	}
	
	
       printf("\n %d\t time=%lf",iteration,t); 
	
	
	return 0;
}

