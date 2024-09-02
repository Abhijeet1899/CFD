#include <stdio.h>
#include <stdlib.h>
#include<math.h>


int main()
{
	int m = 100; //Number of grid points 
	int i,j,l;
	double lhs, rhs[m], residual;

	double t, dt= 0.01;
	double temp;
	//Define two matrices to store the psi values
	double u_old[m], u_new[m];
	
	
	double gamma;
	double error=1.0;
	
	gamma =  1.0;
	
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
		
	

	int iteration = 0;

	char name[50];

	
	while(t<500)
	{
	//	if((iteration+100)%100==0)
	//	{
		//	sprintf(name,"velocity_%.2f.dat",t);
		//	FILE *fp;
		//	fp=fopen(name,"w");
		//	for(i=0;i<m;i++)
		//	{
			//	fprintf(fp, "%lf\t%lf\n",u_new[i],i*.01);
		//	}
		//	fclose(fp);
	
	//	}
		
		for(i = 0; i < m; i++)
		{
			
			
				u_old[i]= u_new[i]; //replace old psi with new psi values
				
			
		}
		

		
		for(i=0;i<m;i++)
		{
			
			if(i==m-2)
			{	
				rhs[i]=-u_new[i]-1.0;
			}
			else
			{
				rhs[i]=-u_new[i];
			}
		}
	
		for(j=0;j<m;j++)
		{
			lhs=0.0;
			if(j==0)
			{
				lhs=lhs+u_new[j+1];
			}
			
			else if(j==m-2)
			{
				lhs=lhs+u_new[j-1]-3.0*u_new[j];
		
			}
		
			else if (j==m-1)
			{
				lhs=lhs+u_new[j-1]-3.0*u_new[j];
			}
			else 
			{
				lhs=lhs+u_new[j-1]-(3.0)*u_new[j]+u_new[j+1];
			}
			
			
			
				residual=rhs[j]-lhs;
			//	printf("residual(%d)=%lf\n",j,residual);
				if(j>=1 && j<=m-2)
				{
					temp=u_new[j];
					u_new[j]=u_new[j]+residual/(-3.0);
				}
			
			error=error+pow(u_new[j]-u_old[j],2);

		}
		
		error=sqrt(error/m);
		fprintf(fp1,"%lf\t %lf\n",t,error);
		fprintf(file,"%lf\t %lf\n",t, error);	
		t=t+dt;
		iteration++;
	
	}
	printf("\nerror=%lf\n",error);

	
		for(i = 0; i < m; i++ )
	{
		
		printf("%lf\n",u_new[i]);
	}
	
	
      printf("\n %d\t time=%lf",iteration,t);  
	
	
	return 0;
}

