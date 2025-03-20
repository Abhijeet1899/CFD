#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
	int i,j; /// initializing grid indices || here i movers in x direction (i.e represents no. of colums of grid matrix) ||  here j moves in y direction (i.e represents no. of rows of the grid matrix
	int m,n; // m-x(i) limit, n-y(j) limit
	printf("enter the total no. of grid point (M(x),N(y))=");
	scanf("%d %d",&m,&n);
	//m=100;n=100;
	double U=1.0; // velocity of the lid in x direction
	printf("M=%d\tN=%d\n",m,n);
	
	double h,l;// h- height of the lid (y)|| l- lenght of the lid(x)
	//printf("enter the total no. of grid point (L(x),H(y))=");
	//scanf("%lf %lf",&l,&h);
	l=1.0;h=1.0;
	printf("L=%lf\tH=%lfU=%lf\n",l,h,U); 
	
	double Re;
	printf("enter the Re value =");
	scanf(" %lf",&Re);
	//Re=400.0;
	
	double w[m][n], w_old[m][n],psi[m][n], psi_old[m][n],u[m][n],v[m][n]; ///initialising worticity, x-velocity, y-velocity, stream function
	// the first indices of the array defines the row that is y 
	// the second indices of the array defines the column that is x
	
	double dx=l/m, dy=h/n, beta=dx/dy ;
	printf("dx=%lf\tdy=%lf\n\tbeta=%lf\n",dx,dy,beta);
	
	int iteration=0;
	double error_psi, error_w;
	double p=(m-2)*(n-2);
	// initiallizing psi values , taking psi = 0.0 initially
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			psi[i][j]=0.0;
			w[i][j]=0.0;
		//	printf("psi=%lf\t w=%lf ",psi[i][j],w[i][j]);
		}
	//	printf("\n");
	}
	

	// velocity and stream function b.c
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			if(i==0)
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
			// left B.C 
				
			}
			else if(j==0)
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
			//  bottom B.C
			
			}
			else if(i==m-1)
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
			// right B.C
				
			}
			else if(j==n-1)
			{
				u[i][j]=1.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
				 // Top B.C 
				
			}
			
			
			else
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
			}
			
		}
		
	
	}
	
	/// omega (w) b.c
		for(i=0;i<m;i++)
		{
			for(j=0;j<n;j++)
			{
				if(i==0)
				{
					w[i][j]=-2.0*(psi[1][j]-psi[0][j])/(pow(dx,2)); // left B.C 
					
				}
				else if(j==0)
				{
					w[i][j]=-2.0*(psi[i][1]-psi[i][0])/(pow(dy,2)); //  bottom B.C
			
				}
				else if(i==m-1)
				{
					w[i][j]=-2.0*(psi[m-2][j]-psi[m-1][j])/(pow(dx,2)); // right B.C
				
				}
				else if(j==n-1)
				{
					w[i][j]=-2.0*(psi[i][n-2]-psi[i][n-1]+U*dy)/(pow(dy,2)); // Top B.C 
				
				}
				else
				{
					w[i][j]=0.0;
				}
			
			}		
		}
	
//	
	do
	{
		for(i=0;i<m;i++)
		{
			for(j=0;j<n;j++)
			{
				psi_old[i][j]=psi[i][j];
				w_old[i][j]=w[i][j];
			}
		}
		
		for(i=1;i<m-1;i++)
		{
			for(j=1;j<n-1;j++)
			{
				psi[i][j]= (0.5/(1.0+pow(beta,2)))*(dx*dx*w[i][j]+(pow(beta,2))*(psi[i][j+1]+psi[i][j-1])+psi[i+1][j]+psi[i-1][j]);
			
			}
		}
		
		for(i=1;i<m-1;i++)
		{
			for(j=1;j<n-1;j++)
			{
				u[i][j]=(psi[i][j+1]-psi[i][j-1])/(2.0*dy);
				v[i][j]=-(psi[i+1][j]-psi[i-1][j])/(2.0*dx);
			
			}
		}
		
		for(i=1;i<m-1;i++)
		{
			for(j=1;j<n-1;j++)
			{
				w[i][j]=(0.5/(1+pow(beta,2)))*((1.0-u[i][j]*dx*Re/2.0)*w[i+1][j]+(1.0+u[i][j]*dx*Re/2.0)*w[i-1][j]
				+ (1.0-v[i][j]*dy*Re/2.0)*w[i][j+1]*pow(beta,2)+(1.0+v[i][j]*dy*Re/2.0)*w[i][j-1]*pow(beta,2));
			
			}
		}
		
		for(i=0;i<m;i++)
		{
			for(j=0;j<n;j++)
			{
				if(i==0)
				{
					w[i][j]=-2.0*(psi[1][j]-psi[0][j])/(pow(dx,2)); // left B.C 
					
				}
				else if(j==0)
				{
					w[i][j]=-2.0*(psi[i][1]-psi[i][0])/(pow(dy,2)); //  bottom B.C
			
				}
				else if(i==m-1)
				{
					w[i][j]=-2.0*(psi[m-2][j]-psi[m-1][j])/(pow(dx,2)); // right B.C
				
				}
				else if(j==n-1)
				{
					w[i][j]=-2.0*(psi[i][n-2]-psi[i][n-1]+U*dy)/(pow(dy,2)); // Top B.C 
				
				}
			
			}		
		}
		
		error_psi=0.0;
		error_w=0.0;
		printf("error_psi=%lf\terror_w=%lf\t",error_psi,error_w);
		for(i=0;i<m;i++)
		{
			for(j=0;j<n;j++)
			{
				error_psi=error_psi+pow(psi[i][j]-psi_old[i][j],2);
				error_w=error_w+pow(w[i][j]-w_old[i][j],2);
			}
		}
		error_psi=sqrt(error_psi/p);
		error_w=sqrt(error_w/p);
		
		printf("iteration=%d\terror_psi=%lf\terror_w=%lf\n",iteration,error_psi,error_w);
		iteration++;
		

	}while(error_psi>1.0e-6 || error_w>1.0e-6);
	
	
	for(i=0;i<m-1;i++)
	{
		for(j=0;j<n-1;j++)
		{
			//psi[i][j]= dx*dx*w[i][j]+psi[i][j+1]+psi[i][j-1]+psi[i+1][j]+psi[i-1][j];
			//psi[i][j]=(1/(1+pow(beta,2)))*(pow(dx,2)*w[i][j]+(pow(beta,2))*(psi[i][j+1]+psi[i][j-1])+psi[i+1][j]+psi[i-1][j]);
		}
		//printf("\n");
	}
	
	double x,y;
	FILE *fp;
	if(Re==100)
	{
		fp=fopen("output_Re100.dat","w");
	}
	else if(Re==400){
		fp=fopen("output_Re400.dat","w");
	}
	fprintf(fp,"ZONE I=%d, J=%d\n",m,n);
	for(j=0;j<n;j++)
	{
		y=dy*j;
		for(i=0;i<m;i++)
		{
			x=dx*i;
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,y,u[i][j],v[i][j],psi[i][j],w[i][j]);
		}
		
	}
	fclose(fp);
	
	FILE *f;
	f=fopen("u_mid1.dat","w");
	//fprintf(f,"ZONE I=%d, J=%d\n",m,n);
	
		for(i=0;i<m;i++)
		{
			x=dx*i;
			
			
			
				fprintf(f,"%lf\t%lf\n",u[i][50],x);	
			
		
	}
	fclose(f);
	
	

	
}

