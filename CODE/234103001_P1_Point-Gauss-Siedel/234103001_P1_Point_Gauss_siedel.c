#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include<conio.h>
int main()
{
    int row,col,k,l;// initializing variables
  //  printf("enter the grid size J , I\n");
   // scanf("%d  %d", &row, &col); // taking grid size input
    //printf("%d %d\n",row,col);
    int i,j;
   	row=21;col=31; // no. of grid initializing
    double del=0.2;
    int c=row,b=col;
	
    double a[row+1][col+1]; //initialising grid points 
    double a1[row+1][col+1]; // initializing new grid points
    double d,e;
    k=0;
    e=1;
    
    FILE *file1,*fp;
    file1= fopen("P1_iteratiom_vs_error_PGS.txt","w");//write the errors in a file
    fp=fopen("log_iteration_vs_log_error.plt","w");
        fprintf(fp, "Variables= \"X\",\"Y\"\n");
        fprintf(fp, "zone T= \"BLOCK1\", I=10, J=10\n\n");
    do
    {
    	k=k+1;
       ;
       d=0;
        for(j=1;j<=row;j=j+1)
        {
            for(i=1;i<=col;i=i+1)
            {
            	if(k==1)
                {
                    a[j][i]=0; // intial value for interior points
                     
                }
				
               
                

                else
                {
					if((i>1 && j>1) && (j<row && i<col) )
					{
					
                		a[j][i]=(a[j-1][i]+a[j][i-1]+a1[j][i+1]+a1[j+1][i])/4; //here beta=1 so the discritisation equation becomes the following
                		
                    	
                	}
            	}
                
                	a[j][1]=0;                            //left size BC                           
            		a[j][b]=(4*a[j][b-1]-a[j][b-2])/3;    //right size BC     
            		if(i<=6)
        			{
            			a[1][i]=0;                        //bottom BC befor 1m
            			a[c][i]=0;						  // top BC 
        			}

        			else if (i>=7)
         			{
         				a[1][i]=100;					// bottom BC after 1.2 meter
       	   				a[c][i]=0;						//top BC 
       		 
		 			}
		 			
				
            
				if(k>=2 && (j>=2 && i>=2) && (j<=row-1 && i<=col-1))
				{
				
					d=d+pow(a[j][i]-a1[j][i],2);  // error finding Root mean square value
					
            	
				
				}
				
			
        	}
        
        }
         for(j=1;j<=row;j=j+1)
        {
            for(i=1;i<=col;i=i+1)
            {
            	a1[j][i]=a[j][i];
            }
        }
        if(k>1){
		e=sqrt(d/551); // error 
		fprintf(file1, "%d\t %.10f\n",k,e);	
		printf( "iteration=%d\t error=%.10f\n",k,e);
		fprintf(fp, "%lf\t %.10f\n",log10(k),log10(e));	
		}
    }while(e > 1e-8);


		printf("\nlast iteration value\n");
     for(j=row;j>=1;j--)
        {
            for(i=1;i<=col;i++)
            {

                printf(" %.2f",a[j][i]); // printing final values
            }
            printf("\n");
        }
        e=sqrt(d/551);
        printf("\n%.10f",e);
        printf("\n no. of iteration= %d",k);
        
        
        //writing the streamfunction in a file to plot the contours
        FILE *file2;
        file2=fopen("stream_PGS.plt","w");
        fprintf(file2, "Variables= \"X\",\"Y\",\"PHI\"\n");
        fprintf(file2, "zone T= \"BLOCK1\", I=31, J=21, F=Point\n\n");
        
        for(j=1;j<=row;j++)
        {
            for(i=1;i<=col;i++)
            {

                fprintf(file2, "%lf \t %lf \t %lf \n", (i-1)*0.2, (j-1)*0.2, a[j][i]);
            }
            
        }
        
     
        getch();
    return 0;

}

