void write(double *y,int iter)
{
	//==============================================================================//
	char fname[300];
	FILE *fp;
	//==============================================================================//


	sprintf(fname,"../DATA/y_at_iter_%d.dat",iter);
	if((fp=fopen(fname,"w"))==NULL)
	{
		printf("y profile can't be written\n");
		exit(1);
	}

	int i,j;
	
	long array_index;	

	for(j=0;j<nodes.y;j++)
	{				
		for(i=0;i<nodes.x;i++)
		{
			array_index=i+nodes.x*j;
		
			fprintf(fp,"%lf\t%lf\t%lf\n",i*d_sp.x,j*d_sp.y,y[array_index]);
		}

		fprintf(fp,"\n");	
	}	
			
	fclose(fp);
}
