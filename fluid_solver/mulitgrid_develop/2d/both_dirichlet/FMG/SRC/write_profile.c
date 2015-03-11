void write(double *y)
{
	//==============================================================================//
	char fname[300];
	FILE *fp;
	//==============================================================================//


	sprintf(fname,"../DATA/y.dat");
	if((fp=fopen(fname,"w"))==NULL)
	{
		printf("y profile can't be written\n");
		exit(1);
	}

	long i,j;
	
	long array_index;
	
	for(j=0;j<nodes[0].y;j++)
	{				
		for(i=0;i<nodes[0].x;i++)
		{
			array_index=i+nodes[0].x*j; //this is not general but will work for the finest grid	
	
			fprintf(fp,"%lf\t%lf\t%lf\n",i*d_sp[0].x,j*d_sp[0].y,y[array_index]);
		}
		fprintf(fp,"\n");
	}	
			
	fclose(fp);
}
