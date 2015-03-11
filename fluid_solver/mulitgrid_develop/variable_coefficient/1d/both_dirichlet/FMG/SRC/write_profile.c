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

	int i;
					
	for(i=0;i<nodes_x[0];i++)
	{
		fprintf(fp,"%lf\t%lf\n",i*dx[0],y[i]);
	}	
			
	fclose(fp);
}
