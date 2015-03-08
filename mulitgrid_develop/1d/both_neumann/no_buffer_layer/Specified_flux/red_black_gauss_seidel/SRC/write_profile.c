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

	int i;
					
	for(i=0;i<nodes_x;i++)
	{
		fprintf(fp,"%lf\t%lf\n",i*dx,y[i]);
	}	
			
	fclose(fp);
}
