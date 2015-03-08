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
					
	for(i=1;i<nodes_x-1;i++)
	{
		fprintf(fp,"%lf\t%lf\n",(i-1)*dx,y[i]);
	}	
			
	fclose(fp);
}
