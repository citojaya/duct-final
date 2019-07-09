#include "common.h"

/*--- Find the reading position in the file--*/
void findRec(FILE *inFile, char* strDest){
	int nbytes = 256;
	char* strSrc;
	strSrc = (char *)malloc(nbytes+1);

	rewind(inFile);
	int n=strlen(strDest);
	while(!feof(inFile)){
		fgets(strSrc, 256, inFile);
		strSrc[n]='\0';
		if (strcmp(strDest, strSrc) == 0){
			break;
		}
	}

	if(strcmp(strDest, strSrc) != 0){
		exit(1);
	}
	free(strSrc);
}

/*---Read input data from a file ----*/
void readInput(char *infile, int *np, real *dens, real *ymod, 
			real *pois, real *sfc, real *rec, real *dmpn, real *rf, real *cyldia,
			 real *dt, int *nW, int *updateDPM, real *maxVel){
	// input file reading
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *InFile = fopen(filename, "rt");

	if (InFile == NULL){
		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
		char c = getchar();
		exit(1);
	}

	// expand size for computing area
	real exComp = 0.0;

	// loading phase
	real DieDepth = 1.0;
	real UnDepth  = 0.0;
	real BdDepth  = 0.0;

	real parDia = 0.0;
	
	findRec(InFile, "PAR_NUMBER");
	fscanf(InFile, "%d",  np);

	findRec(InFile, "MATERIAL");
	fscanf(InFile, "%lf", dens);
	fscanf(InFile, "%lf", ymod);
	fscanf(InFile, "%lf", pois);
	fscanf(InFile, "%lf", sfc);
	fscanf(InFile, "%lf", dmpn);
	fscanf(InFile, "%lf", rf);
	fscanf(InFile, "%lf", &haPP);
	fscanf(InFile, "%lf", &haPW);

	findRec(InFile, "SIMULATION");
	fscanf(InFile, "%lf", dt);	

	findRec(InFile, "WALLS");
	fscanf(InFile, "%d", nW);	

	findRec(InFile, "DPM");
	fscanf(InFile, "%d", updateDPM);	

	findRec(InFile, "MAXFLOWVEL");
	fscanf(InFile, "%lf", maxVel);	

	findRec(InFile, "PERMITIVITY");
	fscanf(InFile, "%lf", &permitivity);
	findRec(InFile, "CAPACITANCEDISTANCE");
	fscanf(InFile, "%lf", &Zs);
	findRec(InFile, "VOLTAGE1");
	fscanf(InFile, "%lf", &V1);
	findRec(InFile, "VOLTAGE2");
	fscanf(InFile, "%lf", &V2);
	findRec(InFile, "IMAGECONSTANT");
	fscanf(InFile, "%lf", &imageConst);
	findRec(InFile, "ALPHA");
	fscanf(InFile, "%lf", &alpha);
	findRec(InFile, "KS");
	fscanf(InFile, "%lf", &ks);
	findRec(InFile, "ESFTRUE");
	fscanf(InFile, "%lf", &esfTrue);
	

	findRec(InFile, "ROUGHSURFACE");
	fscanf(InFile, "%lf", &lamda1);
	fscanf(InFile, "%lf", &lamda2);
	fscanf(InFile, "%lf", &rms1);
	fscanf(InFile, "%lf", &rms2);

	findRec(InFile, "CAPILLARY");
	fscanf(InFile, "%lf", &s_min);
	fscanf(InFile, "%lf", &liq_vol);
	fscanf(InFile, "%lf", &surf_tens);
	fscanf(InFile, "%lf", &cont_ang);
	
    fclose(InFile);
	//fclose(LogFile);
}

void writeLogNum(char *infile, char *line, real num){
	FILE *LogFile = fopen(infile, "a");
	fprintf(LogFile,line);
	fprintf(LogFile,"%lf\n",num);
	fclose(LogFile);
}

void writeLog3Num(char *infile, char *line, real v1, real v2, real v3){
	FILE *LogFile = fopen(infile, "a");
	fprintf(LogFile,line);
	fprintf(LogFile,"%lf %lf %lf\n",v1,v2,v3);
	fclose(LogFile);
}

void writeLogLine(char *infile, char *line){
	FILE *LogFile = fopen(infile, "a");
	fprintf(LogFile,line);
	fclose(LogFile);
}

void writeInjectionFile(char *infile){
	char filename[20];
	strcpy(filename, infile);
	FILE *injFile = fopen(filename, "w");
	Injection *I;
    Injection *Ilist = Get_dpm_injections();

    loop(I,Ilist)
    {
    	Particle *p;
        loop(p,I->p)
        { 
			char line[100];
			strcpy(line, "((");

			char pX[8];
			sprintf(pX, "%f", demPart[p->part_id].pos[0]/lengthFactor);
			strcat(line ,pX); 
			strcat(line ," "); 
			char pY[8];
			sprintf(pY, "%f", demPart[p->part_id].pos[1]/lengthFactor);
			strcat(line ,pY); 
			strcat(line ," ");
			char pZ[8];
			sprintf(pZ, "%f", demPart[p->part_id].pos[2]/lengthFactor);
			strcat(line ,pZ);
			strcat(line ," "); 
			char pVelX[8];
			sprintf(pVelX, "%f", demPart[p->part_id].vel[0]/velocityFactor);
			strcat(line ,pVelX);
			strcat(line ," ");
			char pVelY[8];
			sprintf(pVelY, "%f", demPart[p->part_id].vel[1]/velocityFactor);
			strcat(line ,pVelY);
			strcat(line ," ");
			char pVelZ[8];
			sprintf(pVelZ, "%f", demPart[p->part_id].vel[2]/velocityFactor);
			strcat(line ,pVelZ);
			strcat(line ," ");
			char pDia[8];
			sprintf(pDia, "%f", demPart[p->part_id].dia/lengthFactor);
			strcat(line ,pDia);
			strcat(line ," 0.0 1.0))\n");
			fprintf(injFile, line);
		}
	}

	fclose(injFile);
}

void diaInput(char *infile, struct demParticle *par, int *np){
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *pDiaFile = fopen(filename, "rt");

	if (pDiaFile == NULL){
		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
		char c = getchar();
		exit(1);
	}
	int num = 0;

	real pDia, pX, pY, pZ;
    //printf("No of par %d\n",np);
	findRec(pDiaFile, "PARTICLE");
	for(int i=0; i<*np; i++){
		fscanf(pDiaFile, "%lf", &pDia);
		fscanf(pDiaFile, "%lf", &pX);
		fscanf(pDiaFile, "%lf", &pY);
		fscanf(pDiaFile, "%lf", &pZ);
		demPart[i].dia = pDia;
		demPart[i].pos[0] = pX;
		demPart[i].pos[1] = pY;
		demPart[i].pos[2] = pZ;
		
	}
	fclose(pDiaFile);
}

void readWalls(char *infile, int *walls){
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *wFile = fopen(filename, "rt");

	if (wFile == NULL){
		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
		char c = getchar();
		exit(1);
	}

	findRec(wFile, "WALL_NO");
	for(int i=0; i<noOfWalls; i++){
		int wallNo;
		fscanf(wFile, "%d", &wallNo);
		walls[i] = wallNo;
	}
	fclose(wFile);
}

real readInputVelocity(char *infile){
	char filename[20];
	//real inletVel;
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *f = fopen(filename, "rt");
	findRec(f, "MAXFLOWVEL");
	fscanf(f, "%lf", &inletVel);

	fclose(f);
	return inletVel;
}

void readGeom(char *infile, real *ductxmin, real *ductxmax, real *ductxedge1, real *ductxedge2, real *ductymin, 
            real *ductymax, real *ductzmin, real *ductzmax, real *ductzedge){
	char filename[20];
	//real inletVel;
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *f = fopen(filename, "rt");
	findRec(f, "GEOMETRY");
	
	fscanf(f, "%lf", ductxmin);
	fscanf(f, "%lf", ductxmax);
	fscanf(f, "%lf", ductxedge1);
	fscanf(f, "%lf", ductxedge2);
	fscanf(f, "%lf", ductymin);
	fscanf(f, "%lf", ductymax);
	fscanf(f, "%lf", ductzmin);
	fscanf(f, "%lf", ductzmax);
	fscanf(f, "%lf", ductzedge);

	fclose(f);	
}


void readDomain(char *infile){
	char filename[20];
	//real inletVel;
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *f = fopen(filename, "rt");
	findRec(f, "BOUNDARY");
	
	fscanf(f, "%lf", &xmin);
	fscanf(f, "%lf", &xmax);
	fscanf(f, "%lf", &ymin);
	fscanf(f, "%lf", &ymax);
	fscanf(f, "%lf", &zmin);
	fscanf(f, "%lf", &zmax);

	findRec(f, "REFERENCEVALUES");
	fscanf(f, "%lf", &largestParDia);
	fscanf(f, "%lf", &largestParDensity);

	fclose(f);	
}

void printCPUTime(){
    clock_t CPU_time_1 = clock();
    writeLogNum("logfile2.log","DEM TIME ",demTime/timeFactor);
    writeLogNum("logfile2.log","CPU TIME ",CPU_time_1/1e3);
   
    prevCPUTime = CPU_time_1;
}

void demSave(){
	//FILE *outfile; 
	char filename[20];
	sprintf(filename, "particle.dat");
	FILE *outfile = fopen(filename, "a");
	fprintf(outfile, "TIME = %lf\n",CURRENT_TIME);

  	Injection *I;
  	Injection *Ilist = Get_dpm_injections();
  
  	// Update FLUENT particle postion and velocity 
  	int ip = 0; 
	 loop(I,Ilist)
  	 {
     	Particle *p;
     	loop(p,I->p)
     	{
			if(demPart[p->part_id].active == 1){
		 	fprintf(outfile, "%11.5lf   %11.5lf   %11.5lf   %11.5f  %11.5lf   %11.5lf  %11.5lf %11.5lf %11.5lf %d %d %11.5lf %lld\n",
			demPart[p->part_id].pos[0]/(lengthFactor*conversion), demPart[p->part_id].pos[1]/(lengthFactor*conversion),
			demPart[p->part_id].pos[2]/(lengthFactor*conversion),
			demPart[p->part_id].vel[0]/(velocityFactor),demPart[p->part_id].vel[1]/(velocityFactor),demPart[p->part_id].vel[2]/(velocityFactor),
			demPart[p->part_id].dia/(lengthFactor*conversion),
			demPart[p->part_id].maxWallCollE*1.0e12/energyFactor,
			demPart[p->part_id].maxPartCollE*1.0e12/energyFactor,
			demPart[p->part_id].noOfWallColl,
			demPart[p->part_id].noOfPartColl,
			demPart[p->part_id].eCharge*1.0e15,
			p->part_id);
			}
		 }			
     } 
	fclose(outfile);	
}

void cordSave(){

}

void  writeFluidVelocity(){
	char filename[20];
	sprintf(filename, "cfdvelocity.dat");
	FILE *outfile = fopen(filename, "a");
	fprintf(outfile, "TIME = %lf\n",CURRENT_TIME);

	Domain *d = Get_Domain(1);
	Thread *t;
	cell_t c;
	real x[ND_ND];
	thread_loop_c (t,d)
	begin_c_loop(c,t)
		//int cI = C_UDMI(c,t,0); //initialize cfdcell solid volume
		C_CENTROID(x,c,t);
		fprintf(outfile, "%11.5lf   %11.5lf   %11.5lf   %11.5f  %11.5lf   %11.5lf \n",
			x[0]/conversion,x[1]/conversion,x[2]/conversion,x[2]/conversion,
			C_U(c,t),C_V(c,t),C_W(c,t));
	end_c_loop(c,t)
	fclose(outfile);
}

