
#ifndef COMMON_H
#define COMMON_H

#include "udf.h"
#include "dpm.h"
//#include "dem.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/*--------------Common definitions---------------*/


//#define NUM_MAT 2    // two types material
#define PI 3.1415926f
#define gravity 9.8f //ms^-2
#define conversion 1.0e-3f //convert length values to meters
#define overlapLimit 0.0005f //ovelap cuttoff value

/*------------ Particle information -------------*/
#define XMAX_BOUND 0.014 //Maxmimum x bound for active particles
#define multif3 2.2//multification factor used in bounding box divisions
					//simulation time inversly prortional to value
#define DIM 3 // 3D problem
#define NBSIZE 40 //size of neighbourlist
#define CNTSIZE 15
#define WALLCNTSIZE 8 //Maximum number of walls, a particle can have contact with
#define NO_OF_FACES 20//number of faces contacting with particles
#define NO_OF_PARTICLES_IN_BDCELL 40

// err detection utility
#define FPRINTF(a) fprintf a

/*-----------------------------------------------*/
/*--------------Global declaration---------------*/ 
/*-----------------------------------------------*/

//FILE *LogFile;
real prevCPUTime;
int cfdcycles, demcycles, iter, noOfIteration;
int updateSourceTerm;
int particle_counter;
int noOfWalls, noOfCFDCells, noOfWallFaces;
int *walls;
int np, parArraySize; //number of particles, particle array size 
unsigned int initialized;

real inletVel;
real ductxmin, ductxmax, ductxedge1, ductxedge2, ductymin, ductymax, ductzmin, ductzmax, ductzedge;
real cutGap; //particle parameters
real cellRadius; //Radius defined by bounding box cell given by 0.5*sqrt(dx^2+dy*2)
real refLength,refDensity,lengthFactor,volumeFactor,massFactor,timeFactor,
	densityFactor, forceFactor, pressureFactor, StressFactor, energyFactor, momentFactor,
	powerFactor, velocityFactor, accFactor, angVelFactor, angAccFactor, freqFactor, inertiaFactor, largestParDia,largestParDensity;
real vGapMn; //minimum gap for vanderwal force 
real maxVel; //maximum flow velocity at the inlet 
real lamda1, lamda2, rms1, rms2;
real s_min,liq_vol, surf_tens, cont_ang; 
int xDiv, yDiv, zDiv; //Number of divisions in orthogonal domain cells
real domainDx, domainDy, domainDz; //Domain cell size
real xmax, ymax, zmax; //Max values of domain boundary 
real xmin, ymin, zmin; //Min values of domain boundary 

real Zs, V1, V2, permitivity, imageConst, alpha, ks, esfTrue; //electrostatic constants
real rIn, rOut, allowedDisp; //if particle displacement > allowedDisp -> update neighbourlist

real *uVec, *ipRVec, *jpRVec, *ijVec, *rotVel;
real *ipCntPntVel, *jpCntPntVel, *cntPntVel;
real dens, ymod, pois, sfc, rf, rec, dmpn, elasticMod, haPP, haPW; //particle material property
real dsmaxCff, dti, dd, dsmax, fdt; //used in contact force calculation
real cyldia; //cylinder diameter
real timeStep, demTime, maxTime;
struct BdBox *bdBox;
struct demParticle *demPart;
struct cfdCell *cfdcell;
//Tracked_Particle *dpmList[2];
int updateDPM;
int saveDEM; //counter used for saving DEM particle for TECPLOT
int maxCnt;
real maxCharge, maxES, maxVF;
unsigned int maxChargePart;

//cfdcell information
struct cfdCell{
	real porosity;
	real solidVol;
	real dragFX;
	real dragFY;
	real dragFZ;
	int noOfParts;
};

//Walls face
struct wallFace{
	real node1[DIM];
	real node2[DIM];
	real node3[DIM];
	real centroid[DIM];
};

//Bounding box which holds CFD cells and DEM particles
struct BdBox{
	int totalFaces; //number of contact faces
	int noOfParticles; //number of DEM particles
	int parts[NO_OF_PARTICLES_IN_BDCELL];
	struct wallFace face[NO_OF_FACES];
	real faceCentroid[NO_OF_FACES*DIM];
};


//Particle
struct demParticle{
	real dt; //particle time step
	real currentTime; //current time
	real displacement; //if displacement > rMax update neighbourlist
	real dia, inert, mass, nrmDisp;
	real *pos, *angVel, *vel, *hisDisp, *force, *momentum;
	int neigh[NBSIZE];
	int contList[CNTSIZE];
	int wallCntList[WALLCNTSIZE];
	int noOfCnt,noOfPartColl, noOfWallColl, noOfWallCnt;
	real maxPartCollE, maxWallCollE;
	int incontact;
	int noOfNeigh, cordNo;
	//Particle *cfdp;
	int prevCellIndex, preFluidCell;
	
	short int insertable;
	short int active; 
	real haPp;//Hamarker constant
	real *dragForce;
	real dragFX, dragFY, dragFZ;
	real eCharge;
};


void printCPUTime();
real solidFraction(int ip);

void writeLogNum(char *infile, char *line, real num);
void writeLog3Num(char *infile, char *line, real v1, real v2, real v3);
real readInputVelocity(char *infile);
void readGeom(char *infile, real *ductxmin, real *ductxmax, real *ductxedge1, real *ductxedge2, real *ductymin, 
            real *ductymax, real *ductzmin, real *ductzmax, real *ductzedge);

void readDomain(char *infile);
void test(Tracked_Particle *p, Thread *t);

int *allocateIntArray(int size);
real *allocateDoubleArray(int size);

struct BdBox *allocateBdBoxArray(int size);
struct demParticle *allocatePar(int np);
struct cfdCell *allocateCFDCell(int nc);

real partVol(int p);
void insertToBdBox(int p, int cI);
void addToBdBox();
void addFaceToBdBox();
void readInput(char *infile, int *np, real *dens, real *ymod, 
			real *pois, real *sfc, real *rec, real *dmpn, real *rf,
			real *cyldia, real *dt, int *nW, int *updateDPM, real *maxVel);
void findRec(FILE *inFile, char* strDest);
void diaInput(char *diaFile, struct demParticle *par, int *np);
void readWalls(char *infile, int *walls);

void writeTec();
void demInit();
void buildDEMCFDCellMap();
void copyDEMInfo();

void writeDump();
void demSave();
void  writeFluidVelocity();
void cordSave();
void cellSave();
void allocate();
void deallocate();

void updateForce(Tracked_Particle *p);
void partContactForce(int ip, int jp, real nrmDsp);
void findContactFromDPM(Particle *p);
void findContactFromMesh(Particle *p);
void findContactFromBoundary(Particle *p);
void calculateDragForce(Particle *p);
void boundaryContactForce(int pI, real *n1, real *n2, real *n3, real *uVec);
void ppVWForce(int ip, int jp, real vGap);
void pWallVWForce(int p, real vGap, real *uVec);
void ppElectForce(int ip, int jp, real gap, real *uVec);
void ppCapillaryForce(int ip, int jp, real gap);
void charge(int p, int jp);

void assignGravity();
real getOverlap(real *parPos, real dia, real *n1, real *n2, real *n3, real *uVec);

void neighbourContactForce(int pI);
void surfaceContactForce(int p, real nrmDisp, real *uVec);
void vecAdd(real *v1, real *v2, real *vec);
void crossProd(real *v1, real *v2, real *vec);
void vecSub(real *v1, real *v2, real *vec);
real relVel(int ip, int jp);
void unitVec(real *v, real *vec);
void sclMult(real scl, real *vec);
void sclVecMult(real scl, real *inVec, real *outVec);
real vecMag(real *vec);
void projVec(real *v1, real *n, real *vec, int type);
real dotProduct(real *v1, real *v2);
void getUnitVector(real *v1, real *v2, real *v3, real *uVec);

void initialize(real *sortedList, int *sortedParIndex, int *cellSE, int np,
    real *pos, real *parDia);

int insertable(int ip, int jp);
void addNeighbour(int  ip, int jp);
void updateNeighbourList(int p);
void deleteNeighbour(int ip, int jp);
void addContact(int ip, int cnt);
void deleteContact(int ip, int cnt);
void addWallContact(int ip, int cnt);
void deleteWallContact(int ip, int cnt);
int contactExist(int ip, int cnt);
void update(real *pX, int np);

real getCenterDist(int ip, int jp);
void setReduceUnits();
void updateBBFluidVel();
void insertToBdBox(int p, int cI);
void deleteParticle(int p, int cI);
void forceCalculation(Particle *p);
void updatePosition(Particle *p);

void run();

#endif 
