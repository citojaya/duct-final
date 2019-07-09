
#include "common.h"
#include "surf.h"
#include "mem.h"


// #define C2 100.0
// #define NUM_UDM 3

static int udm_offset = UDM_UNRESERVED;

/* Initialize dem information
Executes only once at the begining of FLUEMT*/
DEFINE_INIT(my_init, d)
{ 
  writeLogLine("logfile3.log", "INIT\n");
  initialized = 0;
  noOfCFDCells = 0;
  noOfWallFaces = 0;

  cell_t c;
  face_t f;
  real x[ND_ND];
  
  Thread *t;
  //Assign cell_id and store in UDMI
  thread_loop_c (t,d)
    begin_c_loop(c,t)
      C_UDMI(c,t,0) = noOfCFDCells; //initialize cfdcell solid volume
      C_CENTROID(x,c,t);
      noOfCFDCells++;

      int n;
      c_face_loop(c,t,n) 
      { 
        Thread *tf = C_FACE_THREAD(c,t,n); 
        if(THREAD_TYPE(tf) == THREAD_F_WALL) 
        { 
            f = C_FACE(c,t,n);
            F_UDMI(f,t,1) = noOfWallFaces;
            noOfWallFaces++;
            //F_CENTROID(x,f,tf);
        }
      }
    
    
    end_c_loop(c,t)

  cfdcell = allocateCFDCell(noOfCFDCells);
  
 for(int i=0; i<noOfCFDCells; i++){
        cfdcell[i].porosity = 0.99; //initial cfdcell porosity is set to 1
        cfdcell[i].solidVol = 0.0;
        cfdcell[i].dragFX = 0.0;
        cfdcell[i].dragFY = 0.0;
        cfdcell[i].dragFZ = 0.0;
        cfdcell[i].noOfParts = 0;
  }

  //initialize DEM
  demInit();
  writeLogNum("logfile4.log","NUM OF WALL FACES ",noOfWallFaces);
}

/*
Initialize DPM particle information
Executes only once on injected particles
*/
DEFINE_DPM_INJECTION_INIT(solid_paritcles, I)
{
  writeLogLine("logfile3.log", "INJECTION INITI\n");

  Injection *I2;
  Injection *Ilist = Get_dpm_injections();
  
  np = 0;
  loop(I2,Ilist)
  {
    Particle *p;
    loop(p,I2->p_init)
    {
      //Set particle mass according to DEM density input
      P_MASS(p) = (4.0/3.0)*PI*pow((0.5*P_DIAM(p)),3.0)*dens;
      np++;
    }
  }

  //Setup DEM scaling 
  setReduceUnits();

  //Insert particles to cell 
  loop(I2,Ilist)
  {
    Particle *p;
    loop(p,I2->p_init)
    {
      //Insert to cell
      demPart[p->part_id].pos[0] =  P_POS(p)[0]*lengthFactor;
      demPart[p->part_id].pos[1] =  P_POS(p)[1]*lengthFactor;
      demPart[p->part_id].pos[2] =  P_POS(p)[2]*lengthFactor;
      demPart[p->part_id].vel[0] = P_VEL(p)[0]*velocityFactor;
      demPart[p->part_id].vel[1] = P_VEL(p)[1]*velocityFactor;
      demPart[p->part_id].vel[2] = P_VEL(p)[2]*velocityFactor;
      demPart[p->part_id].dia = P_DIAM(p)*lengthFactor;
      demPart[p->part_id].mass = P_MASS(p)*massFactor;
      demPart[p->part_id].inert = 2.0*demPart[p->part_id].mass*pow(0.5*demPart[p->part_id].dia,2)/5.0; 
      demPart[p->part_id].vel[0] = 0.0;
      demPart[p->part_id].vel[1] = 0.0;
      demPart[p->part_id].vel[2] = 0.0;
    }
  }

  //insert to cells
  //Do not place this block in the above loop. Particle allocation to cells will not work. 
  //Therefore need a seperate loop. what a shame!
  loop(I2,Ilist)
  {
    Particle *p;
    loop(p,I2->p_init)
    {
      int iIndex = ceil((demPart[p->part_id].pos[0]-xmin)/domainDx); //ceil gives upper value
      int jIndex = ceil((demPart[p->part_id].pos[1]-ymin)/domainDy); //but we need lower value
      int kIndex = ceil((demPart[p->part_id].pos[2]-zmin)/domainDz); // therefore -1
      int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;  

      insertToBdBox(p->part_id, cellIndex);//add particle to the cell
    }
  }
 
  Thread *t;
  cell_t c;
  loop(I2,Ilist)
  {
    Particle *p;
    loop(p,I2->p_init)
    {
      t = P_CELL_THREAD(p);
      c = P_CELL(p);
      int cI = C_UDMI(c,t,0);
      //Update neighbour list
      updateNeighbourList(p->part_id);
    }
  }

  int maxBDSize = 0;
  for(int i=0; i<xDiv*yDiv*zDiv; i++){
    if(bdBox[i].noOfParticles > maxBDSize){
      maxBDSize = bdBox[i].noOfParticles;
    }
  }
  writeLogNum("logfile2.log","maxBDSize ",maxBDSize);

  int maxNBSize = 0;
  for(int i=0; i<np; i++){
    if(maxNBSize, demPart[i].noOfNeigh > maxNBSize){
      maxNBSize = demPart[i].noOfNeigh;
    }
  }

  //initialized = 1;
  writeLogNum("logfile2.log","maxNBSize ",maxNBSize);

  writeLogNum("logfile2.log","XMIN ",xmin/lengthFactor);
  writeLogNum("logfile2.log","YMIN ",ymin/lengthFactor);
  writeLogNum("logfile2.log","ZMIN ",zmin/lengthFactor);
  writeLogLine("logfile2.log","----------------\n");
  writeLogNum("logfile2.log","XMAX ",xmax/lengthFactor);
  writeLogNum("logfile2.log","YMAX ",ymax/lengthFactor);
  writeLogNum("logfile2.log","ZMAX ",zmax/lengthFactor);
  writeLogLine("logfile2.log","----------------\n");
  writeLogNum("logfile2.log","X LENGTH ",(xmax-xmin)/lengthFactor);
  writeLogNum("logfile2.log","Y LENGTH ",(ymax-ymin)/lengthFactor);
  writeLogNum("logfile2.log","Z LENGTH ",(zmax-zmin)/lengthFactor);
  writeLogLine("logfile2.log","----------------\n");
  writeLogNum("logfile2.log","DOMAIN DX ",domainDx/lengthFactor);
  writeLogNum("logfile2.log","DOMAIN DY ",domainDy/lengthFactor);
  writeLogNum("logfile2.log","DOMAIN DZ ",domainDz/lengthFactor);
  writeLogLine("logfile2.log","----------------\n");
  writeLogNum("logfile2.log","DOMAIN XDIV ",xDiv);
  writeLogNum("logfile2.log","DOMAIN YDIV ",yDiv);
  writeLogNum("logfile2.log","DOMAIN ZDIV ",zDiv);
  writeLogLine("logfile2.log","----------------\n");

  //fflush(stdout);
}

DEFINE_SOURCE(xmom_source,c,t,dS,eqn)
{
  real source = 0.;
  if(initialized)
  {
    int cellIndex = C_UDMI(c,t,0);
    source = cfdcell[cellIndex].dragFX/C_VOLUME(c,t);
  }
  dS[eqn] = 0.0;
  return source;
}

DEFINE_SOURCE(ymom_source,c,t,dS,eqn)
{
  real source = 0.;
  if(initialized)
  {
    int cellIndex = C_UDMI(c,t,0);
    source = cfdcell[cellIndex].dragFY/C_VOLUME(c,t);

  }
  dS[eqn] = 0.0;
  return source;
}

DEFINE_SOURCE(zmom_source,c,t,dS,eqn)
{
  real source = 0.;
  if(initialized)
  {
    int cellIndex = C_UDMI(c,t,0);
    source = cfdcell[cellIndex].dragFZ/C_VOLUME(c,t);
  }
  dS[eqn] = 0.0;
  return source;
}


DEFINE_PROFILE(inlet_x_velocity, thread, position)
{
  real x[ND_ND]; /* this will hold the position vector */
  real y, h;
  face_t f;

  real maxVel = readInputVelocity("infile");
  h = 0.01; /* inlet height in m */
  
  //writeLogNum("logfile3.log", "inlet velocity ",maxVel);
  begin_f_loop(f,thread)
  {
    F_CENTROID(x, f, thread);
    y = maxVel*4.*(0.25-pow((0.5-x[2]/h),2));
    F_PROFILE(f, thread, position) = y;
  }
  end_f_loop(f, thread)
}

DEFINE_PROFILE(bed_por,thread, np)
{
  cell_t c;
  real x[ND_ND];
  begin_c_loop(c,thread)
  {

    C_PROFILE(c,thread,np) = 1.0;//cfdcell[cI].porosity;
    if(initialized)
    {
      int cI = C_UDMI(c,thread,0);
      //cfdcell[cI].porosity = 1.0 - cfdcell[cI].solidVol/(C_VOLUME(c,thread));
      cfdcell[cI].porosity = fmax(0.2, 1.0 - cfdcell[cI].solidVol/(C_VOLUME(c,thread)));
      C_PROFILE(c,thread,np) = cfdcell[cI].porosity; 
    }
      
  }end_c_loop(c,thread)
}

DEFINE_PROFILE(viscous_resis,thread,i)
{
  cell_t c;
  real x[ND_ND];

    begin_c_loop(c,thread)
    {
      /***********************************************
        viscous_resistance = 150*(1-por)^2/(Dp^2*por^3)
      ************************************************/
      real instPor = 1.0;
      if(initialized)
      {
        int cI = C_UDMI(c,thread,0);
        instPor = cfdcell[cI].porosity;
        C_PROFILE(c,thread,i) = 0.01*(150.*pow((1-instPor),2)/(pow(largestParDia,2)*pow(instPor,3)));
      }
      
    }end_c_loop(c,thread)
  
}

DEFINE_PROFILE(inertial_resis,thread,i)
{
  cell_t c;
  real x[ND_ND];
    begin_c_loop(c,thread)
    {
      /***********************************************
        inertial_resistance = 3.5*(1-por)/(Dp*por^3)
      ************************************************/
      real instPor = 1.0;
      if(initialized)
      {
        int cI = C_UDMI(c,thread,0);
        instPor = cfdcell[cI].porosity;
        C_PROFILE(c,thread,i) = 3.5*(1-instPor)/(largestParDia*pow(instPor,3));
     }
    }end_c_loop(c,thread)
}

/*
When particles are in contact with surface this macro is executed by FLUENT
Loop through all particles and update contact surface normal vector and contact surface node positions
Only applicable when contact surface is triangular
*/
DEFINE_DPM_BC(bc_abort, p, t, f, f_normal, dim)
{
  long int ip = p->part_id;

  //return PATH_END;
  return PATH_ABORT;
  //return PATH_ACTIVE;
}


DEFINE_EXECUTE_AT_END(execute_at_end)
{
  cfdcycles++;
  initialized = 1;
  updateSourceTerm = 1;

  Injection *I;
  Injection *Ilist = Get_dpm_injections();

  Domain *d = Get_Domain(1);
  Thread *t;
  cell_t c;
  thread_loop_c (t,d)
  begin_c_loop(c,t)
    int cI = C_UDMI(c,t,0); //initialize cfdcell solid volume
    cfdcell[cI].solidVol = 0.0;
    cfdcell[cI].dragFX = 0.0;
    cfdcell[cI].dragFY = 0.0;
    cfdcell[cI].dragFZ = 0.0;
    cfdcell[cI].noOfParts = 0;
  end_c_loop(c,t)

  for(int i=0; i<iter; i++)
  {
    demTime += timeStep;
    loop(I,Ilist)
    {
      Particle *p;
      loop(p,I->p)
      { 
        if(demPart[p->part_id].active == 1){
          forceCalculation(p);
          maxCnt = fmax(maxCnt,demPart[p->part_id].noOfCnt);
        }

      }
    }

    loop(I,Ilist)
    {
      Particle *p;
      loop(p,I->p)
      { 
        if(demPart[p->part_id].active == 1){
          updatePosition(p);
        }
      }
    }

    Thread *tc;
    cell_t pc;
    
    loop(I,Ilist)
    {
      Particle *p;
      loop(p,I->p)
      {
          tc = P_CELL_THREAD(p);
          pc = P_CELL(p);
          int cI = C_UDMI(pc,tc,0);
          updateNeighbourList(p->part_id);
      }
    }
    updateSourceTerm = 0;
  }//end if iteration


  if(cfdcycles%10 == 0)
  {
    demSave();
    printCPUTime();
  }
  if(cfdcycles%100 == 0){
    writeInjectionFile("initial.inj");
    writeFluidVelocity();
  }
}

/*Execute at the begining of every fluid iteration*/
// DEFINE_ADJUST(define_adjust, d)
// {
//     if(updateDPM == 0){
//       Injection *I;
//       Injection *Ilist = Get_dpm_injections();

//       loop(I,Ilist)
//       {
//         Particle *p;
//         loop(p,I->p)
//         { 
//           P_POS(p)[0] = demPart[p->part_id].pos[0]/lengthFactor;
//           P_POS(p)[1] = demPart[p->part_id].pos[1]/lengthFactor;
//           P_POS(p)[2] = demPart[p->part_id].pos[2]/lengthFactor;
//        }
//       }
//     }
// }

DEFINE_ON_DEMAND(demand)
{
  writeInjectionFile("initial.inj");
}


DEFINE_EXECUTE_AT_EXIT(execute_at_exit)
{
  deallocate();
  //fflush(stdout);
}

/* Rotate the capsule for a given speed*/
// DEFINE_CG_MOTION(piston,dt,vel,omega,time,dtime)
//  {
//     omega[2]=50;
//  } 
