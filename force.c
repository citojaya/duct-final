#include "common.h"

/* Calculate particle-particle, particle-wall contact forces and drag force due to fluid*/
void forceCalculation(Particle *p)
{
    real dt = timeStep;
    demPart[p->part_id].force[0] = 0.0; //gravitational force
    demPart[p->part_id].force[1] = 0.0;//gravitational force
    demPart[p->part_id].force[2] = -demPart[p->part_id].mass; 
    demPart[p->part_id].momentum[0] = 0.0;
    demPart[p->part_id].momentum[1] = 0.0;
    demPart[p->part_id].momentum[2] = 0.0;

    findContactFromDPM(p);
    //findContactFromMesh(p);
    //findContactFromBoundary(p);
 
    //Find particle-particle contact force
    neighbourContactForce(p->part_id);

    //Find drag force on particles from fluid 
    calculateDragForce(p);
}

/*Particle-particle vanderwal force*/
void ppVWForce(int ip, int jp, real vGap){
    
    real ipDia = demPart[ip].dia;
    real jpDia = demPart[jp].dia;
    real vGapMn = 1.0e-9*lengthFactor;
   
    real ijHa = sqrt(demPart[ip].haPp*demPart[jp].haPp);
    vGap = fmax(vGap, vGapMn);
    
    real rStar = ipDia*jpDia/(ipDia+jpDia);
    real fv = -ijHa*rStar/(6.*vGap*vGap);
    maxVF = fmax(maxVF,-fv);

    real *uVec = allocateDoubleArray(DIM);
    vecSub(demPart[ip].pos, demPart[jp].pos, uVec);
    unitVec(uVec, uVec);

    demPart[ip].force[0] += uVec[0]*fv;
    demPart[ip].force[1] += uVec[1]*fv;
    demPart[ip].force[2] += uVec[2]*fv;
    free(uVec);
}

/*Particle-wall vanderwal force*/
void pWallVWForce(int p, real vGap, real *uVec){
    real vGapMn = 1.0e-9*lengthFactor;
    vGap = fmax(vGap, vGapMn);
    real fv = -haPW*demPart[p].dia*0.5/(6.*vGap*vGap);
    demPart[p].force[0] += uVec[0]*fv;
    demPart[p].force[1] += uVec[1]*fv;
    demPart[p].force[2] += uVec[2]*fv;   
    // if(cfdcycles%10 == 0){
    //     writeLogNum("logfile5.log","P-WALL FORCE ",fv*1e12);
    // }  
}

/*Partcile-particle capillary force*/
void ppCapillaryForce(int ip, int jp, real gap){
    
    real s_rupture = (1+0.5*cont_ang)*pow(liq_vol,1/3);
    if(gap/lengthFactor < s_rupture){
        real *uVec = allocateDoubleArray(DIM);
        vecSub(demPart[ip].pos, demPart[jp].pos, uVec);
        unitVec(uVec, uVec);
        real sepMin = 5.0e-6; //((Hornbaker, Albert et al. 1997, Nase, Vargas et al. 2001)
        real separation = fmax(sepMin, gap/lengthFactor);
        real rStar = (demPart[ip].dia*demPart[jp].dia/(demPart[ip].dia+demPart[jp].dia))/lengthFactor;
        real capf = -2.*PI*rStar*surf_tens*cos(cont_ang)/(1+1/sqrt(1+2.*liq_vol/(PI*rStar*pow(separation,2))-separation));
        // if(ip == 1)
        //    writeLogNum("logfile4.log","CAP-FORCE ",1e6*fe/forceFactor);
        // demPart[ip].force[0] += uVec[0]*capf*forceFactor;
        // demPart[ip].force[1] += uVec[1]*capf*forceFactor;
        // demPart[ip].force[2] += uVec[2]*capf*forceFactor;
        free(uVec);  
    }         
}

/*
Calculate charge on particle
The charge model is based on ref (Matsusaka et al., 2000 and Pei et al., 2016)
@param
p - int particle number
*/
void charge(int p, int jp){
    real velMag = vecMag(demPart[p].vel);
    real parR = demPart[p].dia*0.5;
    //k0 = z0/(4.*math.pi*var.eps*pow(parR/var.lengthFactor,2)) //constant
    real k0 = imageConst*Zs/permitivity/(4.0*PI*pow(parR/lengthFactor,2)); //(Pei et al., 2016)
    real vDash = vDash = k0*(demPart[p].eCharge);
   
    real emod = alpha*(1.0-pois*pois)/ymod;
    real S = 1.36*pow(emod, 2.0/5.0)*pow(largestParDensity*densityFactor,2.0/5.0)*pow(parR*2.0,2)*pow(velMag,4.0/5.0);
    S = S/(lengthFactor*lengthFactor); //convert to SI units
    
    if(p < 2000){
        V1 = V1*0.3;
    }
    else if(p >2000 && p<4000){
        V1 = V1*0.6;
    }
    real deltaV = V1 - V2 - vDash;
    real deltaQ = ks*S*deltaV; //
    demPart[p].eCharge += deltaQ;
    if(demPart[p].eCharge > maxCharge){
        maxChargePart = p;
    }
    maxCharge = fmax(fabs(demPart[p].eCharge), maxCharge);
   
}

/*Particle-particle electrostatic force*/
void ppElectForce(int ip, int jp, real gap, real *uVec){
    maxCharge = fmax(maxCharge, demPart[ip].eCharge);   
    real fe = 0.0;
    real r1r2 = 0.0;
    //particle-particle
    if(jp > 0){
        r1r2 = getCenterDist(ip,jp)/lengthFactor;
        fe = demPart[ip].eCharge*demPart[jp].eCharge*forceFactor/(4.*PI*permitivity*r1r2*r1r2);
    }
    else{//particle-wall
        r1r2 = 0.5*demPart[ip].dia + gap;
        //Since uVec is always towards center fe should be minus in order to get particle-wall attraction
        fe = -demPart[ip].eCharge*demPart[ip].eCharge*forceFactor/(4.*PI*permitivity*r1r2*r1r2);
    }
    //if(ip == 21027){
    maxES = fmax(maxES, fabs(fe));
    //}
    // if(ip == 1)
    //     writeLogNum("logfile4.log","ES-FORCE ",1e6*maxES/forceFactor);  
 
    // if(cfdcycles%10 == 0){
    //     writeLog3Num("logfile5.log","ES CHARGE ",demPart[ip].eCharge*1e20, maxES*1e12, 0.0);
    // }  

    demPart[ip].force[0] += esfTrue*uVec[0]*fe;
    demPart[ip].force[1] += esfTrue*uVec[1]*fe;
    demPart[ip].force[2] += esfTrue*uVec[2]*fe;      
}

/*
Calculate particle-wall contact forces
param:
ip - ith particle
nrmDisp - overlap
uVec - unit vector normal to contact surface
*/
void surfaceContactForce(int p, real nrmDisp, real *uVec){
    //real rStar = 0.5*P_DIAM(p)*lengthFactor;
    real rStar = 0.5*demPart[p].dia;
    sclVecMult(-0.5*demPart[p].dia,uVec,ipRVec);
    
    crossProd(demPart[p].angVel,ipRVec,rotVel);
    vecAdd(demPart[p].vel,rotVel,ipCntPntVel);

    real *relVel = allocateDoubleArray(DIM);
 
    sclVecMult(1.0,demPart[p].vel,relVel);
    
    real nrmVel = dotProduct(relVel,uVec);
    free(relVel);
    
    sclVecMult(1.0,ipCntPntVel,cntPntVel);

    real *totalForce = allocateDoubleArray(DIM);
    real *momentum = allocateDoubleArray(DIM);

    real nrmCntForce = elasticMod*sqrt(rStar*nrmDisp)*nrmDisp;
    real nrmDampForce = -dmpn*elasticMod*sqrt(rStar*nrmDisp)*nrmVel;

    real *cntPntDisp = allocateDoubleArray(DIM);
    sclVecMult(timeStep,cntPntVel,cntPntDisp);
    real *tipCntPntDisp = allocateDoubleArray(DIM);
    projVec(cntPntDisp, uVec, tipCntPntDisp, 0);
    real slidingDisp = vecMag(tipCntPntDisp);

    real *disp = allocateDoubleArray(DIM);

    //projVec(demPart[p].hisDisp, uVec, disp, 1);
    //vecAdd(disp, tipCntPntDisp, disp);
    sclVecMult(1.0,tipCntPntDisp,disp);
    dsmax = dsmaxCff*nrmDisp;
    dd = vecMag(disp);
    dti = vecMag(tipCntPntDisp);
    fdt = sfc*nrmCntForce;
    real *fdtVec = allocateDoubleArray(DIM);

    if(dd < 1e-6){
        dd = 0.0;
    }
    if(dti < 1e-6){
        dti = 0.0;
    }
    if(dd >= dsmax){
        sclMult(dsmax/dd,disp);
        if(dti != 0){
            sclVecMult(-fdt/dti,tipCntPntDisp,fdtVec);
        }
        else{
            sclMult(0.0,fdtVec);
        }
    }
    else{
        if(dd != 0.0){
            fdt = fdt*(1.0 - pow((1.0 - dd/dsmax),1.5));
            sclVecMult(-fdt/dd,disp,fdtVec);
        }
        else{
            sclMult(0.0,fdtVec);
        }
    }

    //sclVecMult(1.0,disp, demPart[p].hisDisp);

    //sum of forces
    real nrmForce = (nrmCntForce + nrmDampForce);
    sclVecMult(nrmForce, uVec, totalForce);

    vecAdd(totalForce, fdtVec, totalForce);
    crossProd(ipRVec, totalForce, momentum);
    real *rotMom = allocateDoubleArray(DIM);
    sclVecMult(0.5*rf*demPart[p].dia*nrmCntForce, demPart[p].angVel, rotMom);
    vecAdd(momentum, rotMom, momentum);

    //Update force and momentum on particle
    vecAdd(demPart[p].force, totalForce, demPart[p].force);
    vecAdd(demPart[p].momentum, momentum, demPart[p].momentum);

    free(rotMom);   
    free(fdtVec);
    free(tipCntPntDisp);
    free(cntPntDisp);
    free(totalForce);
    free(momentum);
    free(disp);
}

/*
Calculate interparticle forces
param:
ip - ith particle
jp - neighbour particle
nrmDisp - overlap
*/
void partContactForce(int ip, int jp, real nrmDisp){
    real rStar = 0.5*demPart[ip].dia*demPart[jp].dia/(demPart[ip].dia+demPart[jp].dia);
    vecSub(demPart[ip].pos, demPart[jp].pos, ijVec);
    unitVec(ijVec,uVec);

    sclVecMult(-0.5*demPart[ip].dia,uVec,ipRVec);
    sclVecMult(0.5*demPart[jp].dia,uVec,jpRVec);

    crossProd(demPart[ip].angVel,ipRVec,rotVel);
    vecAdd(demPart[ip].vel,rotVel,ipCntPntVel);

    crossProd(demPart[jp].angVel,jpRVec,rotVel);
    vecAdd(demPart[jp].vel,rotVel,jpCntPntVel);

    real *relVel = allocateDoubleArray(DIM);
    vecSub(demPart[ip].vel,demPart[jp].vel,relVel);

    real nrmVel = dotProduct(relVel,uVec);
    vecSub(ipCntPntVel,jpCntPntVel,cntPntVel);
    free(relVel);
 
    real *totalForce = allocateDoubleArray(DIM);
    real *momentum = allocateDoubleArray(DIM);

    real nrmCntForce = elasticMod*sqrt(rStar*nrmDisp)*nrmDisp;
    real nrmDampForce = -dmpn*elasticMod*sqrt(rStar*nrmDisp)*nrmVel;

    real *cntPntDisp = allocateDoubleArray(DIM);
    sclVecMult(timeStep,cntPntVel,cntPntDisp);
    real *tipCntPntDisp = allocateDoubleArray(DIM);
    projVec(cntPntDisp, uVec, tipCntPntDisp, 0);
    real slidingDisp = vecMag(tipCntPntDisp);

    real *disp = allocateDoubleArray(DIM);

    sclVecMult(1.0,tipCntPntDisp,disp);

    //projVec(demPart[ip].hisDisp, uVec, disp, 1);
    //vecAdd(disp, tipCntPntDisp, disp);
    dsmax = dsmaxCff*nrmDisp;
    dd = vecMag(disp);
    dti = vecMag(tipCntPntDisp);
    fdt = sfc*nrmCntForce;

    real *fdtVec = allocateDoubleArray(DIM);
 
    if(dd < 1e-6){
        dd = 0.0;
    }
    if(dti < 1e-6){
        dti = 0.0;
    }
    if(dd >= dsmax){
        sclMult(dsmax/dd,disp);
        if(dti != 0){
            sclVecMult(-fdt/dti,tipCntPntDisp,fdtVec);
        }
        else{
            sclMult(0.0,fdtVec);
        }
    }
    else{
        //sclVecMult(1.0,disp, demPart[ip].hisDisp);
        if(dd != 0.0){
            fdt = fdt*(1.0 - pow((1.0 - dd/dsmax),1.5));
            sclVecMult(-fdt/dd,disp,fdtVec);
        }
        else{
            sclMult(0.0,fdtVec);
        }
    }

    //sum of forces
    real nrmForce = (nrmCntForce + nrmDampForce);
    //writeLog("logfile2.log","nrmForce ",nrmForce);
    sclVecMult(nrmForce, uVec, totalForce);
    vecAdd(totalForce, fdtVec, totalForce);
    crossProd(ipRVec, totalForce, momentum);
    real *rotMom = allocateDoubleArray(DIM);
    sclVecMult(0.5*rf*demPart[ip].dia*nrmCntForce, demPart[ip].angVel, rotMom);
    vecAdd(momentum, rotMom, momentum);

    //Update force and momentum on particle
    vecAdd(demPart[ip].force, totalForce, demPart[ip].force);
    vecAdd(demPart[ip].momentum, momentum, demPart[ip].momentum);

    free(rotMom);
    free(fdtVec);
    free(tipCntPntDisp);
    free(cntPntDisp);
    free(totalForce);
    free(momentum);
    free(disp);
}

void calculateDragForce(Particle *p){
    real x[ND_ND];
    Thread *tc = P_CELL_THREAD(p);
    cell_t c = P_CELL(p);
    //Fluid drag force calculation
    real velFX = C_U(c,tc)*velocityFactor;  
    real velFY = C_V(c,tc)*velocityFactor;  
    real velFZ = C_W(c,tc)*velocityFactor;    
    real pGX = C_P_G(c,tc)[0]*pressureFactor/lengthFactor;
    real pGY = C_P_G(c,tc)[1]*pressureFactor/lengthFactor;
    real pGZ = C_P_G(c,tc)[2]*pressureFactor/lengthFactor;
    real velPX = demPart[p->part_id].vel[0];
    real velPY = demPart[p->part_id].vel[1];
    real velPZ = demPart[p->part_id].vel[2];
    real density = C_R(c,tc)*densityFactor;
    real visc = C_MU_L(c,tc)*massFactor/(lengthFactor*timeFactor);
    
    C_CENTROID(x,c,tc);
    int cellIndex = C_UDMI(c,tc,0);
   
    real instPor = cfdcell[cellIndex].porosity;
    instPor = fmin(instPor,0.99);
    
    real relVelMag = sqrt((velFX-velPX)*(velFX-velPX)+(velFY-velPY)*(velFY-velPY)+(velFZ-velPZ)*(velFZ-velPZ));

    real Re = instPor*demPart[p->part_id].dia*relVelMag*density/visc;
    // writeLog3Num("logfile6.log"," ",instPor,relVelMag/velocityFactor,density/densityFactor);
    // writeLogNum("logfile6.log","VISC ",visc*(lengthFactor*timeFactor/massFactor));
    // writeLogNum("logfile6.log","RE ",Re);
  
    real beeta, dCoeff;

	if (Re < 1000){
		dCoeff = 24.*(1.+0.15*pow(Re, 0.687))/Re;
	}
	else {
		dCoeff = 0.44;
	}
/*
    // -  Second method  Yi He paper-//

    if(instPor <= 0.8){
        beeta = 150.*(1.-instPor)*(1.-instPor)*visc/(instPor*pow(demPart[p->part_id].dia,2)) 
                + 1.75*(1.-instPor)*density*relVelMag/demPart[p->part_id].dia;
    }
    else if(instPor > 0.8){
        beeta = (3.0/4.0)*dCoeff*relVelMag*density*(1.0-instPor)*pow(instPor,-2.7)/demPart[p->part_id].dia;
    }
    real pfFX = (velFX-velPX)*partVol(p->part_id)*beeta/(1.0-instPor); //always instPor != 0
    real pfFY = (velFY-velPY)*partVol(p->part_id)*beeta/(1.0-instPor);
    real pfFZ = (velFZ-velPZ)*partVol(p->part_id)*beeta/(1.0-instPor);
*/

    // if(cfdcycles%2 == 0){
    //     writeLogNum("logfile5.log","POR ",instPor);
    // }

    // (Zhu et al., 2007)
    // real alpha = 3.7-0.65*exp(-0.5*pow((1.5-log10(Re)),2));
    // real pfFX = pow(instPor,-(alpha+1))*dCoeff*PI*density*(velFX-velPX)*fabs(velFX-velPX)*pow(demPart[p->part_id].dia,2)/8.0;
    // real pfFY = pow(instPor,-(alpha+1))*dCoeff*PI*density*(velFY-velPY)*fabs(velFY-velPY)*pow(demPart[p->part_id].dia,2)/8.0;
    // real pfFZ = pow(instPor,-(alpha+1))*dCoeff*PI*density*(velFZ-velPZ)*fabs(velFZ-velPZ)*pow(demPart[p->part_id].dia,2)/8.0;

 
    
    //real pfFX = pow((0.63+4.8/pow()),2)*density*
    // (Chu et al., )

    beeta = 3.7-0.65*exp(-0.5*pow((1.5-log10f(Re)),2));
    real A = pow(0.63+4.8/pow(Re,0.5),2); 
    real B = pow(instPor,(2.0-beeta));
    
    real pfFX = A*B*PI*density*(velFX-velPX)*fabs(velFX-velPX)*pow(demPart[p->part_id].dia,2)/8.0;
    real pfFY = A*B*PI*density*(velFY-velPY)*fabs(velFY-velPY)*pow(demPart[p->part_id].dia,2)/8.0;
    real pfFZ = A*B*PI*density*(velFZ-velPZ)*fabs(velFZ-velPZ)*pow(demPart[p->part_id].dia,2)/8.0;

    //real pfFX = A*B*PI*density*(velFX-velPX)*relVelMag*pow(demPart[p->part_id].dia,2)/8.0;
    //real pfFY = A*B*PI*density*(velFY-velPY)*relVelMag*pow(demPart[p->part_id].dia,2)/8.0;
    //real pfFZ = A*B*PI*density*(velFZ-velPZ)*relVelMag*pow(demPart[p->part_id].dia,2)/8.0;


    real pGFX = -pGX*PI*pow(demPart[p->part_id].dia,3)/6.0;
    real pGFY = -pGY*PI*pow(demPart[p->part_id].dia,3)/6.0;
    real pGFZ = -pGZ*PI*pow(demPart[p->part_id].dia,3)/6.0;

    //Update force on particles
    demPart[p->part_id].force[0]  += (pfFX + pGFX);
    demPart[p->part_id].force[1]  += (pfFY + pGFY); 
    demPart[p->part_id].force[2]  += (pfFZ + pGFZ);

    //Store drag force on particle for later use in calculating source term for fluid
    demPart[p->part_id].dragFX = pfFX + pGFX;
    demPart[p->part_id].dragFY = pfFY + pGFY;
    demPart[p->part_id].dragFZ = pfFZ + pGFZ;

    if(updateSourceTerm == 1){
        cfdcell[cellIndex].solidVol += partVol(p->part_id)/volumeFactor;
        cfdcell[cellIndex].noOfParts += 1;
        cfdcell[cellIndex].dragFX +=  -demPart[p->part_id].dragFX/forceFactor; 
        cfdcell[cellIndex].dragFY +=  -demPart[p->part_id].dragFY/forceFactor; 
        cfdcell[cellIndex].dragFZ +=  -demPart[p->part_id].dragFZ/forceFactor;
    }

}


