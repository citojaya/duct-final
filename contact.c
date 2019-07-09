#include "common.h"


/*Calculate overlap with contact surface
param:
parPos - particle center
dia - particle diameter
node1, node2, node3 - contact surface nodes
uVec - unit vector normal to surface
*/
real getOverlap(real *parPos, real dia, real *n1, real *n2, real *n3, real *uVec){ 
    real *v1 = allocateDoubleArray(DIM); //vector running from node1 to particles center
    real *v2 = allocateDoubleArray(DIM); //vector running from node1 to particles center
    real *ppDash = allocateDoubleArray(DIM); //vector running from particle center to projection
    vecSub(parPos,n1,v1);  
    projVec(v1, uVec, v2, 0);
    vecSub(v1,v2,ppDash);
    uVec[0] = ppDash[0];
    uVec[1] = ppDash[1];
    uVec[2] = ppDash[2];

    unitVec(uVec,uVec);
    real overlap = vecMag(ppDash)-0.5*dia;

    free(v1);
    free(v2);
    free(ppDash);
    return overlap;
}


/*Find contact forces with all neighbour particles*/
void neighbourContactForce(int pI){
    //Loop through all neighbour particles
    int count = 0;

    for(int i=0; i<demPart[pI].noOfCnt; i++){
        int cntPart = demPart[pI].contList[i];
        demPart[cntPart].incontact = 1; //cntPart is incontact with pI
    }

    for(int i=0; i<demPart[pI].noOfNeigh; i++){
        int jp = demPart[pI].neigh[i];
        real gap = getCenterDist(pI,jp)-(demPart[pI].dia+demPart[jp].dia)*0.5;
        //if(pI < jp){
            if(gap < 0.0){
                count += 1;
                partContactForce(pI,jp, -gap);
                deleteContact(pI,jp); //delete existing contact before adding
                addContact(pI,jp);
                //For each new contact calculate electric charge 
                if(demPart[jp].incontact == 0){
                    charge(pI, jp);
                    demPart[pI].noOfPartColl += 1;
                    real *uVec = allocateDoubleArray(DIM);
                    vecSub(demPart[pI].pos, demPart[jp].pos, ijVec);
                    unitVec(ijVec,uVec);
                    real *relVel = allocateDoubleArray(DIM);
                    vecSub(demPart[pI].vel,demPart[jp].vel,relVel);
                    real nrmVel = dotProduct(relVel,uVec);
                    free(uVec);
                    free(relVel);
                    demPart[pI].maxPartCollE = fmax(demPart[pI].maxPartCollE, 0.5*demPart[pI].mass*nrmVel*nrmVel);
                }
                demPart[jp].incontact = 0; //reset
            }
            else{
                demPart[jp].incontact = 0; //reset
                deleteContact(pI,jp);
            }
            if(gap < 100.e-9*lengthFactor)//activate vanderwaal force when gap<100nm
            {
                ppVWForce(pI, jp, gap);
                //ppVWForce(jp, pI, gap);
            }
            if(gap < demPart[pI].dia){
                 ppElectForce(pI, jp, gap, uVec);
            }
        //}
    }
    demPart[pI].cordNo = count;
}

/*Find particle-wall contact and calculate force using DPM cell location*/
void findContactFromDPM(Particle *p){
    int n;
    face_t f;
    Thread *tc = P_CELL_THREAD(p);
    cell_t c = P_CELL(p);
    Thread *tf;
    real x[ND_ND];
    real cell[ND_ND];

    if(demPart[p->part_id].preFluidCell != C_UDMI(c,tc,0)){
        demPart[p->part_id].noOfWallCnt = 0;
    }
    demPart[p->part_id].preFluidCell = C_UDMI(c,tc,0);

    c_face_loop(c,tc,n) 
    { 
        tf = C_FACE_THREAD(c,tc,n);
        if(THREAD_TYPE(tf) == THREAD_F_WALL) 
        { 
            f = C_FACE(c,tc,n);
            //writeLogNum("logfile4.log","FACE ID ",F_UDMI(f,tc,1));
            //F_CENTROID(x,f,tf);

            int m;
            int j=0;
            real *node1 = allocateDoubleArray(DIM);
            real *node2 = allocateDoubleArray(DIM);
            real *node3 = allocateDoubleArray(DIM);
            f_node_loop(f,tf,m)//Loop over face nodes
            {
                Node *node = F_NODE(f,tf,m);
                switch(j){
                    case 0:
                        node1[0] = NODE_X(node)*lengthFactor;
                        node1[1] = NODE_Y(node)*lengthFactor;
                        node1[2] = NODE_Z(node)*lengthFactor;
                    case 1:
                        node2[0] = NODE_X(node)*lengthFactor;
                        node2[1] = NODE_Y(node)*lengthFactor;
                        node2[2] = NODE_Z(node)*lengthFactor;
                    case 2:
                        node3[0] = NODE_X(node)*lengthFactor; 
                        node3[1] = NODE_Y(node)*lengthFactor;
                        node3[2] = NODE_Z(node)*lengthFactor;
                }
                j++; 
            }//end of loop over face nodes
             
            real *uVec = allocateDoubleArray(DIM);
            getUnitVector(node1,node2,node3,uVec);
            real gap = getOverlap(demPart[p->part_id].pos,demPart[p->part_id].dia, node1,node2, node3, uVec);
            
            if(gap < demPart[p->part_id].dia){
                 ppElectForce(p->part_id, 0, gap, uVec);
            
                //Calculate Vanderwaal force
                if(gap < 100.e-9*lengthFactor){pWallVWForce(p->part_id, gap, uVec);}   
                
                if(gap < 0) //If contact exists calculate contact force
                {
                    //writeLog3Num("logfile5.log","Face node ", p->part_id,0.0,1e6*gap/lengthFactor);
                    if(contactExist(p->part_id, F_UDMI(f,tc,1)) == 0){
                        //deleteWallContact(p->part_id, F_UDMI(f,tc,1));
                        //if(p->part_id == 9642){
                        charge(p->part_id, -1);
                        addWallContact(p->part_id, F_UDMI(f,tc,1));
                        demPart[p->part_id].noOfWallColl += 1;
                        //}
                        //real velMag = vecMag(demPart[p->part_id].vel);
                        real *relVel = allocateDoubleArray(DIM);
                        sclVecMult(1.0,demPart[p->part_id].vel,relVel);
                        real nrmVel = dotProduct(relVel,uVec);
                        free(relVel);
                        demPart[p->part_id].maxWallCollE = fmax(demPart[p->part_id].maxWallCollE,0.5*demPart[p->part_id].mass*nrmVel*nrmVel);

                    }
                    surfaceContactForce(p->part_id, -gap, uVec);
                }
            }
            if(gap > 0){
                //if(p->part_id == 9642){
                    
                deleteWallContact(p->part_id,F_UDMI(f,tc,1));
                //}
            }
            free(uVec);
            free(node1);
            free(node2);
            free(node3);
        }
    }
}

/*Find particle-wall contact from solid boundary*/
void findContactFromBoundary(Particle *p){
    
//working code
    checkYContact(p, ductymin*lengthFactor, ductymax*lengthFactor);
    //check top bound contact
    if(demPart[p->part_id].pos[2] > (ductzmax*lengthFactor-0.5*demPart[p->part_id].dia)){
        checkZContact(p, 0, ductzmax*lengthFactor);
    }
    else if(demPart[p->part_id].pos[0] <= ductxedge1*lengthFactor || 
        demPart[p->part_id].pos[0] >= ductxedge2*lengthFactor){
        checkZContact(p, 0,ductzmax*lengthFactor);
    }
    else if(demPart[p->part_id].pos[2] <= 0){
         checkXContact(p, ductxedge1*lengthFactor, ductxedge2*lengthFactor);
         checkZContact(p, ductzmin*lengthFactor, ductzmax*lengthFactor);
    }
 
    else{    
        //check for x=100 edge
        if(demPart[p->part_id].pos[0] > ductxedge1*lengthFactor && 
            demPart[p->part_id].pos[0] < ductxedge1*lengthFactor+0.5*demPart[p->part_id].dia){
            if(demPart[p->part_id].pos[2] > 0 && demPart[p->part_id].pos[2] < 0.5*demPart[p->part_id].dia){
                real *iVec = allocateDoubleArray(DIM);
                real *jVec = allocateDoubleArray(DIM);
                
                iVec[0] = demPart[p->part_id].pos[0];
                iVec[1] = 0.0;
                iVec[2] = demPart[p->part_id].pos[2];
                jVec[0] = ductxedge1*lengthFactor;
                jVec[1] = 0.0;//demPart[p->part_id].pos[1];
                jVec[2] = 0.0;
            
                real *ijVec = allocateDoubleArray(DIM);
                vecSub(iVec, jVec, ijVec);
                real gap = vecMag(ijVec);
                unitVec(ijVec, uVec);

                if(gap < 0.5*demPart[p->part_id].dia){
                    surfaceContactForce(p->part_id, (0.5*demPart[p->part_id].dia-gap), uVec);
                }

                free(ijVec);
                free(iVec);
                free(jVec); 
                
            }
        }
     
        //check for x=150 edge
        else if(demPart[p->part_id].pos[0] < ductxedge2*lengthFactor && 
            demPart[p->part_id].pos[0] > ductxedge2*lengthFactor-0.5*demPart[p->part_id].dia){
            
            if(demPart[p->part_id].pos[2] > 0 && demPart[p->part_id].pos[2] < 0.5*demPart[p->part_id].dia){
                real *iVec = allocateDoubleArray(DIM);
                real *jVec = allocateDoubleArray(DIM);
                
                iVec[0] = demPart[p->part_id].pos[0];
                iVec[1] = 0.0;
                iVec[2] = demPart[p->part_id].pos[2];
                jVec[0] = ductxedge2*lengthFactor;
                jVec[1] = 0.0;//demPart[p->part_id].pos[1];
                jVec[2] = 0.0;
       
                real *ijVec = allocateDoubleArray(DIM);
                vecSub(iVec, jVec, ijVec);
                real gap = vecMag(ijVec);
                unitVec(ijVec, uVec);

                if(gap < 0.5*demPart[p->part_id].dia){
                    surfaceContactForce(p->part_id, (0.5*demPart[p->part_id].dia-gap), uVec);
                 }
                free(ijVec);
                free(iVec);
                free(jVec); 
            }
        }
    }
}

/*Find particle-wall contact and calculate force using trianglusr boundary Mesh*/
void findContactFromMesh(Particle *p){
    int cI = demPart[p->part_id].prevCellIndex;
    int nearestFaceIndex = -1;
    real minDist = 100000;
    real ipX = demPart[p->part_id].pos[0];
    real ipY = demPart[p->part_id].pos[1];
    real ipZ = demPart[p->part_id].pos[2];
    writeLogNum("logfile2.log", "particle cell index ", cI);
    for(int i=0; i<bdBox[cI].totalFaces; i++){
        
        real jpX = bdBox[cI].face[i].centroid[0]*lengthFactor;
        real jpY = bdBox[cI].face[i].centroid[1]*lengthFactor;
        real jpZ = bdBox[cI].face[i].centroid[2]*lengthFactor;

        real dist = sqrt((ipX-jpX)*(ipX-jpX) + (ipY-jpY)*(ipY-jpY) + (ipZ-jpZ)*(ipZ-jpZ));
        if(dist < minDist){
            nearestFaceIndex = i;
            minDist = dist;
        }
    }
    if(nearestFaceIndex != -1){
        bdBox[cI].face[nearestFaceIndex].node1[0] = bdBox[cI].face[nearestFaceIndex].node1[0]*lengthFactor;
        bdBox[cI].face[nearestFaceIndex].node1[1] = bdBox[cI].face[nearestFaceIndex].node1[1]*lengthFactor;
        bdBox[cI].face[nearestFaceIndex].node1[2] = bdBox[cI].face[nearestFaceIndex].node1[2]*lengthFactor;

        //writeLogNum("logfile2.log", "nearestFaceIndex ", cI);
            //real gap = demPart[p->part_id].pos[1] - demPart[p->part_id].dia*0.5;
        real *uVec = allocateDoubleArray(DIM);
        getUnitVector(bdBox[cI].face[nearestFaceIndex].node1,bdBox[cI].face[nearestFaceIndex].node2,
                            bdBox[cI].face[nearestFaceIndex].node3,uVec);
        //real *parPos = allocateDoubleArray(DIM);


        real gap = getOverlap(demPart[p->part_id].pos,demPart[p->part_id].dia, bdBox[cI].face[nearestFaceIndex].node1,  
                                bdBox[cI].face[nearestFaceIndex].node2, bdBox[cI].face[nearestFaceIndex].node3, uVec);
   
            //writeLog("logfile2.log","PAR TIME ",P_TIME(p));
        if(gap < 0) //If contact exists calculate contact force
        {
                //writeLogNum("logfile2.log", "PAR ", p->part_id);
                //
                //if(gap < -demPart[p->part_id].dia*0.05){gap = -demPart[p->part_id].dia*0.002;}
                //writeLogNum("logfile2.log", "gap ", 1e3*gap/lengthFactor);
            surfaceContactForce(p->part_id, -gap, uVec);
        } 
        free(uVec);
            //free(parPos);
            
    }
}

void checkXContact(Particle *p, real xMin, real xMax){
    //Contact with xMin
    real gap = demPart[p->part_id].pos[0] - xMin - demPart[p->part_id].dia*0.5;
    uVec[0] = 1.0;
    uVec[1] = 0.0;
    uVec[2] = 0.0;  
    
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p->part_id, -gap, uVec);
    }  
    if(gap < 100.e-9*lengthFactor){pWallVWForce(p->part_id, gap, uVec);}
    //Contact with xMax
    gap = xMax - demPart[p->part_id].pos[0] - demPart[p->part_id].dia*0.5;
    uVec[0] = -1.0;
    uVec[1] = 0.0;
    uVec[2] = 0.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p->part_id, -gap, uVec);
    }
    if(gap < 100.e-9*lengthFactor){pWallVWForce(p->part_id, gap, uVec);}
    //if(gap < demPart[p->part_id].dia*0.05){pWallVWForce(p->part_id, gap, uVec);}  
}

void checkZContact(Particle *p, real zMin, real zMax){
    //Contact with zMin
    real gap = -(zMin - demPart[p->part_id].pos[2]) - demPart[p->part_id].dia*0.5;
    uVec[0] = 0.0;
    uVec[1] = 0.0;
    uVec[2] = 1.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p->part_id, -gap, uVec);
    }  
    if(gap < 100e-9*lengthFactor){pWallVWForce(p->part_id, gap, uVec);}
    //Contact with z=88mm wall
    gap = zMax - demPart[p->part_id].pos[2] - demPart[p->part_id].dia*0.5;
    uVec[0] = 0.0;
    uVec[1] = 0.0;
    uVec[2] = -1.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p->part_id, -gap, uVec);
    } 
    if(gap < 100.e-9*lengthFactor){pWallVWForce(p->part_id, gap, uVec);} 
}

void checkYContact(Particle *p, real yMin, real yMax){
    // Contact with yMin
    real gap = -(yMin - demPart[p->part_id].pos[1]) - demPart[p->part_id].dia*0.5;
    uVec[0] = 0.0;
    uVec[1] = 1.0;
    uVec[2] = 0.0;  

    if(gap < 0) //If contact exists calculate contact force
    {
        //demPart[p->part_id].force[1] += -demPart[p->part_id].mass;
        surfaceContactForce(p->part_id, -gap, uVec);
    }
    if(gap < 100.e-9*lengthFactor){pWallVWForce(p->part_id, gap, uVec);}  

    // Contact with yMax
    gap = yMax - (demPart[p->part_id].pos[1] + demPart[p->part_id].dia*0.5);
    uVec[0] = 0.0;
    uVec[1] = -1.0;
    uVec[2] = 0.0;  

    if(gap < 0) //If contact exists calculate contact force
    {
        surfaceContactForce(p->part_id, -gap, uVec);
    }
    if(gap < 100.e-9*lengthFactor){pWallVWForce(p->part_id, gap, uVec);}
}





