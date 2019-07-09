#include "common.h"

/*Update particle position and insert them into cells*/
void updatePosition(Particle *p){
    demPart[p->part_id].insertable = 1;
    real dxDot = demPart[p->part_id].force[0]*timeStep/demPart[p->part_id].mass;
    real dyDot = demPart[p->part_id].force[1]*timeStep/demPart[p->part_id].mass;
    real dzDot = demPart[p->part_id].force[2]*timeStep/demPart[p->part_id].mass;
    
    demPart[p->part_id].vel[0] += dxDot;
    demPart[p->part_id].vel[1] += dyDot;
    demPart[p->part_id].vel[2] += dzDot;

    real dx = demPart[p->part_id].vel[0]*timeStep;
    real dy = demPart[p->part_id].vel[1]*timeStep;
    real dz = demPart[p->part_id].vel[2]*timeStep;

    demPart[p->part_id].pos[0] += dx;
    demPart[p->part_id].pos[1] += dy;
    demPart[p->part_id].pos[2] += dz;
   
    demPart[p->part_id].currentTime += timeStep;

    demPart[p->part_id].displacement += sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
    if(demPart[p->part_id].pos[0] > XMAX_BOUND*lengthFactor){
        demPart[p->part_id].active = 0;
    }

    //Update DPM particle position
    if(demPart[p->part_id].active == 1){
        if(updateDPM == 0){
            P_POS(p)[0] = demPart[p->part_id].pos[0]/lengthFactor;
            P_POS(p)[1] = demPart[p->part_id].pos[1]/lengthFactor; 
            P_POS(p)[2] = demPart[p->part_id].pos[2]/lengthFactor;
        }
        //Insert to cell
        //ceil gives upper bound, but we need lower bound. Therefore use -1 for index
        // int iIndex = ceil((demPart[p->part_id].pos[0]-xmin)/domainDx);
        // int jIndex = ceil((demPart[p->part_id].pos[1]-ymin)/domainDy);
        // int kIndex = ceil((demPart[p->part_id].pos[2]-zmin)/domainDz);
        // int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;  
        //if(demPart[p->part_id].active == 1){
       
    }
    // Thread *tc = P_CELL_THREAD(p);
    // cell_t c = P_CELL(p);
    // C_UDMI(c,tc,0) = 0.0;  
}



