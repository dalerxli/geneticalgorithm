/***************************************
 * calculate of force fx,fy,fz of 
 * each particle i at corrent position 
 * calculate potential energy Upotential
 ***************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void Force(int n)
{
  int i,j;
  double rcij,rcij2; // rcutoff
  double sigmaij;//,sigmaij2;
  double rij,rij2,rij6; // avoid sqrt() if possible
  double xij,yij,zij;
 // double dudr; // du/dr
  double fij;
  double fxij,fyij,fzij;
  double uij;

  for(i=0;i<NumberOfParticles;i++)
  {
    force[i].x = 0.;
    force[i].y = 0.;
    force[i].z = 0.;
  }

  Upotential = 0.;
  Ftotal = 0.; // NOT vector sum!

 for(i=0;i<NumberOfParticles-1;i++)
 {
  for(j=i+1;j<NumberOfParticles;j++)
  {
   // sigmaij = SigmaIJ(identity[i],identity[j]);
    
    sigmaij = CoreIJ(identity[i],identity[j]);
    //sigmaij2 = SQR(sigmaij2);

   // rcij = rc/sigmaA*sigmaij;
    rcij = rc/coreA*sigmaij;

    rcij2 = SQR(rcij);

    xij = position[n][i].x - position[n][j].x;
    yij = position[n][i].y - position[n][j].y;
    zij = position[n][i].z - position[n][j].z;
  
//    xij = xij - L*round(xij/L);
//    yij = yij - L*round(yij/L);
//    zij = zij - L*round(zij/L);
     
    rij2 = SQR(xij)+SQR(yij)+SQR(zij);

  //  if(rij2 < rcij2) // try avoid sqrt for r>rc
    if(1>0) // try avoid sqrt for r>rc
    {
      rij = sqrt(rij2); // avoid sqrt if possible
       //dudr =  Potential_dr(sigma[i],sigma[jID],rij);
      fij =  -Potential_dr(identity[i],identity[j],rij);

      uij =  Potential(identity[i],identity[j],rij);

      fxij = fij*xij/rij;
      fyij = fij*yij/rij;
      fzij = fij*zij/rij;

      force[i].x += fxij;
      force[i].y += fyij;
      force[i].z += fzij;

      force[j].x += -fxij;
      force[j].y += -fyij;
      force[j].z += -fzij;

      Upotential += uij;
      Ftotal     += fabs(fij); // not Vector Sum
    } //endif rij < rcutoff
  }//end loop j
 }//end loop i

 return;

} // end function Force()


double Energy(double kappa)
{
  int i,j;
  double rcij,rcij2; // rcutoff
  double rij,rij2,rij6; // avoid sqrt() if possible
  double sigmaij;//,sigmaij2;
  double xij,yij,zij;
  double dxij,dyij,dzij;
  double uij;
  
  double U;

  U = 0.;

 for(i=0;i<NumberOfParticles-1;i++)
 {
  for(j=i+1;j<NumberOfParticles;j++)
  {
   // sigmaij = SigmaIJ(identity[i],identity[j]);
    
    sigmaij = CoreIJ(identity[i],identity[j]);
    //sigmaij2 = SQR(sigmaij2);

   // rcij = rc/sigmaA*sigmaij;
    rcij = rc/coreA*sigmaij;

    rcij2 = SQR(rcij);


    xij = position[globaln][i].x - position[globaln][j].x;
    yij = position[globaln][i].y - position[globaln][j].y;
    zij = position[globaln][i].z - position[globaln][j].z;
    
    dxij = kappa*hstep[i].x - kappa*hstep[j].x;
    dyij = kappa*hstep[i].y - kappa*hstep[j].y;
    dzij = kappa*hstep[i].z - kappa*hstep[j].z;
   
    xij += dxij;
    yij += dyij;
    zij += dzij;
  
//    xij = xij - L*round(xij/L);
//    yij = yij - L*round(yij/L);
//    zij = zij - L*round(zij/L);
     
    rij2 = SQR(xij)+SQR(yij)+SQR(zij);

  //  if(rij2 < rcij2) // try avoid sqrt for r>rc
    if(1>0) // try avoid sqrt for r>rc
    {
      rij = sqrt(rij2); // avoid sqrt if possible

      uij =  Potential(identity[i],identity[j],rij);

      U += uij;
    } //endif rij < rcutoff
  }//end loop j
 }//end loop i

 return U; 

} // end function Energy()


void ForceLattice(int n)
{
  int i,j;
  double rcij,rcij2; // rcutoff
  double sigmaij;//,sigmaij2;
  double rij,rij2,rij6; // avoid sqrt() if possible
  double xij,yij,zij;
 // double dudr; // du/dr
  double fij;
  double fxij,fyij,fzij;
  double uij;

  int check;

  int nx,ny,nz;
  VECTOR dR;

  for(i=0;i<NumberOfParticles;i++)
  {
    force[i].x = 0.;
    force[i].y = 0.;
    force[i].z = 0.;
  }

  Upotential = 0.;
  Ftotal = 0.; // NOT vector sum!
 
for(i=0;i<NumberOfParticles;i++)
{
for(nx=-NCells;nx<=NCells;nx++)
for(ny=-NCells;ny<=NCells;ny++)
for(nz=-NCells;nz<=NCells;nz++)
{
  dR = Add(ScalarProduct(nx,bravais[n][0]),ScalarProduct(ny,bravais[n][1]));
  dR = Add(dR,ScalarProduct(nz,bravais[n][2]));

  for(j=0;j<NumberOfParticles;j++)
  {
   // sigmaij = SigmaIJ(identity[i],identity[j]);
    sigmaij = CoreIJ(identity[i],identity[j]);
    //sigmaij2 = SQR(sigmaij2);

   // rcij = rc/sigmaA*sigmaij;
    rcij = rc/coreA*sigmaij;

    rcij2 = SQR(rcij);

    xij = position[n][i].x - position[n][j].x;
    yij = position[n][i].y - position[n][j].y;
    zij = position[n][i].z - position[n][j].z;

    xij -= dR.x;
    yij -= dR.y;
    zij -= dR.z;
  
    rij2 = SQR(xij)+SQR(yij)+SQR(zij);

    check = 0;
    if(nx == 0 && ny ==0 && nz ==0 && j == i)
    check = 1;

    if(check != 1 && rij2 > 0.0 && rij2 < rcij2) // try avoid sqrt for r>rc
    {
      rij = sqrt(rij2); // avoid sqrt if possible

//if(rij2 < 1.0)
//printf("(i,j)=(%d,%d),(nx,ny,nz) = (%d %d %d) rij = %lf\n",i,j,nx,ny,nz,rij);
       //dudr =  Potential_dr(sigma[i],sigma[jID],rij);
      fij =  -Potential_dr(identity[i],identity[j],rij);

      uij =  Potential(identity[i],identity[j],rij);

      fxij = fij*xij/rij;
      fyij = fij*yij/rij;
      fzij = fij*zij/rij;

      force[i].x += fxij;
      force[i].y += fyij;
      force[i].z += fzij;

      force[j].x += -fxij;
      force[j].y += -fyij;
      force[j].z += -fzij;

      Upotential += uij;
      Ftotal     += fabs(fij); // not Vector Sum
    } //endif rij < rcutoff
  }//end loop j
}//end loop nx ny nz
}//end loop i


 Upotential /= 2.0;
 Ftotal /= 2.0;

 return;

} // end function ForceLattice()


double EnergyLattice(double kappa)
{
  int i,j;
  double rcij,rcij2; // rcutoff
  double rij,rij2,rij6; // avoid sqrt() if possible
  double sigmaij;//,sigmaij2;
  double xij,yij,zij;
  double dxij,dyij,dzij;
  double uij;
  
  double U;

  int check;
  
  int nx,ny,nz;
  VECTOR dR;

  U = 0.;

for(i=0;i<NumberOfParticles;i++)
{
for(nx=-NCells;nx<=NCells;nx++)
for(ny=-NCells;ny<=NCells;ny++)
for(nz=-NCells;nz<=NCells;nz++)
{
  dR = Add(ScalarProduct(nx,bravais[globaln][0]),ScalarProduct(ny,bravais[globaln][1]));
  dR = Add(dR,ScalarProduct(nz,bravais[globaln][2]));

  for(j=0;j<NumberOfParticles;j++)
  {
   // sigmaij = SigmaIJ(identity[i],identity[j]);
    sigmaij = CoreIJ(identity[i],identity[j]);
    //sigmaij2 = SQR(sigmaij2);

   // rcij = rc/sigmaA*sigmaij;
    rcij = rc/coreA*sigmaij;

    rcij2 = SQR(rcij);

    xij = position[globaln][i].x - position[globaln][j].x;
    yij = position[globaln][i].y - position[globaln][j].y;
    zij = position[globaln][i].z - position[globaln][j].z;
    
    xij -= dR.x;
    yij -= dR.y;
    zij -= dR.z;
    
    dxij = kappa*hstep[i].x - kappa*hstep[j].x;
    dyij = kappa*hstep[i].y - kappa*hstep[j].y;
    dzij = kappa*hstep[i].z - kappa*hstep[j].z;
   
    xij += dxij;
    yij += dyij;
    zij += dzij;
    
    rij2 = SQR(xij)+SQR(yij)+SQR(zij);
    
    check = 0;
    if(nx == 0 && ny ==0 && nz ==0 && j == i)
    check = 1;

    if(check != 1 && rij2 > 0.0 && rij2 < rcij2) // try avoid sqrt for r>rc
    {
      rij = sqrt(rij2); // avoid sqrt if possible

      uij =  Potential(identity[i],identity[j],rij);

      U += uij;
    } //endif rij < rcutoff
  }//end loop j
}//end loop nx ny nz
}//end loop i

 U /= 2.0;

 return U; 

} // end function EnergyLattice()


//rotate particles with respect to rotaxis with angle
//globaln and rotaxis need to be specified outside
double EulerEnergyLattice(double angle)
{
  int i,j;
  double rcij,rcij2; // rcutoff
  double rij,rij2,rij6; // avoid sqrt() if possible
  double sigmaij;//,sigmaij2;
  double xij,yij,zij;
  double uij;

  VECTOR positiontemp[MAX_NumberOfParticles]; // position of particles
  
  double U;

  int check;
  
  int nx,ny,nz;
  VECTOR dR;

for(i=0;i<NumberOfParticles;i++)
{
 positiontemp[i] = Equal(position[globaln][i]);
 RotationVector(rotaxis,angle,&(positiontemp[i]));
}

  U = 0.;

for(i=0;i<NumberOfParticles;i++)
{
for(nx=-NCells;nx<=NCells;nx++)
for(ny=-NCells;ny<=NCells;ny++)
for(nz=-NCells;nz<=NCells;nz++)
{
  dR = Add(ScalarProduct(nx,bravais[globaln][0]),ScalarProduct(ny,bravais[globaln][1]));
  dR = Add(dR,ScalarProduct(nz,bravais[globaln][2]));

  for(j=0;j<NumberOfParticles;j++)
  {
   // sigmaij = SigmaIJ(identity[i],identity[j]);
    sigmaij = CoreIJ(identity[i],identity[j]);
    //sigmaij2 = SQR(sigmaij2);

   // rcij = rc/sigmaA*sigmaij;
    rcij = rc/coreA*sigmaij;

    rcij2 = SQR(rcij);

    xij = positiontemp[i].x - positiontemp[j].x;
    yij = positiontemp[i].y - positiontemp[j].y;
    zij = positiontemp[i].z - positiontemp[j].z;
    
    xij -= dR.x;
    yij -= dR.y;
    zij -= dR.z;
    
    rij2 = SQR(xij)+SQR(yij)+SQR(zij);
    
    check = 0;
    if(nx == 0 && ny ==0 && nz ==0 && j == i)
    check = 1;

    if(check != 1 && rij2 > 0.0 && rij2 < rcij2) // try avoid sqrt for r>rc
    {
      rij = sqrt(rij2); // avoid sqrt if possible

      uij =  Potential(identity[i],identity[j],rij);

      U += uij;
    } //endif rij < rcutoff
  }//end loop j
}//end loop nx ny nz
}//end loop i

 U /= 2.0;

 return U; 

} // end function EulerEnergyLattice()


double EulerEnergyLattice1(double angle)
{
  int i,j;
  double rcij,rcij2; // rcutoff
  double rij,rij2,rij6; // avoid sqrt() if possible
  double sigmaij;//,sigmaij2;
  double xij,yij,zij;
  double uij;

  VECTOR positiontemp[MAX_NumberOfParticles]; // position of particles
  
  double U;

  int check;
  
  int nx,ny,nz;
  VECTOR dR;

for(i=0;i<NumberOfParticles;i++)
{
 positiontemp[i] = Equal(position[globaln][i]);
 //if(i>globalk2+1) // only rotate half of the chain
 if(i<=globalk2+1) // only rotate half of the chain
 RotationVector(rotaxis,angle,&(positiontemp[i]));
}

  U = 0.;

for(i=0;i<NumberOfParticles;i++)
{
for(nx=-NCells;nx<=NCells;nx++)
for(ny=-NCells;ny<=NCells;ny++)
for(nz=-NCells;nz<=NCells;nz++)
{
  dR = Add(ScalarProduct(nx,bravais[globaln][0]),ScalarProduct(ny,bravais[globaln][1]));
  dR = Add(dR,ScalarProduct(nz,bravais[globaln][2]));

  for(j=0;j<NumberOfParticles;j++)
  {
   // sigmaij = SigmaIJ(identity[i],identity[j]);
    sigmaij = CoreIJ(identity[i],identity[j]);
    //sigmaij2 = SQR(sigmaij2);

   // rcij = rc/sigmaA*sigmaij;
    rcij = rc/coreA*sigmaij;

    rcij2 = SQR(rcij);

    xij = positiontemp[i].x - positiontemp[j].x;
    yij = positiontemp[i].y - positiontemp[j].y;
    zij = positiontemp[i].z - positiontemp[j].z;
    
    xij -= dR.x;
    yij -= dR.y;
    zij -= dR.z;
    
    rij2 = SQR(xij)+SQR(yij)+SQR(zij);
    
    check = 0;
    if(nx == 0 && ny ==0 && nz ==0 && j == i)
    check = 1;

    if(check != 1 && rij2 > 0.0 && rij2 < rcij2) // try avoid sqrt for r>rc
    {
      rij = sqrt(rij2); // avoid sqrt if possible

      uij =  Potential(identity[i],identity[j],rij);

      U += uij;
    } //endif rij < rcutoff
  }//end loop j
}//end loop nx ny nz
}//end loop i

 U /= 2.0;

 return U; 

} // end function EulerEnergyLattice()


double EulerEnergyLattice2(double angle)
{
  int i,j;
  double rcij,rcij2; // rcutoff
  double rij,rij2,rij6; // avoid sqrt() if possible
  double sigmaij;//,sigmaij2;
  double xij,yij,zij;
  double uij;

  VECTOR positiontemp[MAX_NumberOfParticles]; // position of particles
  
  double U;

  int check;
  
  int nx,ny,nz;
  VECTOR dR;

for(i=0;i<NumberOfParticles;i++)
{
 positiontemp[i] = Equal(position[globaln][i]);
 //if(i<=globalk2+1) // only rotate half of the chain
 if(i>globalk2+1) // only rotate half of the chain
 RotationVector(rotaxis,angle,&(positiontemp[i]));
}

  U = 0.;

for(i=0;i<NumberOfParticles;i++)
{
for(nx=-NCells;nx<=NCells;nx++)
for(ny=-NCells;ny<=NCells;ny++)
for(nz=-NCells;nz<=NCells;nz++)
{
  dR = Add(ScalarProduct(nx,bravais[globaln][0]),ScalarProduct(ny,bravais[globaln][1]));
  dR = Add(dR,ScalarProduct(nz,bravais[globaln][2]));

  for(j=0;j<NumberOfParticles;j++)
  {
   // sigmaij = SigmaIJ(identity[i],identity[j]);
    sigmaij = CoreIJ(identity[i],identity[j]);
    //sigmaij2 = SQR(sigmaij2);

   // rcij = rc/sigmaA*sigmaij;
    rcij = rc/coreA*sigmaij;

    rcij2 = SQR(rcij);

    xij = positiontemp[i].x - positiontemp[j].x;
    yij = positiontemp[i].y - positiontemp[j].y;
    zij = positiontemp[i].z - positiontemp[j].z;
    
    xij -= dR.x;
    yij -= dR.y;
    zij -= dR.z;
    
    rij2 = SQR(xij)+SQR(yij)+SQR(zij);
    
    check = 0;
    if(nx == 0 && ny ==0 && nz ==0 && j == i)
    check = 1;

    if(check != 1 && rij2 > 0.0 && rij2 < rcij2) // try avoid sqrt for r>rc
    {
      rij = sqrt(rij2); // avoid sqrt if possible

      uij =  Potential(identity[i],identity[j],rij);

      U += uij;
    } //endif rij < rcutoff
  }//end loop j
}//end loop nx ny nz
}//end loop i

 U /= 2.0;

 return U; 

} // end function EulerEnergyLattice()

