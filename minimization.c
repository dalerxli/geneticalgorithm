/***************************************************
* energy minimization
* conjugate-gradient
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

void Minimize(int n)
{
 int iter; // iteration
 int i,d;
 double ffnew,ffold,gamma;

// VECTOR position_old[MAX_NumberOfParticles];
 VECTOR force_old[MAX_NumberOfParticles];
 double Uold;

 double ax,bx,cx,fa,fb,fc;
 double tol,kappa;

 tol = 1.0e-10;

 globaln = n;

 ForceLattice(n); // force[i]
 Uold = Upotential;
 fitness[n] = -Upotential;

 for(i=0;i<NumberOfParticles;i++)
 {
  force_old[i] = Equal(force[i]);
  hstep[i] = Equal(force[i]);
 }

 for(iter=0;iter<MaxIter;iter++)
 {
 ax = 0.0;
 bx = 1.0;

  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,&EnergyLattice);
  brent(ax,bx,cx,&EnergyLattice,tol,&kappa);
  
 // mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,&Energy);
 // brent(ax,bx,cx,&Energy,tol,&kappa);

  for(i=0;i<NumberOfParticles;i++)
   position[n][i] = Add(position[n][i],ScalarProduct(kappa,hstep[i])); 
  
  ForceLattice(n);
  
  if(Uold - Upotential < converge)
  {
/////////////////////////////////////
 //  fitness[n] = -Upotential;
 //  LatticeShift(n);
//   return;
/////////////////////////////////////
   break;
  }
  
  ffold = 0.;
  ffnew = 0.;
  for(i=0;i<NumberOfParticles;i++)
  {
   ffold += DotProduct(force_old[i],force_old[i]);
 //  ffnew += DotProduct(force[i],force[i]); // Fletcher-Reeves 100 times slower!!
   ffnew += DotProduct(Subtract(force[i],force_old[i]),force[i]); // Polak-Ribiere
  }
  gamma = ffnew/ffold;
  
  Uold = Upotential;

  for(i=0;i<NumberOfParticles;i++)
   force_old[i] = Equal(force[i]);
  
  for(i=0;i<NumberOfParticles;i++)
   hstep[i] = Add(force[i],ScalarProduct(gamma,hstep[i]));
  

 }//end loop over iteration

/*****************************************/
//for(n=0;n<PopulationSize;n++)
 for(i=0;i<NumberOfParticles-1;i++)
  displacement[n][i] = Subtract(position[n][i+1],position[n][i]); //Bi is vector from i to i+1

//for(n=0;n<PopulationSize;n++)
for(d=0;d<Dimension;d++)
{
 projection[n][d] = Equal(bravais[n][d]);

 for(i=ito[d];i<ifrom[d];i++)
 projection[n][d] = Subtract(projection[n][d],displacement[n][i]);
}
 
 fitness[n] = -Upotential;

//shift first particle to origin
 LatticeShift(n);
/*****************************************/

 return;
}

void MinimizeEuler(int n)
{
 double ax,bx,cx,fa,fb,fc;
 double angle;
 double tol;
 VECTOR rotaxis2;

 int i,d;

 tol = 1.0e-10;

 globaln = n;
 
 rotaxis.x = 0.;
 rotaxis.y = 0.;
 rotaxis.z = 1.;
 ax = 0.0;
 bx = 1.0;
 mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,&EulerEnergyLattice);
 brent(ax,bx,cx,&EulerEnergyLattice,tol,&angle);
 for(i=0;i<NumberOfParticles;i++)
 RotationVector(rotaxis,angle,&(position[n][i]));

 rotaxis.x = cos(angle);
 rotaxis.y = sin(angle);
 rotaxis.z = 0.;
 ax = 0.0;
 bx = 1.0;
 mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,&EulerEnergyLattice);
 brent(ax,bx,cx,&EulerEnergyLattice,tol,&angle);
 for(i=0;i<NumberOfParticles;i++)
 RotationVector(rotaxis,angle,&(position[n][i]));

 rotaxis2.x = 0.;
 rotaxis2.y = 0.;
 rotaxis2.z = 1.;
 RotationVector(rotaxis,angle,&(rotaxis2));
 rotaxis = Equal(rotaxis2);
 ax = 0.0;
 bx = 1.0;
 mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,&EulerEnergyLattice);
 brent(ax,bx,cx,&EulerEnergyLattice,tol,&angle);
 for(i=0;i<NumberOfParticles;i++)
 RotationVector(rotaxis,angle,&(position[n][i]));

 for(i=0;i<NumberOfParticles-1;i++)
  displacement[n][i] = Subtract(position[n][i+1],position[n][i]); //Bi is vector from i to i+1

 for(d=0;d<Dimension;d++)
 {
  projection[n][d] = Equal(bravais[n][d]);

  for(i=ito[d];i<ifrom[d];i++)
  projection[n][d] = Subtract(projection[n][d],displacement[n][i]);
 }
 
 LatticeShift(n);

 return;
}

//rotate separately
void MinimizeEulerTwo(int n,int k2)
{
 double ax,bx,cx,fa,fb,fc;
 double angle;
 double tol;
 VECTOR rotaxis2;

 int i,d;

 tol = 1.0e-10;

 globaln = n;
 globalk2 = k2; // joint point is k2+1
 
 rotaxis.x = 0.;
 rotaxis.y = 0.;
 rotaxis.z = 1.;
 ax = 0.0;
 bx = 1.0;
 mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,&EulerEnergyLattice1);
 brent(ax,bx,cx,&EulerEnergyLattice1,tol,&angle);
 for(i=0;i<NumberOfParticles;i++)
 if(i<=k2+1) RotationVector(rotaxis,angle,&(position[n][i]));

 rotaxis.x = cos(angle);
 rotaxis.y = sin(angle);
 rotaxis.z = 0.;
 ax = 0.0;
 bx = 1.0;
 mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,&EulerEnergyLattice1);
 brent(ax,bx,cx,&EulerEnergyLattice1,tol,&angle);
 for(i=0;i<NumberOfParticles;i++)
 if(i<=k2+1) RotationVector(rotaxis,angle,&(position[n][i]));

 rotaxis2.x = 0.;
 rotaxis2.y = 0.;
 rotaxis2.z = 1.;
 RotationVector(rotaxis,angle,&(rotaxis2));
 rotaxis = Equal(rotaxis2);
 ax = 0.0;
 bx = 1.0;
 mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,&EulerEnergyLattice1);
 brent(ax,bx,cx,&EulerEnergyLattice1,tol,&angle);
 for(i=0;i<NumberOfParticles;i++)
 if(i<=k2+1) RotationVector(rotaxis,angle,&(position[n][i]));
 
 rotaxis.x = 0.;
 rotaxis.y = 0.;
 rotaxis.z = 1.;
 ax = 0.0;
 bx = 1.0;
 mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,&EulerEnergyLattice2);
 brent(ax,bx,cx,&EulerEnergyLattice2,tol,&angle);
 for(i=0;i<NumberOfParticles;i++)
 if(i>k2+1) RotationVector(rotaxis,angle,&(position[n][i]));

 rotaxis.x = cos(angle);
 rotaxis.y = sin(angle);
 rotaxis.z = 0.;
 ax = 0.0;
 bx = 1.0;
 mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,&EulerEnergyLattice2);
 brent(ax,bx,cx,&EulerEnergyLattice2,tol,&angle);
 for(i=0;i<NumberOfParticles;i++)
 if(i>k2+1) RotationVector(rotaxis,angle,&(position[n][i]));

 rotaxis2.x = 0.;
 rotaxis2.y = 0.;
 rotaxis2.z = 1.;
 RotationVector(rotaxis,angle,&(rotaxis2));
 rotaxis = Equal(rotaxis2);
 ax = 0.0;
 bx = 1.0;
 mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,&EulerEnergyLattice2);
 brent(ax,bx,cx,&EulerEnergyLattice2,tol,&angle);
 for(i=0;i<NumberOfParticles;i++)
 if(i>k2+1) RotationVector(rotaxis,angle,&(position[n][i]));

 for(i=0;i<NumberOfParticles-1;i++)
  displacement[n][i] = Subtract(position[n][i+1],position[n][i]); //Bi is vector from i to i+1

 for(d=0;d<Dimension;d++)
 {
  projection[n][d] = Equal(bravais[n][d]);

  for(i=ito[d];i<ifrom[d];i++)
  projection[n][d] = Subtract(projection[n][d],displacement[n][i]);
 }
 
 LatticeShift(n);

 return;
}

void Compress(int n)
{
 int i,d,dd;
 double norm;
 int overlap;

 dd=0;
 norm=0.;
 for(d=0;d<Dimension;d++)
 {
  if(Norm(bravais[n][d]) > norm)
  {
   norm = Norm(bravais[n][d]);
   dd = d;
  }
 }

 do
 {
  bravais[n][dd] = ScalarProduct(compressrate,bravais[n][dd]);
 }
 while(OverlapLattice(n) == 0);
  
 bravais[n][dd] = ScalarProduct(1./compressrate,bravais[n][dd]);
 
 projection[n][dd] = Equal(bravais[n][dd]);
 for(i=ito[dd];i<ifrom[dd];i++)
 projection[n][dd] = Subtract(projection[n][dd],displacement[n][i]);

 ForceLattice(n);
 fitness[n] = -Upotential;
 
 return;
}
