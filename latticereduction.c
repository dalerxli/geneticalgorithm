/***************************************************
* find Bravais lattice vector from projection vector
* updata displacement vectors and projection vectors
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

void LatticeShift(int n) // shift the first particle position to origin
{
 int i;

 position[n][0].x = 0.; 
 position[n][0].y = 0.; 
 position[n][0].z = 0.; 
 for(i=1;i<NumberOfParticles;i++)
 position[n][i] = Add(position[n][i-1],displacement[n][i-1]);
 
 return;
}

void VectorReduction(int n) // reduce the displacement vectors
{
 int nx,ny,nz;
 int i,d;

 double norm;

 VECTOR dR,dr;

for(i=0;i<NumberOfParticles-1;i++)
{
 displacement[n][i] = Subtract(position[n][i+1],position[n][i]);

 norm =  Norm(displacement[n][i]);

 for(nx=-NCells;nx<=NCells;nx++)
 for(ny=-NCells;ny<=NCells;ny++)
 for(nz=-NCells;nz<=NCells;nz++)
 {
  dR = Add(ScalarProduct(nx,bravais[n][0]),ScalarProduct(ny,bravais[n][1]));
  dR = Add(dR,ScalarProduct(nz,bravais[n][2]));

  dr = Subtract(Add(position[n][i+1],dR),position[n][i]);

  if(Norm(dr) < norm)
  {
   displacement[n][i] = Equal(dr);
   norm = Norm(dr);
  }
 }
 position[n][i+1] = Add(position[n][i],displacement[n][i]);
}

//useful here?
 LatticeShift(n); // shift 1st particle to origin

for(d=0;d<Dimension;d++)
{
 projection[n][d] = Equal(bravais[n][d]);

 for(i=ito[d];i<ifrom[d];i++)
 projection[n][d] = Subtract(projection[n][d],displacement[n][i]);
}
 
 return;
}

//void LatticeReduction(void)
void LatticeReduction(int n)
{ 
int d,i;

VECTOR rb[MAX_Dimensions]; // lattice vectors

double slattice,stemp; // lattice area

//for(n=0;n<PopulationSize;n++)
//{
// rb[0] = Equal(projection[n][0]);
// for(i=ito[0];i<ifrom[0];i++)
// rb[0] = Add(rb[0],displacement[n][i]);
 
// rb[1] = Equal(projection[n][1]);
// for(i=ito[1];i<ifrom[1];i++)
// rb[1] = Add(rb[1],displacement[n][i]);

// rb[2] = Equal(projection[n][2]);
// for(i=ito[2];i<ifrom[2];i++)
// rb[2] = Add(rb[2],displacement[n][i]);

 rb[0] = Equal(bravais[n][0]);
 rb[1] = Equal(bravais[n][1]);
 rb[2] = Equal(bravais[n][2]);
 
 stemp = LatticeArea(rb[0],rb[1],rb[2]); 
 slattice = stemp;
 bravais[n][0] = Equal(rb[0]);
 bravais[n][1] = Equal(rb[1]);
 bravais[n][2] = Equal(rb[2]);
 
 stemp = LatticeArea(Add(rb[0],rb[1]),rb[1],rb[2]); 
 if(stemp < slattice)
 {
  bravais[n][0] = Add(rb[0],rb[1]);
  bravais[n][1] = Equal(rb[1]);
  bravais[n][2] = Equal(rb[2]);
  slattice = stemp;
 }
 
 stemp = LatticeArea(Subtract(rb[0],rb[1]),rb[1],rb[2]); 
 if(stemp < slattice)
 {
  bravais[n][0] = Subtract(rb[0],rb[1]);
  bravais[n][1] = Equal(rb[1]);
  bravais[n][2] = Equal(rb[2]);
  slattice = stemp;
 }
 
 stemp = LatticeArea(Add(rb[0],rb[2]),rb[1],rb[2]); 
 if(stemp < slattice)
 {
  bravais[n][0] = Add(rb[0],rb[2]);
  bravais[n][1] = Equal(rb[1]);
  bravais[n][2] = Equal(rb[2]);
  slattice = stemp;
 }
 
 stemp = LatticeArea(Subtract(rb[0],rb[2]),rb[1],rb[2]); 
 if(stemp < slattice)
 {
  bravais[n][0] = Subtract(rb[0],rb[2]);
  bravais[n][1] = Equal(rb[1]);
  bravais[n][2] = Equal(rb[2]);
  slattice = stemp;
 }
 
 stemp = LatticeArea(rb[0],Add(rb[1],rb[0]),rb[2]); 
 if(stemp < slattice)
 {
  bravais[n][0] = Equal(rb[0]);
  bravais[n][1] = Add(rb[1],rb[0]);
  bravais[n][2] = Equal(rb[2]);
  slattice = stemp;
 }
 
 stemp = LatticeArea(rb[0],Subtract(rb[1],rb[0]),rb[2]); 
 if(stemp < slattice)
 {
  bravais[n][0] = Equal(rb[0]);
  bravais[n][1] = Subtract(rb[1],rb[0]);
  bravais[n][2] = Equal(rb[2]);
  slattice = stemp;
 }
 
 stemp = LatticeArea(rb[0],Add(rb[1],rb[2]),rb[2]); 
 if(stemp < slattice)
 {
  bravais[n][0] = Equal(rb[0]);
  bravais[n][1] = Add(rb[1],rb[2]);
  bravais[n][2] = Equal(rb[2]);
  slattice = stemp;
 }
 
 stemp = LatticeArea(rb[0],Subtract(rb[1],rb[2]),rb[2]); 
 if(stemp < slattice)
 {
  bravais[n][0] = Equal(rb[0]);
  bravais[n][1] = Subtract(rb[1],rb[2]);
  bravais[n][2] = Equal(rb[2]);
  slattice = stemp;
 }
 
 stemp = LatticeArea(rb[0],rb[1],Add(rb[2],rb[0])); 
 if(stemp < slattice)
 {
  bravais[n][0] = Equal(rb[0]);
  bravais[n][1] = Equal(rb[1]);
  bravais[n][2] = Add(rb[2],rb[0]);
  slattice = stemp;
 }
 
 stemp = LatticeArea(rb[0],rb[1],Subtract(rb[2],rb[0])); 
 if(stemp < slattice)
 {
  bravais[n][0] = Equal(rb[0]);
  bravais[n][1] = Equal(rb[1]);
  bravais[n][2] = Subtract(rb[2],rb[0]);
  slattice = stemp;
 }

 stemp = LatticeArea(rb[0],rb[1],Add(rb[2],rb[1])); 
 if(stemp < slattice)
 {
  bravais[n][0] = Equal(rb[0]);
  bravais[n][1] = Equal(rb[1]);
  bravais[n][2] = Add(rb[2],rb[1]);
  slattice = stemp;
 }
 
 stemp = LatticeArea(rb[0],rb[1],Subtract(rb[2],rb[1])); 
 if(stemp < slattice)
 {
  bravais[n][0] = Equal(rb[0]);
  bravais[n][1] = Equal(rb[1]);
  bravais[n][2] = Subtract(rb[2],rb[1]);
  slattice = stemp;
 }

 BasisVector(bravais[n][0],bravais[n][1],bravais[n][2],n);

 VectorReduction(n);

// Rotation(n);

//}//end loop over individuals
 
return;
}
