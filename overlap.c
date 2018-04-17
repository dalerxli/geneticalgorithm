/***************************************
 * check overlap
 ***************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

int OverlapLattice(int n)
{
  int i,j;
  double rij2;
  double sigmaij,sigmaij2;
  double xij,yij,zij;
  double dxij,dyij,dzij;
  int check;
  int nx,ny,nz;
  VECTOR dR;

  int overlap;

//NCells = NCells+1;
  overlap = 0;

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
   // sigmaij = CoreIJ(identity[i],identity[j]);
// potential minimum distance
    sigmaij = SigmaIJ(identity[i],identity[j]);
    sigmaij2 = SQR(sigmaij);

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

    if(check != 1 && rij2 < sigmaij2)
    {
     overlap = 1;
//     NCells = NCells-1;
     return overlap;
    } //endif
  }//end loop j
}//end loop nx ny nz
}//end loop i

// NCells = NCells-1;

 return overlap; 

} // end function 

