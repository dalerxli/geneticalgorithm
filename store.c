/***************************************************
* store the lattice information of n to m
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

void Store(int n,int m)
{
 int i,d;

 for(i=0;i<NumberOfParticles;i++)
 {
  position[m][i] = Equal(position[n][i]);
  if(i<NumberOfParticles-1)
  displacement[m][i] = Equal(displacement[n][i]);
 }

 for(d=0;d<Dimension;d++)
 {
  projection[m][d] = Equal(projection[n][d]);
  bravais[m][d] = Equal(bravais[n][d]);
 }

 return;
}
