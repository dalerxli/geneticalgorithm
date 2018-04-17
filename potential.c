/***************************************************
* return distance r, u(r), du/dr, u_tail, w_tail
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

double Distance(int n,int i,int j) // return rij^2
{
  double r2,dx,dy,dz;
  
  // miminum image convention under periodic boundary condition
  dx = position[n][i].x-position[n][j].x;
// dx = dx - L*round(dx/L);

  dy = position[n][i].y-position[n][j].y;
// dy = dy - L*round(dy/L);

  dz = position[n][i].z-position[n][j].z;
// dz = dz - L*round(dz/L);
 
  r2 = SQR(dx)+SQR(dy)+SQR(dz); // reduced unit
  //r = sqrt(SQR(dx)+SQR(dy)+SQR(dz));

  return r2;
}

double SigmaIJ(int iID,int jID)
{
 double sigmaij;

 if(iID == jID)
 {
  if(iID == 1) // A-A
  sigmaij = sigmaA;
  else if(iID == 2) // B-N
  sigmaij = sigmaB;
  else
  sigmaij = sigmaC;
 }
 else
 {
  if((iID == 1 && jID == 2) || (iID == 2 && jID == 1))
   sigmaij = sigmaAB;
  else if((iID == 2 && jID == 3) || (iID == 3 && jID == 2))
   sigmaij = sigmaBC;
  else
   sigmaij = sigmaCA;
 }

 return sigmaij;
}

double CoreIJ(int iID,int jID)
{
 double coreij;

 if(iID == jID)
 {
  if(iID == 1) // A-A
  coreij = coreA;
  else if(iID == 2) // B-N
  coreij = coreB;
  else
  coreij = coreC;
 }
 else
 {
  if((iID == 1 && jID == 2) || (iID == 2 && jID == 1))
   coreij = coreAB;
  else if((iID == 2 && jID == 3) || (iID == 3 && jID == 2))
   coreij = coreBC;
  else
   coreij = coreCA;
 }

 return coreij;
}

// A: ID = 1 B:ID = 2 C:ID = 3
double Potential(int iID,int jID,double r) // u(rij)
{
  double u;
 
  if(PotentialType == 0) // L-J
  {
    if(iID == jID)
    {
     if(iID == 1) // A-A
      if(r < rcA)
      u = epsilonfactor*epsilonA*(pow(coreA/r,n_rep)-pow(coreA/r,n_att));
      else
      u = 0.;
     else if(iID == 2) // B-B
      if(r < rcB)
      u = epsilonfactor*epsilonB*(pow(coreB/r,n_rep)-pow(coreB/r,n_att));
      else
      u = 0.;
     else
      if(r < rcC) // C-C
      u = epsilonfactor*epsilonC*(pow(coreC/r,n_rep)-pow(coreC/r,n_att));
      else
      u = 0.;
    }
    else  
    {
     if((iID == 1 && jID == 2) || (iID == 2 && jID == 1)) // A-B
      if(r < rcAB)
      u = epsilonfactor*epsilonAB*(pow(coreAB/r,n_rep)-pow(coreAB/r,n_att));
      else
      u = 0.;
     else if((iID == 2 && jID == 3) || (iID == 3 && jID == 2)) // B-C
      if(r < rcBC)
      u = epsilonfactor*epsilonBC*(pow(coreBC/r,n_rep)-pow(coreBC/r,n_att));
      else
      u = 0.;
     else // C-A
      if(r < rcCA)
      u = epsilonfactor*epsilonCA*(pow(coreCA/r,n_rep)-pow(coreCA/r,n_att));
      else
      u = 0.;
    }
  } //endif L-J

  return u;
}

double Potential_dr(int iID,int jID,double r) // du(rij)/dr
{
  double dudr;
 
  if(PotentialType == 0) // L-J m-n
  {
    if(iID == jID)
    {
     if(iID == 1) // A-A
      if(r < rcA)
      dudr = -epsilonfactor*n_rep*epsilonA/r*(pow(coreA/r,n_rep)-pow(coreA/r,n_att)*n_att/n_rep);
      else
      dudr = 0.;
     else if(iID == 2) // B-B
      if(r < rcB)
      dudr = -epsilonfactor*n_rep*epsilonB/r*(pow(coreB/r,n_rep)-pow(coreB/r,n_att)*n_att/n_rep);
      else
      dudr = 0.;
     else
      if(r < rcC)
      dudr = -epsilonfactor*n_rep*epsilonC/r*(pow(coreC/r,n_rep)-pow(coreC/r,n_att)*n_att/n_rep);
      else
      dudr = 0.;
    }
    else // A-B
    {
     if((iID == 1 && jID == 2) || (iID == 2 && jID == 1))
      if(r < rcAB)
      dudr = -epsilonfactor*n_rep*epsilonAB/r*(pow(coreAB/r,n_rep)-pow(coreAB/r,n_att)*n_att/n_rep);
      else
      dudr = 0.;
     else if((iID == 2 && jID == 3) || (iID == 3 && jID == 2))
      if(r < rcBC)
      dudr = -epsilonfactor*n_rep*epsilonBC/r*(pow(coreBC/r,n_rep)-pow(coreBC/r,n_att)*n_att/n_rep);
      else
      dudr = 0.;
     else
      if(r < rcCA)
      dudr = -epsilonfactor*n_rep*epsilonCA/r*(pow(coreCA/r,n_rep)-pow(coreCA/r,n_att)*n_att/n_rep);
      else
      dudr = 0.;
    }
  } //endif L-J

  return dudr;
}
