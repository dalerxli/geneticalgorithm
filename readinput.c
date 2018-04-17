/*********************************************************************
 * input simulation parameters from file "input"
 *********************************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

void ReadInput(void)
{
 int i;
 FILE *fp;
 
 fp=fopen("input","r"); 
 fscanf(fp,"%d %d",&InitialType,&IDtype);
 fscanf(fp,"%d %d %d",&NGeneration, &PopulationSize,&OffspringNumber);
 fscanf(fp,"%d %d %d %d",&NumberOfParticles,&NA,&NB,&NC);
 fscanf(fp,"%lf %lf %lf",&sigmaA,&sigmaB,&sigmaC);
 fscanf(fp,"%lf %lf %lf",&sigmaAB,&sigmaBC,&sigmaCA);
 fscanf(fp,"%lf",&randomseed);
 fscanf(fp,"%d",&PotentialType);
 fscanf(fp,"%lf %lf",&n_rep,&n_att);
 fscanf(fp,"%lf %lf %lf",&epsilonA,&epsilonB,&epsilonC);
 fscanf(fp,"%lf %lf %lf",&epsilonAB,&epsilonBC,&epsilonCA);
 fscanf(fp,"%lf",&rc);
 fscanf(fp,"%lf %lf %d",&ministep,&converge,&MaxIter);
 fscanf(fp,"%lf %lf %lf %lf %d",&mutationrate,&similaritythreshold,&distancemax,&compressrate,&NCells);
 fscanf(fp,"%d %d",&RotationSwitch,&CompressSwitch);
 fclose(fp);
 
 Dimension = 3;

 distancemax = sigmaA*NumberOfParticles*distancemax;

 //compressrate = 0.99;

 if(PotentialType == 0)
 {
  // make the mimium of the potential coincides with hard sphere diameter
  coreA = (sigmaA)*pow(n_rep/n_att,1./(n_att-n_rep));
  coreB = (sigmaB)*pow(n_rep/n_att,1./(n_att-n_rep));
  coreC = (sigmaC)*pow(n_rep/n_att,1./(n_att-n_rep));
  coreAB = (sigmaAB)*pow(n_rep/n_att,1./(n_att-n_rep));
  coreBC = (sigmaBC)*pow(n_rep/n_att,1./(n_att-n_rep));
  coreCA = (sigmaCA)*pow(n_rep/n_att,1./(n_att-n_rep));
  epsilonfactor = 1.0/(pow(n_rep/n_att,n_att/(n_att-n_rep)) - pow(n_rep/n_att,n_rep/(n_att-n_rep)));
 }
 
 rcA = rc*coreA;
 rcB = rc*coreB;
 rcC = rc*coreC;
 rcAB = rc*coreAB;
 rcBC = rc*coreBC;
 rcCA = rc*coreCA;
 
// TimeBig = 1.0E10;
// tolerance = 1.0E-10;
// timetol = 1.0E-9;

 fA = 1.*NA/NumberOfParticles;
 fB = 1.*NB/NumberOfParticles;
 fC = 1.*NC/NumberOfParticles;

 //rho = NumberOfParticles/V*CUBIC(sigmaA);
 //rhoA = NA/V*CUBIC(sigmaA);
 //rhoB = NB/V*CUBIC(sigmaA);
 //rhoC = NC/V*CUBIC(sigmaA);
 //L = pow(V,1.0/3.0); 

 printf("\n");
 printf("population size = %d\n",PopulationSize);
 printf("offspring number = %d\n",OffspringNumber);

 printf("\n");
 printf("N = %d\tNA:NB:NC = %d:%d:%d\tfA:fB:fC = %lf:%lf:%lf\n",NumberOfParticles,NA,NB,NC,fA,fB,fC);
 printf("sigmaA = %lf\tsigmaB = %lf\tsigmaC = %lf\n",sigmaA,sigmaB,sigmaC);
 printf("sigmaAB = %lf\tsigmaBC = %lf\tsigmaCA = %lf\n",sigmaAB,sigmaBC,sigmaCA);
 printf("\n");
 
 printf("coreA = %lf\tcoreB = %lf\tcoreC = %lf\n",coreA,coreB,coreC);
 printf("m = %lf\tn = %lf\tepsilon factor = %lf\n",n_rep,n_att,epsilonfactor);
 printf("A rm = %lf\tum = %lf\tf0 = %lf\n",pow(n_rep/n_att,1.0/(n_rep-n_att))*coreA,Potential(1,1,pow(n_rep/n_att,1.0/(n_rep-n_att))*coreA),-Potential_dr(1,1,coreA));
 printf("B rm = %lf\tum = %lf\tf0 = %lf\n",pow(n_rep/n_att,1.0/(n_rep-n_att))*coreB,Potential(2,2,pow(n_rep/n_att,1.0/(n_rep-n_att))*coreB),-Potential_dr(1,1,coreB));
 printf("C rm = %lf\tum = %lf\tf0 = %lf\n",pow(n_rep/n_att,1.0/(n_rep-n_att))*coreC,Potential(3,3,pow(n_rep/n_att,1.0/(n_rep-n_att))*coreC),-Potential_dr(1,1,coreC));
 printf("AB rm = %lf\tum = %lf\tf0 = %lf\n",pow(n_rep/n_att,1.0/(n_rep-n_att))*coreAB,Potential(1,2,pow(n_rep/n_att,1.0/(n_rep-n_att))*coreAB),-Potential_dr(1,1,coreAB));
 printf("BC rm = %lf\tum = %lf\tf0 = %lf\n",pow(n_rep/n_att,1.0/(n_rep-n_att))*coreBC,Potential(3,2,pow(n_rep/n_att,1.0/(n_rep-n_att))*coreBC),-Potential_dr(1,1,coreBC));
 printf("CA rm = %lf\tum = %lf\tf0 = %lf\n",pow(n_rep/n_att,1.0/(n_rep-n_att))*coreCA,Potential(3,1,pow(n_rep/n_att,1.0/(n_rep-n_att))*coreCA),-Potential_dr(1,1,coreCA));

 printf("max iterations in conjugate gradient %d\n",MaxIter);
 printf("mutationrate = %lf\n",mutationrate);
 printf("similarity threshold = %lf\n",similaritythreshold);
 printf("maximum distance between particles = %lf\n",distancemax);
 printf("compress rate = %lf\n",compressrate);
 printf("number of unit cells per dimension = %d\n",NCells);
 return;
}
