/******************************************
 * headfile containing most global variables
 ******************************************/
#include <stdio.h>

#define MAX_NumberOfParticles 10000 // per unit cell
#define MAX_PopulationSize 100 // number of individuals
#define MAX_OffspringNumber 100 // number of offsprings
#define MAX_Dimensions 5

#define SQR(x) ((x)*(x))
#define CUBIC(x) ((x)*(x)*(x))
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define MIN(x,y) ((x)<(y) ? (x) : (y))

extern double randomseed;
extern int JobIndex;
extern int Dimension;

extern int InitialType,IDtype;

extern int NGeneration;
extern int OffspringNumber;
extern int PopulationSize; // number of individuals
extern int NumberOfParticles; //number of particles per unit cell
extern int NA,NB,NC; // NA + NB + NC = N
extern double fA,fB,fC;// fraction

extern double sigma[MAX_NumberOfParticles]; // hardcore diameter  sigma_i
extern double identity[MAX_NumberOfParticles]; // particle identity 1:A 2:B
extern double sigmaA,sigmaB,sigmaC,sigmaAB,sigmaBC,sigmaCA;

extern double rho,rhoA,rhoB,rhoC; // number density rho = N/V
extern double packingfraction; // phi = pi/6 rho in 3D
extern double V,L;  // Volume
extern int NumberOfLatticeSites; 
extern int NCells;
extern double mutationrate,similaritythreshold,distancemax,compressrate;
extern int RotationSwitch,CompressSwitch;

extern int PotentialType;
extern double n_rep,n_att,epsilonfactor; // replulsion and attraction index
extern double epsilonA,epsilonB,epsilonC,epsilonAB,epsilonBC,epsilonCA;
extern double coreA,coreB,coreC,coreAB,coreBC,coreCA;

extern double rc;// rcutoff
extern double rcA,rcB,rcC,rcAB,rcBC,rcCA;
extern double ucA,ucB,ucC,ucAB,ucBC,ucCA; // u(r=rc) at cutoff
extern double ducdrA,ducdrB,ducdrC,ducdrAB,ducdrBC,ducdrCA; // du(r=rc)/dr at cutoff

typedef struct
{
	double x;
	double y;
	double z;
} VECTOR;

typedef struct
{
	double xx;
	double xy;
	double xz;
	double yx;
	double yy;
	double yz;
	double zx;
	double zy;
	double zz;
} MATRIX;

extern MATRIX A;// rotation matrix

extern VECTOR position[3*MAX_PopulationSize][MAX_NumberOfParticles]; // position of particles
extern VECTOR displacement[3*MAX_PopulationSize][MAX_NumberOfParticles]; // displacement vectors
extern VECTOR projection[3*MAX_PopulationSize][MAX_Dimensions]; // lattice vectors

extern VECTOR bravais[3*MAX_PopulationSize][MAX_Dimensions]; // bravais lattice vectors

extern double Lx[3*MAX_PopulationSize],Ly[3*MAX_PopulationSize],Lz[3*MAX_PopulationSize],Anglex[3*MAX_PopulationSize],Angley[3*MAX_PopulationSize],Anglez[3*MAX_PopulationSize];

extern double fitness[3*MAX_PopulationSize]; // -energy
extern double fitnesshistory[100][3*MAX_PopulationSize]; // -energy Congen generations before
extern double volume[3*MAX_PopulationSize]; // volume of unit cell
extern double packfrac[3*MAX_PopulationSize]; // packing fraction

extern int ifrom[MAX_Dimensions];
extern int ito[MAX_Dimensions];

extern double ministep; // energy minimization step size
extern int MaxIter;// max iterations
extern int Congen; // number of generations beyond with no improvement then converge
extern double converge; // convergence criterion
extern VECTOR force[MAX_NumberOfParticles]; // force on particle i
extern VECTOR hstep[MAX_NumberOfParticles]; // displacement in conjugate gradient
extern double Upotential,Ftotal;

extern int globaln,globalk2;
extern VECTOR rotaxis;

extern VECTOR vertice[3*MAX_PopulationSize][8];

void ReadInput(void);
void Initialization(int job);
double BoxMuller(double mm, double ss);
void Writemovie(FILE *FilePtr,int n);

double Distance(int n,int i,int j); // return rij^2
double SigmaIJ(int iID,int jID);
double CoreIJ(int iID,int jID);
double Potential(int iID,int jID,double r); // u(rij)
double Potential_dr(int iID,int jID,double r); // du(rij)/dr
double LatticeArea(VECTOR r1,VECTOR r2,VECTOR r3);
double LatticeVolume(VECTOR r1,VECTOR r2,VECTOR r3);
VECTOR Equal(VECTOR r);
double Norm(VECTOR r);
VECTOR ScalarProduct(double lambda, VECTOR r);
VECTOR Add(VECTOR r1,VECTOR r2);
VECTOR Subtract(VECTOR r1,VECTOR r2);
double CrossNorm(VECTOR r1,VECTOR r2);
double DotProduct(VECTOR r1, VECTOR r2);
VECTOR CrossProduct(VECTOR r1,VECTOR r2);
double VectorAngle(VECTOR r1,VECTOR r2);

void BasisVector(VECTOR a,VECTOR b, VECTOR c, int n);
void LatticeReduction(int n);
void VectorReduction(int nn);

void Compress(int n);
void Minimize(int n);
void MinimizeEuler(int n);
void MinimizeEulerTwo(int n,int k2);
void Force(int n);
double Energy(double kappa);
void ForceLattice(int n);
double EnergyLattice(double kappa);
double EulerEnergyLattice(double angle);
double EulerEnergyLattice1(double angle);
double EulerEnergyLattice2(double angle);

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,double (*func)(double));
double brent(double ax, double bx, double cx, double (*f)(double), double tol,double *xmin);
int OverlapLattice(int n);
void GetVertice(int n);
void WriteVertice(FILE *FilePtr,int n);
void Writegene(FILE *FilePtr,int n);
void Store(int n,int m);
void LatticeShift(int n); // shift the first particle position to origin
void RotationVector(VECTOR rot,double theta,VECTOR *r);
void Rotation(int n,VECTOR rot,double theta);
void RotationXY(int n);
