#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<string.h>

static double PI = 3.1415926535898;

int N = 2;
int z1 = 1;		//H
int z2 = 1;	//H
int max_iterations;

double convergence_threshold;
double a;
double b;
double rab;
double rab2;
double zeta1 = 1.24;		//H
double zeta2 = 1.24;	//H
double d1;
double d2;
double N_A;
double N_B;
double rp;
double rq;
double rpa;
double rpb;
double raq;
double rbq;
double rpq;
double raq2;
double rbq2;
double rpq2;
double rpb;
double rpb2;
double t;
double V1111;
double V1112;
double V1122;
double V1222;
double V1212;
double V2222;
double R = 1.4;
double R2;
double E_total;

double coefs1[3] = {0.444635,0.535328,0.154329};
double alphas1[3] = {0.109818,0.405771,2.22766};
double coefs2[3] = {0.444635,0.535328,0.154329};
double alphas2[3] = {0.109818,0.405771,2.22766};
double r[2];

double **S;				//overlap matrix
double **T;				//kinetic energy
double **V1;			//external potential from atom one
double **V2;			//external potential from atom two
double **H_core;	//one body part of the hamiltonian
double **G;				//mean field electron-electron interactions
double **X;				//canonical orthogonalization matrix
double **Xt;			//transpose of X
double **P;				//population matrix
double **P_old;		//prior population matrix, used to check for convergence
double **P_prime;	//population matrix, diagonalized
double **F;				//Fock matrix
double **F_prime;	//Fock matrix in terms of orthogonalized orbitals
double **F_temp;	//intermediate matrix used to transform between F and F'
double **E;				//orbital energies
double **C;				//
double **C_prime;	//C in terms of orthogonalized orbitals
double ****mat2e;	//two electron integral matrix

double F0(double arg) {
	if (arg < 1e-6) {
		return 1-arg/3;
	}
	else {
		return 0.5*sqrt(PI/arg)*erf(sqrt(arg));
	}
}

double integral(double a,double b,double c,double d,double rab2,double rcd2, double rpq2) {
	return 2*pow(PI,2.5)/((a+b)*(c+d)*sqrt(a+b+c+d))
	*exp(-a*b/(a+b)*rab2-c*d/(c+d)*rcd2)
	*F0((a+b)*(c+d)/(a+b+c+d)*rpq2);
}

void make_G()	{
	for(int i=0;i<N;i++)	{
		for(int j=0;j<N;j++)	{
			G[i][j] = 0;
			for(int k=0;k<N;k++)	{
				for(int l=0;l<N;l++)	{
					G[i][j] += P[k][l] * (mat2e[i][j][k][l] - 0.5*mat2e[i][j][l][k]);
				}
			}
		}
	}
}

void make_F()	{
	for(int i=0;i<N;i++)	{
		for(int j=0;j<N;j++)	{
			F[i][j] = H_core[i][j] + G[i][j];
		}
	}
}

void mult(double **A, double **B, double **C)	{
	for(int i=0;i<N;i++)	{
		for(int j=0;j<N;j++)	{
			C[i][j] = 0;
			for(int k=0;k<N;k++)	{
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

// A = X, B = F, C = F_prime
void transform_XF(double **A, double **B, double **C)	{
	for(int i=0;i<N;i++)	{
		for(int j=0;j<N;j++)	{
			Xt[i][j] = A[j][i];
		}
	}
	mult(B,A,F_temp);
	mult(Xt,F_temp,C);
}

void diag(double **F, double **C, double **E)	{
	double theta=0.5*atan(2*F[0][1]/(F[0][0]-F[1][1]));
	double temp;
	
	C[0][0] = cos(theta);
	C[0][1] = sin(theta);
	C[1][0] = sin(theta);
	C[1][1] = -cos(theta);
	
	E[0][0] = F[0][0]*cos(theta)*cos(theta) + F[1][1]*sin(theta)*sin(theta) + F[0][1]*sin(2*theta);
	E[0][1] = 0;
	E[1][0] = 0;
	E[1][1] = F[1][1]*cos(theta)*cos(theta) + F[0][0]*sin(theta)*sin(theta) - F[0][1]*sin(2*theta);
	
	if(E[0][0]>E[1][1])	{
		temp = E[1][1];
		E[1][1] = E[0][0];
		E[0][0] = temp;
		
		temp = C[0][1];
		C[0][1] = C[0][0];
		C[0][0] = temp;
		
		temp = C[1][1];
		C[1][1] = C[1][0];
		C[1][0] = temp;
	}
}

void make_P()	{		//3.145
	for(int i=0;i<N;i++)	{
		for(int j=0;j<N;j++)	{
			P_old[i][j] = P[i][j];
			P[i][j] = 0;
			for(int k=0;k<N/2;k++)	{
				P[i][j] += 2*C[i][k]*C[j][k];
			}
		}
	}
}

double calculate_energy()	{
	double energy = 0;
	for(int i=0;i<N;i++)	{
		for(int j=0;j<N;j++)	{
			energy +=0.5*P[i][j]*(H_core[i][j]+F[i][j]);
		}
	}
	energy += z1*z2/R;
	return energy;
}