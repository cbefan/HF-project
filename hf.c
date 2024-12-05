#include "hf.h"
// 🤓
void initialization(int argc,char* argv[]) {
	
	// allocate memory
	S = malloc(sizeof(*S));
	T = malloc(sizeof(*T));
	V1 = malloc(sizeof(*V1));
	V2 = malloc(sizeof(*V2));
	H_core = malloc(sizeof(*H_core));
	G = malloc(sizeof(*G));
	X = malloc(sizeof(*X));
	Xt = malloc(sizeof(*Xt));
	P = malloc(sizeof(*P));
	P_old = malloc(sizeof(*P_old));
	P_prime = malloc(sizeof(*P_prime));
	F = malloc(sizeof(*F));
	F_prime = malloc(sizeof(*F_prime));
	F_temp = malloc(sizeof(*F_temp));
	E = malloc(sizeof(*E));
	C = malloc(sizeof(*C));
	C_prime = malloc(sizeof(*C_prime));
	mat2e = malloc(sizeof(***mat2e));
	for(int i=0;i<N;i++) {
		S[i] = malloc(N*sizeof(double));
		T[i] = malloc(N*sizeof(double));
		V1[i] = malloc(N*sizeof(double));
		V2[i] = malloc(N*sizeof(double));
		H_core[i] = malloc(N*sizeof(double));
		G[i] = malloc(N*sizeof(double));
		X[i] = malloc(N*sizeof(double));
		Xt[i] = malloc(N*sizeof(double));
		P[i] = malloc(N*sizeof(double));
		P_old[i] = malloc(N*sizeof(double));
		P_prime[i] = malloc(N*sizeof(double));
		F[i] = malloc(N*sizeof(double));
		F_prime[i] = malloc(N*sizeof(double));
		F_temp[i] = malloc(N*sizeof(double));
		E[i] = malloc(N*sizeof(double));
		C[i] = malloc(N*sizeof(double));
		C_prime[i] = malloc(N*sizeof(double));
		mat2e[i] = malloc(sizeof(**mat2e));
		for(int j=0;j<N;j++)	{
			mat2e[i][j] = malloc(sizeof(*mat2e));
			for(int k=0;k<N;k++)	{
				mat2e[i][j][k] = malloc(N*sizeof(double));
			}
		}
	}
	
	if(argc>1)	{
		printf("Reading arguments\n");
		if(argc != 4)	{
			printf("Expecting three arguments: identity of atom 1, identity of atom 2, and separation distance\n");
			exit(1);
		}
		if(strcmp(argv[1],"H")==0)	{
			printf("Atom 1\t= Hydrogen\n");
			z1 = 1;
			zeta1 = 1.24;
		}
		else if(strcmp(argv[1],"He")==0)	{
			printf("Atom 1\t= Helium\n");
			z1 = 2;
			zeta1 = 2.0925;
		}
		else	{
			printf("Atom 1 not recognized, exiting\n");
			exit(1);
		}
		
		if(strcmp(argv[2],"H")==0)	{
			printf("Atom 2\t= Hydrogen\n");
			z2 = 1;
			zeta2 = 1.24;
		}
		else if(strcmp(argv[2],"He")==0)	{
			printf("Atom 2\t= Helium\n");
			z2 = 2;
			zeta2 = 2.0925;
		}
		else	{
			printf("Atom 2 not recognized, exiting\n");
			exit(1);
		}
		
		sscanf(argv[3],"%lf",&R);
		printf("R\t= %lf\n",R);
	}
	else	{
		printf("Using default settings\n");
		printf("Atom 1	= Hydrogen\n");
		printf("Atom 2	= Hydrogen\n");
		printf("R	= %lf\n",R);
	}
	printf("\n");
	
	r[0] = 0;
	r[1] = R;
	
	// adjust coefs and alphas from zeta_ref=1 to actual zeta
	for(int i=0;i<3;i++) {
		alphas1[i] *= zeta1*zeta1;
		coefs1[i] *= pow(2*alphas1[i]/PI,0.75);
		alphas2[i] *= zeta2*zeta2;
		coefs2[i] *= pow(2*alphas2[i]/PI,0.75);
	}
	
	// calculate little h pieces
	for (int i = 0; i<N; i++) {		// go through pairs of atoms
		for(int j = 0;j<N; j++) {
			rab = r[i]-r[j];
			rab2 =rab*rab;
			for(int k=0;k<3;k++) {	// go through pairs of basis functions
				if(i==0)	{		// assign parameters based on atoms
					a = alphas1[k];
					d1 = coefs1[k];
				}
				else	{
					a = alphas2[k];
					d1 = coefs2[k];
				}
				
				for(int l=0;l<3;l++) {
					if(j==0)	{
						b = alphas1[l];
						d2 = coefs1[l];
					}
					else	{
						b = alphas2[l];
						d2 = coefs2[l];
					}
					rp = (a*r[i]+b*r[j])/(a+b);
					rpa = rp - r[0];
					rpb = rp - r[1];
					
					S[i][j] += d1*d2 * pow((PI/(a+b)),1.5)*exp(-a*b/(a+b)*rab2);
					T[i][j] += d1*d2 * a*b/(a+b)*(3-2*a*b/(a+b)*rab2)*pow((PI/(a+b)),1.5)*exp(-a*b/(a+b)*rab2);
					t = (a+b)*rpa*rpa;
					V1[i][j] += d1*d2 * -2*PI/(a+b)*z1*exp(-a*b/(a+b)*rab2)*F0(t);
					t = (a+b)*rpb*rpb;
					V2[i][j] += d1*d2 * -2*PI/(a+b)*z2*exp(-a*b/(a+b)*rab2)*F0(t);
				}
			}
			H_core[i][j] = T[i][j] + V1[i][j] + V2[i][j];
		}
	}
	
	// calculate two electron integrals
	for(int i=0;i<3;i++)	{	// go through all permutations of the 3 basis functions
		for(int j=0;j<3;j++)	{
			for(int k=0;k<3;k++)	{
				for(int l=0;l<3;l++)	{
					R2 = R*R;

					rq = (alphas1[k]*r[0]+alphas2[l]*r[1])/(alphas1[k]+alphas2[l]);
					rp = (alphas1[i]*r[0]+alphas2[j]*r[1])/(alphas1[i]+alphas2[j]);
					raq = r[0]-rq;
					rbq = r[1] - rq;
					rpq = rp - rq;
					rpb = rp - r[1];
					
					raq2 = raq*raq;
					rbq2 = rbq*rbq;
					rpq2 = rpq*rpq;
					rpb2 = rpb*rpb;
					
					//printf("%lf %lf %lf\n",raq2,rbq2,rpq2);
					
					V1111 += coefs1[i]*coefs1[j]*coefs1[k]*coefs1[l] * integral(alphas1[i],alphas1[j],alphas1[k],alphas1[l],0,0,0);
					V1112 += coefs1[i]*coefs1[j]*coefs1[k]*coefs2[l] * integral(alphas1[i],alphas1[j],alphas1[k],alphas2[l],0,R2,raq2);
					V1122 += coefs1[i]*coefs1[j]*coefs2[k]*coefs2[l] * integral(alphas1[i],alphas1[j],alphas2[k],alphas2[l],0,0,R2);
					V1222 += coefs1[i]*coefs2[j]*coefs2[k]*coefs2[l] * integral(alphas1[i],alphas2[j],alphas2[k],alphas2[l],R2,0,rpb2);
					V1212 += coefs1[i]*coefs2[j]*coefs1[k]*coefs2[l] * integral(alphas1[i],alphas2[j],alphas1[k],alphas2[l],R2,R2,rpq2);
					V2222 += coefs2[i]*coefs2[j]*coefs2[k]*coefs2[l] * integral(alphas2[i],alphas2[j],alphas2[k],alphas2[l],0,0,0);
				}
			}
		}
	}
	
	// calculate canonical orthogonalization matrix
	double s1 = 1+S[0][1];
	double s2 = 1-S[0][1];
	X[0][0] = 1/sqrt(2*s1);
	X[0][1] = 1/sqrt(2*s2);
	X[1][0] = 1/sqrt(2*s1);
	X[1][1] = -1/sqrt(2*s2);
	
	// set up two electron integral matrix
	mat2e[0][0][0][0] = V1111;
	mat2e[0][0][0][1] = V1112;
	mat2e[0][0][1][0] = V1112;
	mat2e[0][1][0][0] = V1112;
	mat2e[1][0][0][0] = V1112;
	mat2e[0][0][1][1] = V1122;
	mat2e[0][1][0][1] = V1212;
	mat2e[1][0][0][1] = V1212;
	mat2e[0][1][1][0] = V1212;
	mat2e[1][0][1][0] = V1212;
	mat2e[1][1][0][0] = V1122;
	mat2e[0][1][1][1] = V1222;
	mat2e[1][0][1][1] = V1222;
	mat2e[1][1][0][1] = V1222;
	mat2e[1][1][1][0] = V1222;
	mat2e[1][1][1][1] = V2222;
	
	// initialize population matrix
	P[0][0] = 0;
	P[0][1] = 0;
	P[1][0] = 0;
	P[1][1] = 0;
}

void scf()	{
	int iter = 0;
	double delta;
	while(iter<max_iterations)	{
		++iter;
		make_G();
		make_F();
		transform_XF(X,F,F_prime);
		diag(F_prime,C_prime,E);
		mult(X,C_prime,C);
		make_P();
		
		delta = 0;
		for(int i=0;i<N;i++)	{
			for(int j=0;j<N;j++)	{
				delta += (P[i][j]-P_old[i][j])*(P[i][j]-P_old[i][j]);
			}
		}
		delta = sqrt(delta/4);
		printf("delta	%lf\n",delta);
		if(delta<convergence_threshold)	{
			printf("convergence reached after %d iterations, exiting scf loop\n\n",iter);
			break;
		}
	}
	if(iter==max_iterations)	{
		printf("convergence NOT reached after %d iterations, exiting scf loop\n\n",iter);
	}
	E_total = calculate_energy();
	
	// Mulliken population analysis
	mult(P,S,P_prime);
	N_A = P_prime[0][0];
	N_B =  P_prime[1][1];
}

int main(int argc, char* argv[]) {
	printf("Read %d command line arguments\n\n",argc-1);
	
	initialization(argc,argv);
	printf("S	%lf %lf\n	%lf %lf\n",S[0][0],S[0][1],S[1][0],S[1][1]);	//3.239
	printf("T	%lf %lf\n	%lf %lf\n",T[0][0],T[0][1],T[1][0],T[1][1]);	//3.230
	printf("V1	%lf %lf\n	%lf %lf\n",V1[0][0],V1[0][1],V1[1][0],V1[1][1]);			//3.231
	printf("V2	%lf %lf\n	%lf %lf\n",V2[0][0],V2[0][1],V2[1][0],V2[1][1]);		//3.232
	printf("H_core	%lf %lf\n	%lf %lf\n",H_core[0][0],H_core[0][1],H_core[1][0],H_core[1][1]);	//3.233
	printf("X	%lf %lf\n	%lf %lf\n",X[0][0],X[0][1],X[1][0],X[1][1]);	//3.262
	printf("V1111	%lf\nV1112	%lf\nV1122	%lf\nV1222	%lf\nV1212	%lf\nV2222	%lf\n\n",V1111,V1112,V1122,V1222,V1212,V2222);	//3.235
	
	max_iterations = 100;
	convergence_threshold = 1e-6;
	
	scf();
	printf("G	%lf %lf\n	%lf %lf\n",G[0][0],G[0][1],G[1][0],G[1][1]);
	printf("F	%lf %lf\n	%lf %lf\n",F[0][0],F[0][1],F[1][0],F[1][1]);
	printf("F'	%lf %lf\n	%lf %lf\n",F_prime[0][0],F_prime[0][1],F_prime[1][0],F_prime[1][1]);
	printf("E	%lf %lf\n	%lf %lf\n",E[0][0],E[0][1],E[1][0],E[1][1]);
	printf("C	%lf %lf\n	%lf %lf\n",C[0][0],C[0][1],C[1][0],C[1][1]);	//3.238
	printf("P	%lf %lf\n	%lf %lf\n",P[0][0],P[0][1],P[1][0],P[1][1]);		//3.239
	printf("P_old	%lf %lf\n	%lf %lf\n",P_old[0][0],P_old[0][1],P_old[1][0],P_old[1][1]);
	printf("P_prime	%lf %lf\n	%lf %lf\n",P_prime[0][0],P_prime[0][1],P_prime[1][0],P_prime[1][1]);
	printf("N_1	%lf\n",N_A);
	printf("N_2	%lf\n",N_B);
	printf("E_tot	%lf\n",E_total);		//ex 3.27
	return 0;
}
