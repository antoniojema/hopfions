#include <iostream>
#include <cmath>
#include <fstream>
//#include <H5Cpp.h>

using namespace std;
//using namespace H5;

double A, Xe, Ye, Ze;
double B, Xm, Ym, Zm;
const int N=100, iterations=100;
double Ex[(N+1)*(N+1)*(N+1)],  Ey[(N+1)*(N+1)*(N+1)],  Ez[(N+1)*(N+1)*(N+1)];
double Hx[(N+1)*(N+1)*(N+1)],  Hy[(N+1)*(N+1)*(N+1)],  Hz[(N+1)*(N+1)*(N+1)];
/*double Ex1[(N+1)*(N+1)*(N+1)], Ey1[(N+1)*(N+1)*(N+1)], Ez1[(N+1)*(N+1)*(N+1)];
double Hx1[(N+1)*(N+1)*(N+1)], Hy1[(N+1)*(N+1)*(N+1)], Hz1[(N+1)*(N+1)*(N+1)];
double Ex2[(N+1)*(N+1)*(N+1)], Ey2[(N+1)*(N+1)*(N+1)], Ez2[(N+1)*(N+1)*(N+1)];
double Hx2[(N+1)*(N+1)*(N+1)], Hy2[(N+1)*(N+1)*(N+1)], Hz2[(N+1)*(N+1)*(N+1)];
*/

int ind(int i, int j, int k){
	return (N+1)*(N+1)*i+(N+1)*j+k;
}

int main(){
	int n, i, j, k;
	double sigma, sigmam, Dx, Dy, Dz, Dt, epsilon, mu, c, L, x;
	ofstream fout, fout2;
	
	L=10;
	sigma = sigmam = 0;
	epsilon = mu = 1;
	c=1;
	Dx = Dy = Dz = L/N;
	Dt = 0.8*Dx/(c*sqrt(3.0));
	
	A = (1.-sigma*Dt/(2.*epsilon))/(1.+sigma*Dt/(2.*epsilon));
	Xe = (Dt/(epsilon*Dx))/(1.+sigma*Dt/(2.*epsilon));
	Ye = (Dt/(epsilon*Dy))/(1.+sigma*Dt/(2.*epsilon));
	Ze = (Dt/(epsilon*Dz))/(1.+sigma*Dt/(2.*epsilon));
	B = (1.-sigmam*Dt/(2.*mu))/(1.+sigmam*Dt/(2.*mu));
	Xm = (Dt/(mu*Dx))/(1.+sigmam*Dt/(2.*mu));
	Ym = (Dt/(mu*Dy))/(1.+sigmam*Dt/(2.*mu));
	Zm = (Dt/(mu*Dz))/(1.+sigmam*Dt/(2.*mu));
	
	//Open output file
	fout.open("3D.txt");
	fout2.open("3Dcm.txt");
	
	//Initial conditions
	for(i=0; i<=N; i++){
		for(j=0; j<=N; j++){
			for(k=0; k<=N; k++){
				Ex[ind(i,j,k)] = 0;
				Ey[ind(i,j,k)] = 0;
				Ez[ind(i,j,k)] = 0;
				Hx[ind(i,j,k)] = 0;
				Hy[ind(i,j,k)] = 0;
				Hz[ind(i,j,k)] = 0;
			}
		}
	}
	/*
	for(i=0; i<=N; i++){
		for(j=0; j<=N; j++){
			for(k=0; k<=N; k++){
				Ex1[ind(i,j,k)] = Ex2[ind(i,j,k)] = Ey1[ind(i,j,k)] = Ey2[ind(i,j,k)] = Ez1[ind(i,j,k)] = Ez2[ind(i,j,k)]=0;
				Hx1[ind(i,j,k)] = Hx2[ind(i,j,k)] = Hy1[ind(i,j,k)] = Hy2[ind(i,j,k)] = Hz1[ind(i,j,k)] = Hz2[ind(i,j,k)]=0;
			}
		}
	}
	*/
	
	for(n=0; n<=iterations; n++){
		//Calculate E [LACKS ALL BOUNDING CONDITIONS]
		//Ex
		for(i=0; i<=N-1; i++){
			for(j=0; j<=N-1; j++){
				for(k=0; k<=N-1; k++){
					Ex[ind(i,j,k)] = A*Ex[ind(i,j,k)] + Ze*(Hy[ind(i,j,k)]-Hy[ind(i,j,k+1)]) + Ye*(Hz[ind(i,j+1,k)]-Hz[ind(i,j,k)]);
					Ey[ind(i,j,k)] = A*Ey[ind(i,j,k)] + Xe*(Hz[ind(i,j,k)]-Hz[ind(i+1,j,k)]) + Ze*(Hx[ind(i,j,k+1)]-Hx[ind(i,j,k)]);
					Ez[ind(i,j,k)] = A*Ez[ind(i,j,k)] + Ye*(Hx[ind(i,j,k)]-Hx[ind(i,j+1,k)]) + Xe*(Hy[ind(i+1,j,k)]-Hy[ind(i,j,k)]);
					//Corriente
					if (i>=4*N/10 && i<6*N/10 && j>=4*N/10 && j<6*N/10 && k>=4*N/10 && k<6*N/10 && n<=31){
						Ex[ind(i,j,k)] -= exp(-0.5*((i-N/2)*(i-N/2)+(j-N/2)*(j-N/2)+(k-N/2)*(k-N/2)))*sin(0.2*n);
					}
				}
				Ez[ind(i,j,N)] = A*Ez[ind(i,j,N)] + Ye*(Hx[ind(i,j,N)]-Hx[ind(i,j+1,N)]) + Xe*(Hy[ind(i+1,j,N)]-Hy[ind(i,j,N)]);
			}
			for(k=0; k<=N-1; k++){
				Ey[ind(i,N,k)] = A*Ey[ind(i,N,k)] + Xe*(Hz[ind(i,N,k)]-Hz[ind(i+1,N,k)]) + Ze*(Hx[ind(i,N,k+1)]-Hx[ind(i,N,k)]);
			}
		}
		for(j=0; j<=N-1; j++){
			for(k=0; k<=N-1;k++){
				Ex[ind(N,j,k)] = A*Ex[ind(N,j,k)] + Ze*(Hy[ind(N,j,k)]-Hy[ind(N,j,k+1)]) + Ye*(Hz[ind(N,j+1,k)]-Hz[ind(N,j,k)]);
			}
		}
		
		/*
		//E -> E1 -> E2
		for(i=0; i<=N; i++){
			for(j=0; j<=N; j++){
				for(k=0; k<=N; k++){
					Ex2[ind(i,j,k)] = Ex1[ind(i,j,k)];
					Ex1[ind(i,j,k)] = Ex[ind(i,j,k)];
					Ey2[ind(i,j,k)] = Ey1[ind(i,j,k)];
					Ey1[ind(i,j,k)] = Ey[ind(i,j,k)];
					Ez2[ind(i,j,k)] = Ez1[ind(i,j,k)];
					Ez1[ind(i,j,k)] = Ez[ind(i,j,k)];
				}
			}
		}
		*/
		
		//Calculate H [LACKS ALL BOUNDING CONDITIONS]
		for(i=1; i<=N; i++){
			for(j=1; j<=N; j++){
				for(k=1; k<=N; k++){
					Hx[ind(i,j,k)] = B*Hx[ind(i,j,k)] + Zm*(Ey[ind(i,j,k)]-Ey[ind(i,j,k-1)]) + Ym*(Ez[ind(i,j-1,k)]-Ez[ind(i,j,k)]);
					Hy[ind(i,j,k)] = B*Hy[ind(i,j,k)] + Xm*(Ez[ind(i,j,k)]-Ez[ind(i-1,j,k)]) + Zm*(Ex[ind(i,j,k-1)]-Ex[ind(i,j,k)]);
					Hz[ind(i,j,k)] = B*Hz[ind(i,j,k)] + Ym*(Ex[ind(i,j,k)]-Ex[ind(i,j-1,k)]) + Xm*(Ey[ind(i-1,j,k)]-Ey[ind(i,j,k)]);
				}
				Hz[ind(i,j,0)] = B*Hz[ind(i,j,0)] + Ym*(Ex[ind(i,j,0)]-Ex[ind(i,j-1,0)]) + Xm*(Ey[ind(i-1,j,0)]-Ey[ind(i,j,0)]);
			}
			for(k=1; k<=N; k++){
				Hy[ind(i,0,k)] = B*Hy[ind(i,0,k)] + Xm*(Ez[ind(i,0,k)]-Ez[ind(i-1,0,k)]) + Zm*(Ex[ind(i,0,k-1)]-Ex[ind(i,0,k)]);
			}
		}
		for(j=1; j<=N; j++){
			for(k=1; k<=N;k++){
				Hx[ind(0,j,k)] = B*Hx[ind(0,j,k)] + Zm*(Ey[ind(0,j,k)]-Ey[ind(0,j,k-1)]) + Ym*(Ez[ind(0,j-1,k)]-Ez[ind(0,j,k)]);
			}
		}
		
		/*
		//H -> H1 -> H2
		for(i=0; i<=N; i++){
			for(j=0; j<=N; j++){
				for(k=0; k<=N; k++){
					Hx2[ind(i,j,k)] = Hx1[ind(i,j,k)];
					Hx1[ind(i,j,k)] = Hx[ind(i,j,k)];
					Hy2[ind(i,j,k)] = Hy1[ind(i,j,k)];
					Hy1[ind(i,j,k)] = Hy[ind(i,j,k)];
					Hz2[ind(i,j,k)] = Hz1[ind(i,j,k)];
					Hz1[ind(i,j,k)] = Hz[ind(i,j,k)];
				}
			}
		}
		*/
		
		//Data -> txt
		for(i=0; i<=N; i++){
			fout << i << "	" << abs(Ex[ind(N/2,N/2,i)]) << endl;
			for(j=0; j<=N; j++){
				//This is to make a colormap of a cross section of the entire scene in order to see the radiation scheme
				fout2 << abs(Ex[ind(i,N/2,j)]) << "	";
			}
			fout2 << endl;
		}
	
	}
	fout.close();
	fout2.close();
	
	return 0;
}
