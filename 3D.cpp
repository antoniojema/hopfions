#include <iostream>
#include <cmath>
#include <fstream>
#include <H5Cpp.h>

using namespace std;
using namespace H5;

double A, Xe, Ye, Ze;
double B, Xm, Ym, Zm;
double C1, C2, C3;
const int N=100, iterations=200;
double Ex[(N+1)*(N+1)*(N+1)],  Ey[(N+1)*(N+1)*(N+1)],  Ez[(N+1)*(N+1)*(N+1)];
double Hx[(N+1)*(N+1)*(N+1)],  Hy[(N+1)*(N+1)*(N+1)],  Hz[(N+1)*(N+1)*(N+1)];
double Ex1[(N+1)*(N+1)*(N+1)], Ey1[(N+1)*(N+1)*(N+1)], Ez1[(N+1)*(N+1)*(N+1)];
double Hx1[(N+1)*(N+1)*(N+1)], Hy1[(N+1)*(N+1)*(N+1)], Hz1[(N+1)*(N+1)*(N+1)];
double Ex2[(N+1)*(N+1)*(N+1)], Ey2[(N+1)*(N+1)*(N+1)], Ez2[(N+1)*(N+1)*(N+1)];
double Hx2[(N+1)*(N+1)*(N+1)], Hy2[(N+1)*(N+1)*(N+1)], Hz2[(N+1)*(N+1)*(N+1)];

int ind(int i, int j, int k){
	return (N+1)*(N+1)*i+(N+1)*j+k;
}

int main(){
	int n, i, j, k;
	double sigma, sigmam, Dx, Dy, Dz, Dt, epsilon, mu, c, L, x;
	
	L=10;
	sigma = sigmam = 0;
	epsilon = mu = 1;
	c=1;
	Dx = Dy = Dz = L/N;
	Dt = 0.8*Dx/(c*sqrt(3.0));
	
	//Constants for the FDTD method
	A = (1.-sigma*Dt/(2.*epsilon))/(1.+sigma*Dt/(2.*epsilon));
	Xe = (Dt/(epsilon*Dx))/(1.+sigma*Dt/(2.*epsilon));
	Ye = (Dt/(epsilon*Dy))/(1.+sigma*Dt/(2.*epsilon));
	Ze = (Dt/(epsilon*Dz))/(1.+sigma*Dt/(2.*epsilon));
	B = (1.-sigmam*Dt/(2.*mu))/(1.+sigmam*Dt/(2.*mu));
	Xm = (Dt/(mu*Dx))/(1.+sigmam*Dt/(2.*mu));
	Ym = (Dt/(mu*Dy))/(1.+sigmam*Dt/(2.*mu));
	Zm = (Dt/(mu*Dz))/(1.+sigmam*Dt/(2.*mu));
	
	//Constants for the MUR boundary conditions
	C1 = (c*Dt-Dx)/(c*Dt+Dx);
	C2 = 2*Dx/(c*Dt+Dx);
	C3 = c*c*Dt*Dt/(2*Dx*(c*Dt+Dx));
	
	/* Initial conditions */
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
	
	for(i=0; i<=N; i++){
		for(j=0; j<=N; j++){
			for(k=0; k<=N; k++){
				Ex1[ind(i,j,k)] = 0;
				Ex2[ind(i,j,k)] = 0;
				Ey1[ind(i,j,k)] = 0;
				Ey2[ind(i,j,k)] = 0;
				Ez1[ind(i,j,k)] = 0;
				Ez2[ind(i,j,k)] = 0;
				Hx1[ind(i,j,k)] = 0;
				Hx2[ind(i,j,k)] = 0;
				Hy1[ind(i,j,k)] = 0;
				Hy2[ind(i,j,k)] = 0;
				Hz1[ind(i,j,k)] = 0;
				Hz2[ind(i,j,k)] = 0;
			}
		}
	}
	
	try
	{
	Exception::dontPrint();
	
	/* Prepare output file */
	H5File fout("3D.h5", H5F_ACC_TRUNC);
	hsize_t dimsf[3];
	dimsf[0] = N+1;
	dimsf[1] = N+1;
	dimsf[2] = N+1;
	DataSpace dspace(3, dimsf);
	Group grp;
	DataSet dset;
	
	for(n=0; n<=iterations; n++){
		/***** E *****/ //LACKS EDGES CONDITIONS
		for(i=0; i<=N-1; i++){
			for(j=0; j<=N-1; j++){
				for(k=0; k<=N-1; k++){
					Ex[ind(i,j,k)] = A * Ex[ind(i,j,k)] + Ze * (Hy[ind(i,j,k)] - Hy[ind(i,j,k+1)]) +
									 Ye * (Hz[ind(i,j+1,k)] - Hz[ind(i,j,k)]);
					Ey[ind(i,j,k)] = A * Ey[ind(i,j,k)] + Xe * (Hz[ind(i,j,k)] - Hz[ind(i+1,j,k)]) +
									 Ze * (Hx[ind(i,j,k+1)] - Hx[ind(i,j,k)]);
					Ez[ind(i,j,k)] = A * Ez[ind(i,j,k)] + Ye * (Hx[ind(i,j,k)] - Hx[ind(i,j+1,k)]) +
									 Xe * (Hy[ind(i+1,j,k)] - Hy[ind(i,j,k)]);
					/** Current **/
					if (i>=4*N/10 && i<6*N/10 && j>=4*N/10 && j<6*N/10 && k>=4*N/10 && k<6*N/10 /*&& n<=31*/){
						Ex[ind(i,j,k)] -= exp(-0.5*((i-N/2)*(i-N/2)+(j-N/2)*(j-N/2)+(k-N/2)*(k-N/2)))*sin(0.2*n);
					}
					/***************/
				}
				Ez[ind(i,j,N)] = A * Ez[ind(i,j,N)] + Ye * (Hx[ind(i,j,N)] - Hx[ind(i,j+1,N)]) +
								 Xe * (Hy[ind(i+1,j,N)] - Hy[ind(i,j,N)]);
				Ey[ind(i,N,j)] = A * Ey[ind(i,N,j)] + Xe * (Hz[ind(i,N,j)] - Hz[ind(i+1,N,j)]) +
								 Ze * (Hx[ind(i,N,j+1)] - Hx[ind(i,N,j)]);
				Ex[ind(N,i,j)] = A * Ex[ind(N,i,j)] + Ze * (Hy[ind(N,i,j)] - Hy[ind(N,i,j+1)]) +
								 Ye * (Hz[ind(N,i+1,j)] - Hz[ind(N,i,j)]);
			}
		}
		
		/* MUR */
		for(i=1; i<=N-1; i++){
			for(j=1; j<=N-1; j++){
				/* Face XY */
				Ex[ind(i,j,N)] = -1.*Ex2[ind(i,j,N-1)] + C1 * (Ex[ind(i,j,N-1)] + Ex2[ind(i,j,N)]) +
								 C2 * (Ex1[ind(i,j,N)] + Ex1[ind(i,j,N-1)]) +
								 C3 * (Ex1[ind(i+1,j,N)] + Ex1[ind(i-1,j,N)] + Ex1[ind(i,j+1,N)] + Ex1[ind(i,j-1,N)] +
									   Ex1[ind(i+1,j,N-1)] + Ex1[ind(i-1,j,N-1)] + Ex1[ind(i,j+1,N-1)] + Ex1[ind(i,j-1,N-1)] -
									   4 * Ex1[ind(i,j,N)] - 4 * Ex1[ind(i,j,N-1)]);
				Ey[ind(i,j,N)] = -1.*Ey2[ind(i,j,N-1)] + C1 * (Ey[ind(i,j,N-1)] + Ey2[ind(i,j,N)]) +
								 C2 * (Ey1[ind(i,j,N)] + Ey1[ind(i,j,N-1)]) +
								 C3 * (Ey1[ind(i+1,j,N)] + Ey1[ind(i-1,j,N)] + Ey1[ind(i,j+1,N)] + Ey1[ind(i,j-1,N)] +
									   Ey1[ind(i+1,j,N-1)] + Ey1[ind(i-1,j,N-1)] + Ey1[ind(i,j+1,N-1)] + Ey1[ind(i,j-1,N-1)] -
									   4 * Ey1[ind(i,j,N)] - 4 * Ey1[ind(i,j,N-1)]);
				/* Face XZ */
				Ex[ind(i,N,j)] = -1.*Ex2[ind(i,N-1,j)] + C1 * (Ex[ind(i,N-1,j)] + Ex2[ind(i,N,j)]) +
								 C2 * (Ex1[ind(i,N,j)] + Ex1[ind(i,N-1,j)]) +
								 C3 * (Ex1[ind(i+1,N,j)] + Ex1[ind(i-1,N,j)] + Ex1[ind(i,N,j+1)] + Ex1[ind(i,N,j-1)] +
									   Ex1[ind(i+1,N-1,j)] + Ex1[ind(i-1,N-1,j)] + Ex1[ind(i,N-1,j+1)] + Ex1[ind(i,N-1,j-1)] -
									   4 * Ex1[ind(i,N,j)] - 4 * Ex1[ind(i,N-1,j)]);
				Ez[ind(i,N,j)] = -1.*Ez2[ind(i,N-1,j)] + C1 * (Ez[ind(i,N-1,j)] + Ez2[ind(i,N,j)]) +
								 C2 * (Ez1[ind(i,N,j)] + Ez1[ind(i,N-1,j)]) +
								 C3 * (Ez1[ind(i+1,N,j)] + Ez1[ind(i-1,N,j)] + Ez1[ind(i,N,j+1)] + Ez1[ind(i,N,j-1)] +
									   Ez1[ind(i+1,N-1,j)] + Ez1[ind(i-1,N-1,j)] + Ez1[ind(i,N-1,j+1)] + Ez1[ind(i,N-1,j-1)] -
									   4 * Ez1[ind(i,N,j)] - 4 * Ez1[ind(i,N-1,j)]);
				/* Face YZ */
				Ey[ind(N,i,j)] = -1.*Ey2[ind(N-1,i,j)] + C1 * (Ey[ind(N-1,i,j)] + Ey2[ind(N,i,j)]) +
								 C2 * (Ey1[ind(N,i,j)] + Ey1[ind(N-1,i,j)]) +
								 C3 * (Ey1[ind(N,i+1,j)] + Ey1[ind(N,i-1,j)] + Ey1[ind(N,i,j+1)] + Ey1[ind(N,i,j-1)] +
									   Ey1[ind(N-1,i+1,j)] + Ey1[ind(N-1,i-1,j)] + Ey1[ind(N-1,i,j+1)] + Ey1[ind(N-1,i,j-1)] -
									   4 * Ey1[ind(N,i,j)] - 4 * Ey1[ind(N-1,i,j)]);
				Ez[ind(N,i,j)] = -1.*Ez2[ind(N-1,i,j)] + C1 * (Ez[ind(N-1,i,j)] + Ez2[ind(N,i,j)]) +
								 C2 * (Ez1[ind(N,i,j)] + Ez1[ind(N-1,i,j)]) +
								 C3 * (Ez1[ind(N,i+1,j)] + Ez1[ind(N,i-1,j)] + Ez1[ind(N,i,j+1)] + Ez1[ind(N,i,j-1)] +
									   Ez1[ind(N-1,i+1,j)] + Ez1[ind(N-1,i-1,j)] + Ez1[ind(N-1,i,j+1)] + Ez1[ind(N-1,i,j-1)] -
									   4 * Ez1[ind(N,i,j)] - 4 * Ez1[ind(N-1,i,j)]);
			}
		}
		
		/* E -> E1 -> E2 */
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
		
		/***** H *****/ //LACKS EDGES CONDITIONS
		for(i=1; i<=N; i++){
			for(j=1; j<=N; j++){
				for(k=1; k<=N; k++){
					Hx[ind(i,j,k)] = B * Hx[ind(i,j,k)] + Zm * (Ey[ind(i,j,k)] - Ey[ind(i,j,k-1)]) +
									 Ym * (Ez[ind(i,j-1,k)] - Ez[ind(i,j,k)]);
					Hy[ind(i,j,k)] = B * Hy[ind(i,j,k)] + Xm * (Ez[ind(i,j,k)] - Ez[ind(i-1,j,k)]) +
									 Zm * (Ex[ind(i,j,k-1)] - Ex[ind(i,j,k)]);
					Hz[ind(i,j,k)] = B * Hz[ind(i,j,k)] + Ym * (Ex[ind(i,j,k)] - Ex[ind(i,j-1,k)]) +
									 Xm * (Ey[ind(i-1,j,k)] - Ey[ind(i,j,k)]);
				}
				Hz[ind(i,j,0)] = B * Hz[ind(i,j,0)] + Ym * (Ex[ind(i,j,0)] - Ex[ind(i,j-1,0)]) +
								 Xm * (Ey[ind(i-1,j,0)] - Ey[ind(i,j,0)]);
				Hy[ind(i,0,j)] = B * Hy[ind(i,0,j)] + Xm * (Ez[ind(i,0,j)] - Ez[ind(i-1,0,j)]) +
								 Zm * (Ex[ind(i,0,j-1)] - Ex[ind(i,0,j)]);
				Hx[ind(0,i,j)] = B * Hx[ind(0,i,j)] + Zm * (Ey[ind(0,i,j)] - Ey[ind(0,i,j-1)]) +
								 Ym * (Ez[ind(0,i-1,j)] - Ez[ind(0,i,j)]);
			}
		}
		
		/* MUR */
		for(i=1; i<=N-1; i++){
			for(j=1; j<=N-1; j++){
				/* Face XY */
				Hx[ind(i,j,0)] = -1.*Hx2[ind(i,j,1)] + C1 * (Hx[ind(i,j,1)] + Hx2[ind(i,j,0)]) +
								 C2 * (Hx1[ind(i,j,0)] + Hx1[ind(i,j,1)]) +
								 C3 * (Hx1[ind(i+1,j,0)] + Hx1[ind(i-1,j,0)] + Hx1[ind(i,j+1,0)] + Hx1[ind(i,j-1,0)] +
									   Hx1[ind(i+1,j,1)] + Hx1[ind(i-1,j,1)] + Hx1[ind(i,j+1,1)] + Hx1[ind(i,j-1,1)] -
									   4 * Hx1[ind(i,j,0)] - 4 * Hx1[ind(i,j,1)]);
				Hy[ind(i,j,0)] = -1.*Hy2[ind(i,j,1)] + C1 * (Hy[ind(i,j,1)] + Hy2[ind(i,j,0)]) +
								 C2 * (Hy1[ind(i,j,0)] + Hy1[ind(i,j,1)]) +
								 C3 * (Hy1[ind(i+1,j,0)] + Hy1[ind(i-1,j,0)] + Hy1[ind(i,j+1,0)] + Hy1[ind(i,j-1,0)] +
									   Hy1[ind(i+1,j,1)] + Hy1[ind(i-1,j,1)] + Hy1[ind(i,j+1,1)] + Hy1[ind(i,j-1,1)] -
									   4 * Hy1[ind(i,j,0)] - 4 * Hy1[ind(i,j,1)]);
				/* Face XZ */
				Hx[ind(i,0,j)] = -1.*Hx2[ind(i,1,j)] + C1 * (Hx[ind(i,1,j)] + Hx2[ind(i,0,j)]) +
								 C2 * (Hx1[ind(i,0,j)] + Hx1[ind(i,1,j)]) +
								 C3 * (Hx1[ind(i+1,0,j)] + Hx1[ind(i-1,0,j)] + Hx1[ind(i,0,j+1)] + Hx1[ind(i,0,j-1)] +
									   Hx1[ind(i+1,1,j)] + Hx1[ind(i-1,1,j)] + Hx1[ind(i,1,j+1)] + Hx1[ind(i,1,j-1)] -
									   4 * Hx1[ind(i,0,j)] - 4 * Hx1[ind(i,1,j)]);
				Hz[ind(i,0,j)] = -1.*Hz2[ind(i,1,j)] + C1 * (Hz[ind(i,1,j)] + Hz2[ind(i,0,j)]) +
								 C2 * (Hz1[ind(i,0,j)] + Hz1[ind(i,1,j)]) +
								 C3 * (Hz1[ind(i+1,0,j)] + Hz1[ind(i-1,0,j)] + Hz1[ind(i,0,j+1)] + Hz1[ind(i,0,j-1)] +
									   Hz1[ind(i+1,1,j)] + Hz1[ind(i-1,1,j)] + Hz1[ind(i,1,j+1)] + Hz1[ind(i,1,j-1)] -
									   4 * Hz1[ind(i,0,j)] - 4 * Hz1[ind(i,1,j)]);
				/* Face YZ */
				Hy[ind(0,i,j)] = -1.*Hy2[ind(1,i,j)] + C1 * (Hy[ind(1,i,j)] + Hy2[ind(0,i,j)]) +
								 C2 * (Hy1[ind(0,i,j)] + Hy1[ind(1,i,j)]) +
								 C3 * (Hy1[ind(0,i+1,j)] + Hy1[ind(0,i-1,j)] + Hy1[ind(0,i,j+1)] + Hy1[ind(0,i,j-1)] +
									   Hy1[ind(1,i+1,j)] + Hy1[ind(1,i-1,j)] + Hy1[ind(1,i,j+1)] + Hy1[ind(1,i,j-1)] -
									   4 * Hy1[ind(0,i,j)] - 4 * Hy1[ind(1,i,j)]);
				Hz[ind(0,i,j)] = -1.*Hz2[ind(1,i,j)] + C1 * (Hz[ind(1,i,j)] + Hz2[ind(0,i,j)]) +
								 C2 * (Hz1[ind(0,i,j)] + Hz1[ind(1,i,j)]) +
								 C3 * (Hz1[ind(0,i+1,j)] + Hz1[ind(0,i-1,j)] + Hz1[ind(0,i,j+1)] + Hz1[ind(0,i,j-1)] +
									   Hz1[ind(1,i+1,j)] + Hz1[ind(1,i-1,j)] + Hz1[ind(1,i,j+1)] + Hz1[ind(1,i,j-1)] -
									   4 * Hz1[ind(0,i,j)] - 4 * Hz1[ind(1,i,j)]);
			}
		}
		
		/* H -> H1 -> H2 */
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
		
		/* Data -> h5 */
		grp = fout.createGroup("/" + to_string(n));
		dset = fout.createDataSet("/"+to_string(n)+"/Ex", PredType::NATIVE_DOUBLE, dspace);
		dset.write(Ex, PredType::NATIVE_DOUBLE);
		dset = fout.createDataSet("/"+to_string(n)+"/Ey", PredType::NATIVE_DOUBLE, dspace);
		dset.write(Ey, PredType::NATIVE_DOUBLE);
		dset = fout.createDataSet("/"+to_string(n)+"/Ez", PredType::NATIVE_DOUBLE, dspace);
		dset.write(Ez, PredType::NATIVE_DOUBLE);
		dset = fout.createDataSet("/"+to_string(n)+"/Hx", PredType::NATIVE_DOUBLE, dspace);
		dset.write(Hx, PredType::NATIVE_DOUBLE);
		dset = fout.createDataSet("/"+to_string(n)+"/Hy", PredType::NATIVE_DOUBLE, dspace);
		dset.write(Hy, PredType::NATIVE_DOUBLE);
		dset = fout.createDataSet("/"+to_string(n)+"/Hz", PredType::NATIVE_DOUBLE, dspace);
		dset.write(Hz, PredType::NATIVE_DOUBLE);
		
		cout << 100.*n/iterations << " %" << endl;
	}
	}catch(FileIException error){
	error.printError();
	return -1;
	}catch(DataSetIException error){
	error.printError();
	return -1;
	}catch(DataSpaceIException error){
	error.printError();
	return -1;
	}
	
	return 0;
}
