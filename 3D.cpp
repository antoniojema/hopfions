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

void EdgeCondE();
void EdgeCondH();
double Edge(double i, double a, double b, double d, double W[], int p, bool p1, bool p2);

int main(){
	int n, i, j, k;
	double sigma, sigmam, Dx, Dy, Dz, Dt, epsilon, mu, c, L, x;
	
	L=10;
	sigma = sigmam = 0;
	epsilon = mu = 1;
	c=1;
	Dx = Dy = Dz = L/N;
	Dt = Dx/(2*c);
	
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
	hsize_t dimsf[1];
	dimsf[0] = (N+1)*(N+1)*(N+1);
	DataSpace dspace(1, dimsf);
	Group grp;
	DataSet dset;
	hsize_t dimsa[1] = { 1 };
	DataSpace attrds = DataSpace(1, dimsa);
	Attribute attr = fout.createAttribute("N",PredType::NATIVE_INT, attrds);
	//Attribute attr = dset.createAttribute("N",PredType::NATIVE_INT, attrds); 
	int attr_N[1] = {N};
	attr.write(PredType::NATIVE_INT, attr_N);
	attr = fout.createAttribute("iterations",PredType::NATIVE_INT, attrds);
	//attr = dset.createAttribute("iterations",PredType::NATIVE_INT, attrds);
	int attr_it[1] = {iterations};
	attr.write(PredType::NATIVE_INT, attr_it);
	
	for(n=0; n<=iterations; n++){
		/***** E *****/
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
					if (i>=4*N/10 && i<6*N/10 && j>=4*N/10 && j<6*N/10 && k>=4*N/10 && k<6*N/10 && n<=31){
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
		
		/* Edges conditions */
		EdgeCondE();
		
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
		
		/***** H *****/
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
		
		/* Edge conditions */
		EdgeCondH();
		
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
	
void EdgeCondE(){ //This only works for cubic grid
	double a = 0.25 * M_PI, b, d, f;
	for(int i=0; i<=N-1; i++){
		b = atan(1.*N/(sqrt(2)*abs(0.5*N-i)));
		d = sqrt(0.5*N*N + (0.5*N-i)*(0.5*N-i));
		f = sqrt((d-1.)/d);
		Ex[ind(i,0,0)] = Edge(i, a, b, f, Ex2, 1, 0, 0);
		Ey[ind(i,0,0)] = Edge(i, a, b, f, Ey2, 1, 0, 0);
		Ez[ind(i,0,0)] = Edge(i, a, b, f, Ez2, 1, 0, 0);
		Ex[ind(i,0,N)] = Edge(i, a, b, f, Ex2, 1, 0, N);
		Ey[ind(i,0,N)] = Edge(i, a, b, f, Ey2, 1, 0, N);
		Ez[ind(i,0,N)] = Edge(i, a, b, f, Ez2, 1, 0, N);
		Ex[ind(i,N,0)] = Edge(i, a, b, f, Ex2, 1, N, 0);
		Ey[ind(i,N,0)] = Edge(i, a, b, f, Ey2, 1, N, 0);
		Ez[ind(i,N,0)] = Edge(i, a, b, f, Ez2, 1, N, 0);
		Ex[ind(i,N,N)] = Edge(i, a, b, f, Ex2, 1, N, N);
		Ey[ind(i,N,N)] = Edge(i, a, b, f, Ey2, 1, N, N);
		Ez[ind(i,N,N)] = Edge(i, a, b, f, Ez2, 1, N, N);
		
		Ex[ind(0,i,0)] = Edge(i, a, b, f, Ex2, 2, 0, 0);
		Ey[ind(0,i,0)] = Edge(i, a, b, f, Ey2, 2, 0, 0);
		Ez[ind(0,i,0)] = Edge(i, a, b, f, Ez2, 2, 0, 0);
		Ex[ind(0,i,N)] = Edge(i, a, b, f, Ex2, 2, 0, N);
		Ey[ind(0,i,N)] = Edge(i, a, b, f, Ey2, 2, 0, N);
		Ez[ind(0,i,N)] = Edge(i, a, b, f, Ez2, 2, 0, N);
		Ex[ind(N,i,0)] = Edge(i, a, b, f, Ex2, 2, N, 0);
		Ey[ind(N,i,0)] = Edge(i, a, b, f, Ey2, 2, N, 0);
		Ez[ind(N,i,0)] = Edge(i, a, b, f, Ez2, 2, N, 0);
		Ex[ind(N,i,N)] = Edge(i, a, b, f, Ex2, 2, N, N);
		Ey[ind(N,i,N)] = Edge(i, a, b, f, Ey2, 2, N, N);
		Ez[ind(N,i,N)] = Edge(i, a, b, f, Ez2, 2, N, N);
		
		Ex[ind(0,0,i)] = Edge(i, a, b, f, Ex2, 3, 0, 0);
		Ey[ind(0,0,i)] = Edge(i, a, b, f, Ey2, 3, 0, 0);
		Ez[ind(0,0,i)] = Edge(i, a, b, f, Ez2, 3, 0, 0);
		Ex[ind(0,N,i)] = Edge(i, a, b, f, Ex2, 3, 0, N);
		Ey[ind(0,N,i)] = Edge(i, a, b, f, Ey2, 3, 0, N);
		Ez[ind(0,N,i)] = Edge(i, a, b, f, Ez2, 3, 0, N);
		Ex[ind(N,0,i)] = Edge(i, a, b, f, Ex2, 3, N, 0);
		Ey[ind(N,0,i)] = Edge(i, a, b, f, Ey2, 3, N, 0);
		Ez[ind(N,0,i)] = Edge(i, a, b, f, Ez2, 3, N, 0);
		Ex[ind(N,N,i)] = Edge(i, a, b, f, Ex2, 3, N, N);
		Ey[ind(N,N,i)] = Edge(i, a, b, f, Ey2, 3, N, N);
		Ez[ind(N,N,i)] = Edge(i, a, b, f, Ez2, 3, N, N);
	}
	
	return;
}

void EdgeCondH(){
	double a = 0.25 * M_PI, b, d, f;
	for(int i=0; i<=N-1; i++){
		b = atan(1.*N/(sqrt(2)*abs(0.5*N-i)));
		d = sqrt(0.5*N*N + (0.5*N-i)*(0.5*N-i));
		f = sqrt((d-1.)/d);
		Hx[ind(i,0,0)] = Edge(i, a, b, f, Hx2, 1, 0, 0);
		Hy[ind(i,0,0)] = Edge(i, a, b, f, Hy2, 1, 0, 0);
		Hz[ind(i,0,0)] = Edge(i, a, b, f, Hz2, 1, 0, 0);
		Hx[ind(i,0,N)] = Edge(i, a, b, f, Hx2, 1, 0, N);
		Hy[ind(i,0,N)] = Edge(i, a, b, f, Hy2, 1, 0, N);
		Hz[ind(i,0,N)] = Edge(i, a, b, f, Hz2, 1, 0, N);
		Hx[ind(i,N,0)] = Edge(i, a, b, f, Hx2, 1, N, 0);
		Hy[ind(i,N,0)] = Edge(i, a, b, f, Hy2, 1, N, 0);
		Hz[ind(i,N,0)] = Edge(i, a, b, f, Hz2, 1, N, 0);
		Hx[ind(i,N,N)] = Edge(i, a, b, f, Hx2, 1, N, N);
		Hy[ind(i,N,N)] = Edge(i, a, b, f, Hy2, 1, N, N);
		Hz[ind(i,N,N)] = Edge(i, a, b, f, Hz2, 1, N, N);
		
		Hx[ind(0,i,0)] = Edge(i, a, b, f, Hx2, 2, 0, 0);
		Hy[ind(0,i,0)] = Edge(i, a, b, f, Hy2, 2, 0, 0);
		Hz[ind(0,i,0)] = Edge(i, a, b, f, Hz2, 2, 0, 0);
		Hx[ind(0,i,N)] = Edge(i, a, b, f, Hx2, 2, 0, N);
		Hy[ind(0,i,N)] = Edge(i, a, b, f, Hy2, 2, 0, N);
		Hz[ind(0,i,N)] = Edge(i, a, b, f, Hz2, 2, 0, N);
		Hx[ind(N,i,0)] = Edge(i, a, b, f, Hx2, 2, N, 0);
		Hy[ind(N,i,0)] = Edge(i, a, b, f, Hy2, 2, N, 0);
		Hz[ind(N,i,0)] = Edge(i, a, b, f, Hz2, 2, N, 0);
		Hx[ind(N,i,N)] = Edge(i, a, b, f, Hx2, 2, N, N);
		Hy[ind(N,i,N)] = Edge(i, a, b, f, Hy2, 2, N, N);
		Hz[ind(N,i,N)] = Edge(i, a, b, f, Hz2, 2, N, N);
		
		Hx[ind(0,0,i)] = Edge(i, a, b, f, Hx2, 3, 0, 0);
		Hy[ind(0,0,i)] = Edge(i, a, b, f, Hy2, 3, 0, 0);
		Hz[ind(0,0,i)] = Edge(i, a, b, f, Hz2, 3, 0, 0);
		Hx[ind(0,N,i)] = Edge(i, a, b, f, Hx2, 3, 0, N);
		Hy[ind(0,N,i)] = Edge(i, a, b, f, Hy2, 3, 0, N);
		Hz[ind(0,N,i)] = Edge(i, a, b, f, Hz2, 3, 0, N);
		Hx[ind(N,0,i)] = Edge(i, a, b, f, Hx2, 3, N, 0);
		Hy[ind(N,0,i)] = Edge(i, a, b, f, Hy2, 3, N, 0);
		Hz[ind(N,0,i)] = Edge(i, a, b, f, Hz2, 3, N, 0);
		Hx[ind(N,N,i)] = Edge(i, a, b, f, Hx2, 3, N, N);
		Hy[ind(N,N,i)] = Edge(i, a, b, f, Hy2, 3, N, N);
		Hz[ind(N,N,i)] = Edge(i, a, b, f, Hz2, 3, N, N);
	}
	
	return;
}

double Edge(double i, double a, double b, double f, double W[], int p, bool p1, bool p2){
	/********************/
	if(p==1){
		if(p1==false && p2==false){
			if(i<1.*N/2) return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(i,0,0)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(i+1,0,0)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(i,1,0)]		+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(i+1,1,0)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(i,0,1)]		+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(i+1,0,1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(i,1,1)]		+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(i+1,1,1)] );
			
			else		 return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(i,0,0)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(i-1,0,0)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(i,1,0)]		+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(i-1,1,0)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(i,0,1)]		+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(i-1,0,1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(i,1,1)]		+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(i-1,1,1)] );
		}else if(p1==false && p2==true){
			if(i<1.*N/2) return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(i,0,N)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(i+1,0,N)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(i,1,N)]		+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(i+1,1,N)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(i,0,N-1)]	+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(i+1,0,N-1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(i,1,N-1)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(i+1,1,N-1)] );
			
			else		 return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(i,0,N)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(i-1,0,N)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(i,1,N)]		+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(i-1,1,N)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(i,0,N-1)]	+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(i-1,0,N-1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(i,1,N-1)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(i-1,1,N-1)] );
		}else if(p1==true && p2==false){
			if(i<1.*N/2) return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(i,N,0)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(i+1,N,0)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(i,N-1,0)]	+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(i+1,N-1,0)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(i,N,1)]		+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(i+1,N,1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(i,N-1,1)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(i+1,N-1,1)] );
			
			else		 return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(i,N,0)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(i-1,N,0)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(i,N-1,0)]	+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(i-1,N-1,0)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(i,N,1)]		+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(i-1,N,1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(i,N-1,1)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(i-1,N-1,1)] );
		}else{
			if(i<1.*N/2) return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(i,N,N)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(i+1,N,N)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(i,N-1,N)]	+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(i+1,N-1,N)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(i,N,N-1)]	+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(i+1,N,N-1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(i,N-1,N-1)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(i+1,N-1,N-1)] );
			
			else		 return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(i,N,N)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(i-1,N,N)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(i,N-1,N)]	+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(i-1,N-1,N)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(i,N,N-1)]	+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(i-1,N,N-1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(i,N-1,N-1)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(i-1,N-1,N-1)] );
		}
	/********************/
	}else if(p==2){
		if(p1==false && p2==false){
			if(i<1.*N/2) return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(0,i,0)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(0,i+1,0)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(1,i,0)]		+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(1,i+1,0)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(0,i,1)]		+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(0,i+1,1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(1,i,1)]		+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(1,i+1,1)] );
			
			else		 return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(0,i,0)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(0,i-1,0)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(1,i,0)]		+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(1,i-1,0)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(0,i,1)]		+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(0,i-1,1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(1,i,1)]		+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(1,i-1,1)] );
		}else if(p1==false && p2==true){
			if(i<1.*N/2) return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(0,i,N)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(0,i+1,N)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(1,i,N)]		+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(1,i+1,N)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(0,i,N-1)]	+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(0,i+1,N-1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(1,i,N-1)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(1,i+1,N-1)] );
			
			else		 return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(0,i,N)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(0,i-1,N)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(1,i,N)]		+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(1,i-1,N)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(0,i,N-1)]	+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(0,i-1,N-1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(1,i,N-1)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(1,i-1,N-1)] );
		}else if(p1==true && p2==false){
			if(i<1.*N/2) return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(N,i,0)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(N,i+1,0)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(N-1,i,0)]	+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(N-1,i+1,0)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(N,i,1)]		+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(N,i+1,1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(N-1,i,1)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(N-1,i+1,1)] );
			
			else		 return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(N,i,0)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(N,i-1,0)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(N-1,i,0)]	+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(N-1,i-1,0)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(N,i,1)]		+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(N,i-1,1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(N-1,i,1)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(N-1,i-1,1)] );
		}else{
			if(i<1.*N/2) return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(N,i,N)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(N,i+1,N)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(N-1,i,N)]	+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(N-1,i+1,N)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(N,i,N-1)]	+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(N,i+1,N-1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(N-1,i,N-1)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(N-1,i+1,N-1)] );
			
			else		 return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(N,i,N)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(N,i-1,N)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(N-1,i,N)]	+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(N-1,i-1,N)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(N,i,N-1)]	+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(N,i-1,N-1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(N-1,i,N-1)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(N-1,i-1,N-1)] );
		}
	/********************/
	}else if(p==3){
		if(p1==false && p2==false){
			if(i<1.*N/2) return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(0,0,i)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(0,0,i+1)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(0,1,i)]		+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(0,1,i+1)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(1,0,i)]		+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(1,0,i+1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(1,1,i)]		+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(1,1,i+1)] );
			
			else		 return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(0,0,i)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(0,0,i-1)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(0,1,i)]		+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(0,1,i-1)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(1,0,i)]		+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(1,0,i-1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(1,1,i)]		+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(1,1,i-1)] );
		}else if(p1==false && p2==true){
			if(i<1.*N/2) return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(0,N,i)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(0,N,i+1)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(0,N-1,i)]	+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(0,N-1,i+1)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(1,N,i)]		+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(1,N,i+1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(1,N-1,i)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(1,N-1,i+1)] );
			
			else		 return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(0,N,i)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(0,N,i-1)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(0,N-1,i)]		+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(0,N-1,i-1)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(1,N,i)]		+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(1,N,i-1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(1,N-1,i)]		+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(1,N-1,i-1)] );
		}else if(p1==true && p2==false){
			if(i<1.*N/2) return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(N,0,i)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(N,0,i+1)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(N,1,i)]		+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(N,1,i+1)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(N-1,0,i)]	+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(N-1,0,i+1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(N-1,1,i)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(N-1,1,i+1)] );
			
			else		 return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(N,0,i)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(N,0,i-1)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(N,1,i)]		+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(N,1,i-1)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(N-1,0,i)]	+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(N-1,0,i-1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(N-1,1,i)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(N-1,1,i-1)] );
		}else{
			if(i<1.*N/2) return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(N,N,i)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(N,N,i+1)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(N,N-1,i)]	+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(N,N-1,i+1)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(N-1,N,i)]	+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(N-1,N,i+1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(N-1,N-1,i)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(N-1,N-1,i+1)] );
			
			else		 return f * (
						 (1-sin(b))*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))	* W[ind(N,N,i)]		+
						 (1-sin(b))*(1-cos(b)*sin(a))*cos(b)*cos(a)		* W[ind(N,N,i-1)]	+
						 (1-sin(b))*cos(b)*sin(a)*(1-cos(b)*cos(a))		* W[ind(N,N-1,i)]	+
						 (1-sin(b))*cos(b)*cos(b)*sin(a)*cos(a)			* W[ind(N,N-1,i-1)]	+
						 sin(b)*(1-cos(b)*sin(a))*(1-cos(b)*cos(a))		* W[ind(N-1,N,i)]	+
						 sin(b)*(1-cos(b)*sin(a))*cos(b)*cos(a)			* W[ind(N-1,N,i-1)]	+
						 sin(b)*cos(b)*(1-cos(b)*cos(a))				* W[ind(N-1,N-1,i)]	+
						 sin(b)*cos(b)*cos(b)*sin(a)*cos(a)				* W[ind(N-1,N-1,i-1)] );
		}
	}
	return 0;
}
