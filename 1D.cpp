#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double A[2], Ze[2];
double B[2], Zm[2];

const int N=2000;
double Ex[N+1], Hy[N+1];

int main(){
	int i, j;
	double Ex1[2], Ex2[2], Hy1[2], Hy2[2];
	double sigma[2], sigmam[2], Dz, Dt, epsilon[2], mu[2], c[2], L, x;
	double front=6;
	int indfront=(N/10)*front;
	ofstream fout;
	
	L=10;
	sigma[0] = 0; sigma[1] = 1;
	sigmam[0] = sigmam[1] = 0;
	Dz = Dt = L/N;
	epsilon[0] = 1; epsilon[1] = 2;
	c[0]=1; c[1]=1/*/(sqrt(epsilon[1]*mu[1]))*/;
	mu[0] = mu[1]=1;
	
	for(i=0;i<=1;i++){
		A[i] = (1.-sigma[i]*Dt/(2.*epsilon[i]))/(1.+sigma[i]*Dt/(2.*epsilon[i]));
		Ze[i] = (Dt/(epsilon[i]*Dz))/(1.+sigma[i]*Dt/(2.*epsilon[i]));
		B[i] = (1.-sigmam[i]*Dt/(2.*mu[i]))/(1.+sigmam[i]*Dt/(2.*mu[i]));
		Zm[i] = (Dt/(mu[i]*Dz))/(1.+sigmam[i]*Dt/(2.*mu[i]));
	}
	
	fout.open("1D.txt");
	
	//Condiciones iniciales
	for(j=0; j<=N; j++){
		x=(10./N)*j;
		Ex[j]=f(2*(x-4.));
		Hy[j]=f(2*(x-4.));
		fout << x << "	" << Ex[j] << "	" << Hy[j] << endl;
	}
	for(j=0; j<=1; j++){
		Ex1[j]=Ex2[j]=Hy1[j]=Hy2[j]=0;
	}
		
	for(i=0; i<=10000; i++){
		//Calculo el nuevo Ex a partir del viejo y de los Hy
		for(j=0; j<indfront; j++){
			Ex[j] = A[0]*Ex[j]+Ze[0]*(Hy[j]-Hy[j+1]);
		}
		for(j=indfront; j<=N-1; j++){
			Ex[j] = A[1]*Ex[j]+Ze[1]*(Hy[j]-Hy[j+1]);
		}
		//CONDICION PERIODICA:
		//Ex[N] = A[1]*Ex[N]+Ze[1]*(Hy[N]-Hy[0]);
		
		//CONDICION MUR (la primera es de internete, la segunda es del libro):
		//Ex[N] = Ex1[1]+((c[1]*Dt-Dz)/(c[1]*Dt+Dz))*(Ex[N-1]-Ex1[0]);
		Ex[N] = -1.*Ex2[1] + ((c[1]*Dt-Dz)/(c[1]*Dt+Dz))*(Ex[N-1]+Ex2[0]) + (2.*Dz/(c[1]*Dt+Dz))*(Ex1[0]+Ex1[1]);
		
		//Ex->Ex1->Ex2
		Ex2[0]=Ex1[0];
		Ex1[0]=Ex[N];
		Ex2[1]=Ex1[1];
		Ex1[1]=Ex[N-1];
		
		//Hago lo propio con Hy
		for(j=1; j<indfront; j++){
			Hy[j] = B[0]*Hy[j]+Zm[0]*(Ex[j-1]-Ex[j]);
		}
		for(j=indfront; j<=N; j++){
			Hy[j] = B[1]*Hy[j]+Zm[1]*(Ex[j-1]-Ex[j]);
		}
		//CONDICION PERIODICA:
		//Hy[0] = B[0]*Hy[0]+Zm[0]*(Ex[N]-Ex[0]);
		
		//CONDICION MUR:
		//Hy[0] = Hy1[1]+((c[0]*Dt-Dz)/(c[0]*Dt+Dz))*(Hy[1]-Hy1[0]);
		Hy[0] = -1.*Hy2[1] + ((c[0]*Dt-Dz)/(c[0]*Dt+Dz))*(Hy[1]+Hy2[0]) + (2.*Dz/(c[0]*Dt+Dz))*(Hy1[0]+Hy1[1]);
		
		//Hy->Hy1->Hy2
		Hy2[0]=Hy1[0];
		Hy1[0]=Hy[0];
		Hy2[1]=Hy1[1];
		Hy1[1]=Hy[1];
		
		//Ploteo los datos en el txt
		if(i%10==0){
			for(j=0; j<=N; j++){
				//QUIERO HACERLO EN H5 ASI QUE LO CAMBIARE (o no)
				fout << (10./N)*j << "	" << Ex[j] << "	" << Hy[j] << endl;
			}
		}
	}
	fout.close();
	
	return 0;
}
