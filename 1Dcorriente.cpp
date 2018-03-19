#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

const int N=2000;
double A, B, Ze, Zm;
double Ex[N+1], Hy[N+1];

int main(){
	int i, j;
	double Ex1[2], Ex2[2], Hy1[2], Hy2[2];
	double sigma, sigmam, Dz, Dt, epsilon, mu, c, L;
	ofstream fout;
	
	L=10;
	sigma = 0; sigmam = 0;
	Dz = Dt = L/N;
	epsilon = mu = 1;
	c=1;
	
	A = (1.-sigma*Dt/(2.*epsilon))/(1.+sigma*Dt/(2.*epsilon));
	Ze = (Dt/(epsilon*Dz))/(1.+sigma*Dt/(2.*epsilon));
	B = (1.-sigmam*Dt/(2.*mu))/(1.+sigmam*Dt/(2.*mu));
	Zm = (Dt/(mu*Dz))/(1.+sigmam*Dt/(2.*mu));
	
	fout.open("1Dalt.txt");
	
	//Condiciones iniciales
	for(j=0; j<=N; j++){
		Ex[j] = Hy[j]= 0;
		fout << (10./N)*j << "	" << Ex[j] << "	" << Hy[j] << endl;
	}
	for(j=0; j<=1; j++){
		Ex1[j] = Ex2[j] = Hy1[j] = Hy2[j] = 0;
	}
		
	for(i=0; i<=10000; i++){
		//Calculo el nuevo Ex a partir del viejo y de los Hy
		for(j=0; j<=N-1; j++){
			Ex[j] = A*Ex[j] + Ze*(Hy[j]-Hy[j+1]) - 0.01*exp(-0.001*(j-1000)*(j-1000))*sin(0.01*i);
			//if(j==1000) Ex[j] -= sin(0.01*i);
		}
		
		//CONDICION MUR:
		Ex[N] = -1.*Ex2[1] + ((c*Dt-Dz)/(c*Dt+Dz))*(Ex[N-1]+Ex2[0]) + (2.*Dz/(c*Dt+Dz))*(Ex1[0]+Ex1[1]);
		
		//Ex->Ex1->Ex2
		Ex2[0]=Ex1[0];
		Ex1[0]=Ex[N];
		Ex2[1]=Ex1[1];
		Ex1[1]=Ex[N-1];
		
		//Hago lo propio con Hy
		for(j=1; j<=N; j++){
			Hy[j] = B*Hy[j] + Zm*(Ex[j-1]-Ex[j]);
		}
		
		//CONDICION MUR:
		Hy[0] = -1.*Hy2[1] + ((c*Dt-Dz)/(c*Dt+Dz))*(Hy[1]+Hy2[0]) + (2.*Dz/(c*Dt+Dz))*(Hy1[0]+Hy1[1]);
		
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
