#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin);

void step_one(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx, const int N);

void step_two(cmplx* const f, const double dt, const double dx, const int N);

//-----------------------------------
int main(){

	const int Nx = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (Nx - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double eta = 0.2;
	double t = 0;

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h;
	
	init(psi0, eta, dx, dt,Nx);

	writeToFile(psi0,"psi_0", dx,Nx,xmin);


	for (int i = 1; i <= Na; i++) {
	  
// 	    Strang-Splitting
 	  step_one(psi1, psi0, dt/2.0, dx, Nx);
	  
		for (int j = 1; j <= Nk-1; j++) {
		 
		  step_two(psi1, dt, dx, Nx);
		  
		  h = psi0;
		  psi0 = psi1;
		  psi1 = h;
		  
		  
		  step_one(psi1, psi0, dt, dx, Nx);
		  	  
		  
		  t += dt;
		  
		}
 	  step_two(psi1, dt, dx, Nx);
	  h = psi0;
	  psi0 = psi1;
	  psi1 = h;
	  step_one(psi1, psi0, dt/2.0, dx, Nx);
	  
	  h = psi0;
	  psi0 = psi1;
	  psi1 = h;
	  
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin);
	}
	  
	delete[] psi0;
	delete[] psi1;
	  
	return 0;
}
//-----------------------------------

// Linear part
void step_one(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx, const int N){
  
  cmplx* d=new cmplx[N];
  cmplx* u=new cmplx[N];
  cmplx* l=new cmplx[N];

  for(int i=0;i<N;i++) d[i] = 1.0 - 2.0*cmplx(0.0, dt/(dx*dx));
  for(int i=0;i<N;i++) u[i] = cmplx(0.0, dt/(dx*dx));
  for(int i=0;i<N;i++) l[i] = cmplx(0.0, dt/(dx*dx));

  
  // Forward
  
  for(int i=1;i<N;i++){ 
    d[i] = d[i] - l[i]*u[i-1]/d[i-1]; // einzig d0 bleibt unverändert also starten wir bei i=1
    f0[i] = f0[i] - l[i]*f0[i-1]/d[i-1];
  }
  // Backward
  
  f1[N-1] = f0[N-1] / d[N-1];
  
  for(int i=N-1;i>=1;i--) 
    f1[i-1] = (f0[i-1]-u[i-1]*f1[i])/d[i-1];
  
  
  delete[] d;
  delete[] u;
  delete[] l;
  
}
// Non-Linear part

void step_two(cmplx* const f, const double dt, const double dx, const int N){
  
  for(int i=0;i<N;i++){
    f[i] = f[i]*exp(-cmplx(0.0, dt*norm(f[i])));
  }
  
}


void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] = 2*f/cosh(eta * x);
	}
}
