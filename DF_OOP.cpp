#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

typedef std::vector<double> vct;
typedef std::vector<vct> matrix;




//Produit matrice vecteur d'une matrice pentadiagonale

vct matvec(double a, double b, double c, int Nx, int Ny, vector<double> x)  {
      int n = Nx*Ny;
      vct y;
      y.resize(n);

      // cas 1
      y[0] = a*x[0] + b*x[1] + c*x[Nx];
      // cas 2
      for (int s = 1; s < Nx-1; s++) {
        y[s] = b*(x[s-1] + x[s+1]) + a*x[s] + c*x[Nx+s];
      }
      // cas 3
      y[Nx-1] = b*x[Nx-2] + a*x[Nx-1] + c*x[2*Nx-1];
      // cas 4 et 5 et 6
      for (int k=1; k < Ny-1; k++) {
        // cas 4
        y[k*Nx] = c*(x[(k-1)*Nx] + x[(k+1)*Nx]) + a*x[k*Nx] + b*x[k*Nx+1];
        // cas 5
        for (int s=k*Nx+1; s < (k+1)*Nx - 1; s++) {
          y[s] = c*(x[s-Nx] + x[s+Nx]) + a*x[s] + b*(x[s-1] + x[s+1]);
        }
        // cas 6
        y[(k+1)*Nx-1] = c*(x[k*Nx-1] + x[(k+2)*Nx-1]) + a*x[(k+1)*Nx-1] + b*x[(k+1)*Nx-2];
      }
      // cas 7
      y[Nx*(Ny-1)] = c*x[Nx*(Ny-2)] + a*x[Nx*(Ny-1)] + b*x[Nx*(Ny-1)+1];
      // cas 8
      for (int s = Nx*(Ny-1)+1; s < Nx*Ny-1; s++) {
        y[s] = b*(x[s-1] + x[s+1]) + a*x[s] + c*x[s-Nx];
      }
      // cas 9
      y[Nx*Ny-1] = c*x[Nx*(Ny-1)-1] + a*x[Nx*Ny-1] + b*x[Nx*Ny-2];
      // sortie
      return y;
    }


double sc(vct v, vct w) {    // produit scalaire de deux vecteurs
  double result = 0.0;
  for (int i = 0; i < v.size(); i++) {
      result += v[i] * w[i];
  }
  return result;
}

vct hm(double a, vct x) {    // scalaire * vecteur
  int n = x.size();
  vct y(n);
  for (int i = 0; i < n; i++) {
    y[i] = a*x[i];
  }
  return y;
}

vct add(vct x, vct y) {   // somme de deux vecteurs, pour faire une soustractione on appelle add(x, -y)
  vct z(x.size());
  for (int i = 0; i < x.size(); i++) {
    z[i] = x[i] + y[i];
  }
  return z;
}

vct GCON(double a, double b, double c, int Nx, int Ny, vct Y, vector<double> x0, double epsilon, int kmax){   // solves AX = Y
    int k,n; double beta, alpha, gamma;
    
    n= Y.size();
    x0.resize(n);  //
    vct r0, p,z,x,r ;
    cout << "test" << endl;
    r0= add(Y,hm(-1,matvec(a, b, c, Nx, Ny, x0)));            // does matvec work?
    p=r0;
    beta= sqrt(sc(r0,r0));
    k=0;
    while ( beta>epsilon && k<=kmax)  {
        z=matvec(a,b,c,Nx,Ny,p);
        alpha=sc(r0,r0)/sc(z,p);
        x= add(x0,hm(alpha,p));
        r=add(r0,hm(-alpha,z));
        gamma=sc(r,r)/sc(r0,r0);
        p=add(r,hm(gamma,p));
        beta=sqrt(sc(r,r));

        k=k+1;
        r0=r;
        x0=x   ;
    }
               
    return x;
}
  

// functions of boundaries and sources
double f(double x, double y, double t, int cas) {
  double L=1.0;
  switch (cas) {
  case 1:
    return 2*(y-y*y+x-x*x);

  case 2:
    return sin(x)+cos(y);

  case 3:
    return exp(-pow(x-L*x/2,2)-pow(y-L*y/2,2))*cos(M_PI*t/2);

  default:
    return 0.0;
  }
}

double g(double x, double y, double t, int cas) {
  switch (cas) {
  case 1:
    return 0.0;
  case 2:
    return sin(x)+cos(y);
  case 3:
    return 0.0;
  default:
    return 0.0;
  }
}

double h(double x, double y, double t, int cas) {
  switch (cas) {
  case 1:
    return 0.0;
  case 2:
    return sin(x)+cos(y);
  case 3:
    return 1.0;
  default:
    return 0.0;
  }
}

// function to calculate the Laplace operator
/*double laplacian(const std::vector<double>& u, int i, int j) {
    double value = 0.0;
    if (i > 0) value += u[(i-1)*N + j];
    if (i < N-1) value += u[(i+1)*N + j];
    if (j > 0) value += u[i*N + (j-1)];
    if (j < N-1) value += u[i*N + (j+1)];
    value -= 4.0 * u[i*N + j];
    return value / (dx*dx + dy*dy);
}

*/


int main() {

int Nx, Ny;

Nx =5; Ny= 5;
double dx, dy, dt(1.), Lx(1.), Ly(1.), D(1.);

int tmax(10);
dx = Lx/(Nx+1);
dy = Ly/(Ny+1);

double a, b,c;
b = -1./dx*dx;
c = -1./dy*dy;
a = -2*(b+c);

int k;
// définit le second membre 
vector<double> F;
F = hm(0,F);      // sans terme source
F.resize(Nx*Ny);


// gauche
  for(k=0;k<=Ny-1; k++) F[k*Nx] =  F[k*Nx] + D*100./(dx*dx);  // 100 ici est la condition de dirichlet
    

// droit
  for(k=1;k<=Ny; k++) F[k*Nx- 1 ]= F[k*Nx- 1] + D*100./(dx*dx);
     

// bas
for(k=0;k<Nx; k++) F[k]=F[k] + D*100./(dy*dy);
  

// haut
for(k=Nx*Ny - Nx; k<=Nx*Ny-1; k++) F[k]= F[k] + D*100./(dy*dy);


// initialisation du vecteur solution T
vector<double> T;
T.resize(Nx*Ny);
T = hm(0,T);
int m(0);
for(int j=0;j<Ny; j++){
  for(int i=0; i<Nx; i++){
    T[j*Nx + i] = 1;
  }
        
}

int kmax; 
double epsilon;
kmax=1000;
epsilon=pow(10,-5);

int time(0);

vector<double> temp,x0;
temp.resize(Nx*Ny); x0.resize(Nx*Ny);
temp = hm(0,temp); x0 = hm(0,x0);
//add(T,hm(dt,F))
  while (time<=tmax){
    T = GCON(a, b, c, Nx, Ny,temp,x0,epsilon,kmax);       // the issue comes from here 
    
    time= time + dt;
    cout<< "----------------------------------------------" << endl; 
     if(time == tmax) cout << "end of loop" << endl;
     cout << 1+1. << endl;
  }

cout << "La solution du système est : ";
for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
	     std::cout << T[i*Nx + j] << " ";
      }
      std::cout << std::endl;
  }





return 0;
}
