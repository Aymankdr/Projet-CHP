#include <iostream>
#include <vector>
#include <cmath>

typedef std::vector<double> vct;
typedef std::vector<vct> matrix;

// Laplacian
class Pentadiagonal {
private:
  double _alpha, _beta, _gamma;
  int _Nx, _Ny;
public:
  Pentadiagonal(double alpha, double beta, double gamma, int Nx, int Ny) {
  _alpha = alpha; _beta = beta; _gamma = gamma; _Nx = Nx; _Ny = Ny;
  }
  vct matmul(const vct x) const {
      int n = _Nx*_Ny;
      vct y(n);
      double a = _alpha, b = _beta, c = _gamma;
      // cas 1
      y[0] = a*x[0] + b*x[1] + c*x[_Nx];
      // cas 2
      for (int s = 1; s < _Nx-1; s++) {
        y[s] = b*(x[s-1] + x[s+1]) + a*x[s] + c*x[_Nx+s];
      }
      // cas 3
      y[_Nx-1] = b*x[_Nx-2] + a*x[_Nx-1] + c*x[2*_Nx-1];
      // cas 4 et 5 et 6
      for (int k=1; k < _Ny-1; k++) {
        // cas 4
        y[k*_Nx] = c*(x[(k-1)*_Nx] + x[(k+1)*_Nx]) + a*x[k*_Nx] + b*x[k*_Nx+1];
        // cas 5
        for (int s=k*_Nx+1; s < (k+1)*_Nx - 1; s++) {
          y[s] = c*(x[s-_Nx] + x[s+_Nx]) + a*x[s] + b*(x[s-1] + x[s+1]);
        }
        // cas 6
        y[(k+1)*_Nx-1] = c*(x[k*_Nx-1] + x[(k+2)*_Nx-1]) + a*x[(k+1)*_Nx-1] + b*x[(k+1)*_Nx-2];
      }
      // cas 7
      y[_Nx*(_Ny-1)] = c*x[_Nx*(_Ny-2)] + a*x[_Nx*(_Ny-1)] + b*x[_Nx*(_Ny-1)+1];
      // cas 8
      for (int s = _Nx*(_Ny-1)+1; s < _Nx*_Ny-1; s++) {
        y[s] = b*(x[s-1] + x[s+1]) + a*x[s] + c*x[s-_Nx];
      }
      // cas 9
      y[_Nx*_Ny-1] = c*x[_Nx*(_Ny-1)-1] + a*x[_Nx*_Ny-1] + b*x[_Nx*_Ny-2];
      // sortie
      return y;
    }
    int get_size() {return _Nx*_Ny;}
  // fin public
};

double sc(vct v, vct w) {
  double result = 0.0;
  for (int i = 0; i < v.size(); i++) {
      result += v[i] * w[i];
  }
  return result;
}

vct hm(double a, vct x) {
  int n = x.size();
  vct y(n);
  for (int i = 0; i < n; i++) {
    y[i] = a*x[i];
  }
  return y;
}

vct add(vct x, vct y) {
  vct z(x.size());
  for (int i = 0; i < x.size(); i++) {
    z[i] = x[i] + y[i];
  }
  return z;
}

vct conjugateGradientMethod(Pentadiagonal* A, vct b, vct x0, double tol, int maxIter) {
    int n = A->get_size();
    vct x(n, 0.0);
    vct r(n, 0.0);
    vct d(n, 0.0);
    double alpha = 0.0;
    double beta = 0.0;
    double rnorm = 0.0;
    double rnorm0 = 0.0;

    x = x0;
    r = b;
    d = r;
    rnorm = sqrt(sc(r, r));
    rnorm0 = rnorm;

    int k = 0;
    while (k < maxIter && rnorm / rnorm0 > tol) {
        alpha = sc(r, r) / sc(d, A->matmul(d));
        x = add( x , hm(alpha,d) );
        r = add( x , hm(-alpha, A->matmul(d)) );
        beta = sc(r, r) / sc(d, A->matmul(d));
        d = add( r , hm(beta,d) );
        rnorm = sqrt(sc(r, r));
        k++;
    }

    return x;
}

vct conjugate_gradient(Pentadiagonal A, vct y, vct x) {
    int n = A.get_size();
    vct r(n);
    vct p(n);
    vct Ap(n);
    vct alpha_beta(2);

    // Calculate the initial residual and search direction
    r = y;
    Ap = A.matmul(p);
    for (int i = 0; i < n; i++) {
        r[i] = y[i] - Ap[i];
        p[i] = r[i];
    }

    for (int k = 0; k < n; k++) {
        // Compute alpha_k
        Ap = A.matmul(p);
        alpha_beta[0] = 0;
        alpha_beta[1] = 0;
        for (int i = 0; i < n; i++) {
            alpha_beta[0] += r[i] * r[i];
            alpha_beta[1] += p[i] * Ap[i];
        }
        double alpha_k = alpha_beta[0] / alpha_beta[1];

        // Update solution and residual vectors
        for (int i = 0; i < n; i++) {
            x[i] += alpha_k * p[i];
            r[i] -= alpha_k * Ap[i];
        }

        // Check for convergence
        double r_norm = 0;
        for (int i = 0; i < n; i++) {
            r_norm += r[i] * r[i];
        }
        if (sqrt(r_norm) < 1e-6) {
            break;
        }

        // Compute beta_k
        alpha_beta[0] = 0;
        alpha_beta[1] = 0;
        vct Ap_old = Ap;
        Ap = A.matmul(p);
        for (int i = 0; i < n; i++) {
            alpha_beta[0] += r[i] * r[i];
            alpha_beta[1] += Ap[i] * Ap[i];
            p[i] = r[i] - alpha_k * Ap_old[i];
        }
        double beta_k = alpha_beta[0] / alpha_beta[1];

        // Update search direction
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta_k * p[i];
        }
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
/*
double laplacian(const vct& u, int i, int j) {
    double value = 0.0;
    if (i > 0) value += u[(i-1)*Nx + j];
    if (i < Nx-1) value += u[(i+1)*Nx + j];
    if (j > 0) value += u[i*Nx + (j-1)];
    if (j < Ny-1) value += u[i*Nx + (j+1)];
    value -= 4.0 * u[i*Nx + j];
    return value / (dx*dx + dy*dy);
} */

int main() {
  // Constants
  const int Nx = 100; // number of grid points in x direction
  const int Ny = 100; // number of grid points in y direction
  const int M  = 100; // number of time steps
  const double D = 1.0; // thermal diffusivity
  const double dt = 0.001; // time step size
  const double dx = 1.0 / (Nx + 1); // grid spacing
  const double dy = 1.0 / (Ny + 1); // grid spacing
  const double t0 = 0.0; // first instant

  // other constants
  const int maxIter = 500;
  const double tol  = 0.001;

  // Case
  int cas = 1;

  // initialize the temperature vector
  vct u(Nx*Ny, 0.0), F(Nx*Ny, 0.0), u_new(Nx*Ny);

  // set the boundary conditions
  for (int i = 0; i < Nx; i++) {
    u[i] = g(i*dx, 0.0, 0.0, cas); // bottom edge
    u[(Ny-1)*Nx + i] = g(i*dx, 1.0, 0.0, cas); // top edge
    u[i*Ny] = h(0.0, i*dy, 0.0,cas); // left edge
    u[i*Ny + Nx-1] = h(1.0, i*dy, 0.0, cas); // right edge
  }

  // New constants
  const double alpha = 2*(1./pow(dx,2)+1./pow(dy,2));
  const double beta  = -1./pow(dx,2);
  const double gamma = -1./pow(dy,2);


  // create the matrix A and B
  Pentadiagonal* A = new Pentadiagonal(alpha, beta, gamma, Nx, Ny);
  Pentadiagonal* B = new Pentadiagonal(1.0 + dt*alpha, dt*beta, dt*gamma, Nx, Ny);

  /*
  // time loop
  for (int k = 1; k <= M; k++) {
    // create the right-hand side vector F
    for (int i = 1; i < Nx-1; i++) {
       for (int j = 1; j < Ny-1; j++) {
    	    int index = i*Nx + j;
    	    F[index] = u[index] + dt*f(i*dx, j*dy, k*dt, cas);
    	    F[index] += D*dt * laplacian(u, i, j);
	     }
    }
  } */


  // solve the linear system AU=F
  double t = t0;
  for (int k = 1; k <= M; k++) {
    // adding timestep
    t  += dt;
    // generate F^k+1 without boundary conditions
    for (int i = 1; i < Nx; i++) {
    	for (int j = 1; j < Ny; j++) {
		F[i*Nx+j] = f(i*dx, j*dy, t, cas);
	}
    }
    // adding to F^k+1 boundary conditions
    for (int j=0; j <= Ny; j++) { // verticals 
	    
    }
    
    // resolution
    u_new = conjugateGradientMethod(B, add(u,hm(dt,F)), u, tol, maxIter); // u + dt*F^n+1
    u     = u_new;
  }


  // print the final temperature values
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
	     std::cout << u[i*Nx + j] << " ";
      }
      std::cout << std::endl;
  }

  return 0;
}
