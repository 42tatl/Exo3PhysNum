#include <iostream>
#include <fstream>
#include <iomanip>
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <valarray>
#include <cmath> 
#include <numeric>
using namespace std;

class Exercice3
{

private:
  double t, dt, tFin;
  double x0,a;
  double GM=6.674e-11;
  double mtot=0e0;
  double Omega;
  double mJ;
  double mS;
  valarray<double> v0 = std::valarray<double>(0.e0, 2);
  int N_excit, nsteps;
  int sampling;
  int last;
  int  nsel_physics;
  bool adapt;
  double alpha = 0e0;
  double beta = 0e0;
  double tol= 0e0;
  valarray<double> y0 = std::valarray<double>(0.e0, 4); // Correctly initialized
  valarray<double> y  = std::valarray<double>(0.e0, 4); // Correctly initialized
  ofstream *outputFile;

  void printOut(bool write)
  {
    if((!write && last>=sampling) || (write && last!=1))
    {
      double Energy = compute_energy(y[0],y[1],y[2],y[3]);
      *outputFile << t << " "<< y[0] << " " << y[1] << " "<< y[2] << " " << y[3] << " " \
      << Energy<< " "<< nsteps<< endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }

  std::valarray<double> get_f(double t, const std::valarray<double>& y) {
    std::valarray<double> ydot(0.0, 4);
    //TO DO
    if (nsel_physics == 1) { //two body problem 
      
      double r = sqrt(pow(y[2],2)+ pow(y[3],2)); //y = (vx, vy, x, y) 
    
      ydot[0] = -(GM*mS*y[2])/pow(r,3); 
      ydot[1] = -(GM*mS*y[3])/pow(r,3);
    }
    else if (nsel_physics == 2) { //reduced three body problem
      double x_prime_S = -alpha*a; 
      double x_prime_J = beta*a;

      double rS = sqrt( pow(y[2] - x_prime_S,2)+pow( y[3],2) );
      double rJ = sqrt( pow(y[2] - x_prime_J,2)+pow( y[3],2) ); //y = (vx, vy, x, y)
            
      ydot[0] = (-GM * mS * (y[2] - x_prime_S) ) /pow(rS, 3) - (GM * mJ *(y[2] - x_prime_J))/pow(rJ, 3) +pow(Omega,2) * y[2] + 2*Omega*y[1];
      ydot[1] = - (GM * y[3] * mS) /pow(rJ, 3) - (GM * mJ * y[3]) /pow(rJ, 3) + pow(Omega,2)*y[3] - 2 * Omega * y[0];
    }
    else{
        cerr << "No dynamics corresponds to this index" << endl;
        return ydot;
    }

    // Velocities
    ydot[2] = y[0];
    ydot[3] = y[1];

    return ydot;
  }  


// Function to compute potential energy per mass in R (nsel_physics=1) or in R'(nsel_physics=2)
double get_Epot(double xx, double yy) {
    double V=0;
    if (nsel_physics==1) {
      double r = sqrt(pow(xx,2)+pow(yy,2));

      V= -GM * mS / r; 
      return V;
    }

    else if (nsel_physics==2){

      double x_prime_S = -alpha*a;
      double x_prime_J = beta*a;

      double rS = sqrt( pow((xx - x_prime_S ),2) + pow( yy,2) ); 
      double rJ = sqrt( pow((xx- x_prime_J ),2) + pow( yy,2) );  

      V = - GM * (mS/rS + mJ/rJ) - 0.5*pow(Omega,2)*(pow(xx,2)+pow(yy,2));
      return V;
    }
    else{
      cerr << "No dynamics corresponds to this index" << endl;
      return 0;
    }
}

// Function to compute mechanical energy per mass in R'
double compute_energy(double vx, double vy, double xx, double yy) { //PER MASS
    double T = 0.5*(pow(vx,2)+pow(vy,2));
    if (nsel_physics==1 or nsel_physics==2) {
      return T + get_Epot(xx,yy);
    }
    else{
      cerr << "No dynamics corresponds to this index" << endl;
      return 0;
    }
}
void initial_condition(void){
  if(nsel_physics==1){
    x0 = 2*a;
    y[0] = v0[0]; //vx0
    y[1] = v0[1]; //vy0
    y[2] = x0; //x0
    y[3] = 0;
  }
  else{
    double xS = alpha*a; 
    x0 = 2*a + xS; 
    y[0] = v0[0]; //vx'0
    y[1] = v0[1] - Omega*x0; //vy'0
    y[2] = x0; //x'0
    y[3] = 0; //y'0
  }
}

std::valarray<double> RK4_do_onestep(const std::valarray<double>& yold, double ti, double dt) {
    std::valarray<double> k1, k2, k3, k4, ynew;
    //TO DO
      
    k1=dt*get_f(ti , yold);
    k2=dt*get_f( ti+(dt/2.0) , yold + (k1/2.0));
    k3=dt*get_f( ti+dt/2.0 , yold + (k2/2.0) );
    k4=dt*get_f( ti+dt , yold+k3 );
    
    ynew= yold + ( k1+2.0*k2+2.0*k3+k4 )/6.0;  //Final computation

    return ynew;
}


public:
    Exercice3(int argc, char* argv[])
    {
        const double pi = 3.1415926535897932384626433832795028841971e0;
        string inputPath("configuration.in"); // Default input file

        if (argc > 1)
            inputPath = argv[1];

        ConfigFile configFile(inputPath);
        for (int i = 2; i < argc; ++i){
          std::cout << "[CLI ARG] " << argv[i] << std::endl;
          configFile.process(argv[i]);}

        tFin         = configFile.get<double>("tFin");
        mS         = configFile.get<double>("m1"); //sun mass
        mJ         = configFile.get<double>("m2"); //jupiter mass
        v0[0]          = configFile.get<double>("v0x");
        v0[1]         = configFile.get<double>("v0y");
        a            = configFile.get<double>("a");
        nsel_physics = configFile.get<int>("nsel_physics");
        adapt        = configFile.get<bool>("adapt");
        tol          = configFile.get<double>("tol");
        sampling     = configFile.get<int>("sampling");
        nsteps       = configFile.get<int>("nsteps"); 
        mtot = mS + mJ; // Total mass of the system
        alpha = mJ / mtot;
        beta  = mS / mtot;
        Omega = sqrt(GM*mtot/pow(a,3));
        dt = tFin / nsteps;
        outputFile = new ofstream(configFile.get<string>("output").c_str());
        outputFile->precision(15);
        
    }

    ~Exercice3()
    {
        outputFile->close();
        delete outputFile;
    }

    void run()
{
    t = 0.0;
    initial_condition();
    last = 0;
    printOut(true); // First output

    if (!adapt) {
        cout << "Non-adaptive case" << endl;

        while (t < tFin - 0.5 * dt) {
            y = RK4_do_onestep(y, t, dt);
            t += dt;
            printOut(false); // Regular sampling
        }
    } 
    else {
        cout << "Adaptive case" << endl;
        nsteps = 0;
        const int n = 4;
        const double f = 0.9;

        while (t < tFin) {
            dt = min(dt, tFin - t);

            valarray<double> y_full  = RK4_do_onestep(y, t, dt);
            valarray<double> y_half1 = RK4_do_onestep(y, t, dt / 2.0);
            valarray<double> y_half2 = RK4_do_onestep(y_half1, t + dt / 2.0, dt / 2.0);
            double d = abs(y_full - y_half2).sum();
            cout << "[ACCEPTED] t = " << t << ", dt = " << dt << ", d = " << d << ", tol = " << tol << endl;
            int refinements = 0;
            if (d <= tol) {
                y = y_half2;
                t += dt;
                dt = f * dt * pow(tol / d, 1.0 / (n + 1));
            } 
            else {
                while (d > tol) {
                    ++refinements;
                    dt = f * dt * pow(tol / d, 1.0 / (n + 1));
                    y_full  = RK4_do_onestep(y, t, dt);
                    y_half1 = RK4_do_onestep(y, t, dt / 2.0);
                    y_half2 = RK4_do_onestep(y_half1, t + dt / 2.0, dt / 2.0);
                    d = abs(y_full - y_half2).sum();
                    cout << "[REFINING] t = " << t << ", dt = " << dt << ", d = " << d << ", tol = " << tol << endl;
                }
                y = y_half2;
                t += dt;
            }
            ++nsteps;
            printOut(false); // Regular sampling
        }
        cout << "Total accepted steps: " << nsteps << endl;
        cout << "Final time: " << t << endl;
    }

    printOut(true); // Final output
}
};

int main(int argc, char* argv[])
{
  Exercice3 engine(argc, argv);
  engine.run();

  return 0;
}
