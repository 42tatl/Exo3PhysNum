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
  double Omega = sqrt(GM*(mtot)/pow(a,3));
  valarray<double> m = std::valarray<double>(0.e0, 2);
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
      
      double r = sqrt(pow(y[2],2)+ pow(y[3],2)); //y = (vx, vy, x, y) PROBLEME ICI: PQ ON PREND PAS X_S
    
      ydot[0] = -(GM*m[0]*y[2])/pow(r,3); 
      ydot[1] = -(GM*m[0]*y[3])/pow(r,3);
    }
    else if (nsel_physics == 2) { //reduced three body problem
      double x_prime_S = -alpha*a; 
      double x_prime_J = beta*a;

      double rS = sqrt( pow(y[2] - x_prime_S,2)+pow( y[3],2) );
      double rJ = sqrt( pow(y[2] - x_prime_J,2)+pow( y[3],2) ); //y = (vx, vy, x, y)
            
      ydot[0] = -GM * (m[0]/pow(rS, 3)*(y[2] - x_prime_S) - m[1]/pow(rJ, 3)*(y[2] - x_prime_J)) +pow(Omega,2)*y[2] + 2*Omega*y[1];
      ydot[1] =  -GM * y[3] * (m[0]/pow(rJ, 3) - m[1]/pow(rJ, 3)) + pow(Omega,2)*y[3] - 2*Omega*y[0];
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

      V= -GM * m[0] / r; //PB: POURQUOI ON AURAIT M[0]+M[1] DANS LE CAS 1
      return V;
    }

    else if (nsel_physics==2){

      double x_prime_S = -alpha*a;
      double x_prime_J = beta*a;

      double rS = sqrt( pow((xx - x_prime_S ),2) + pow( yy,2) ); 
      double rJ = sqrt( pow((xx- x_prime_J ),2) + pow( yy,2) );  

      V = - GM * (m[0]/rS + m[1]/rJ) - 0.5*pow(Omega,2)*(pow(xx,2)+pow(yy,2));
      return V;
    }
    else{
      cerr << "No dynamics corresponds to this index" << endl;
      return 0;
    }
}

// Function to compute mechanical energy per mass in R'
double compute_energy(double xx, double yy, double vx, double vy) { //PER MASS
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
    y[0] = v0[0];
    y[1] = v0[1];
    y[2] = x0;
    y[3] = 0;
  }
  else{
    y[0] = v0[0];
    y[1] = v0[1] - Omega*x0;
    y[2] = x0 - alpha*a;
    y[3] = 0;
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
        for (int i = 2; i < argc; ++i)
            configFile.process(argv[i]);

        tFin         = configFile.get<double>("tFin");
        m[0]         = configFile.get<double>("m1"); //sun mass
        m[1]         = configFile.get<double>("m2"); //jupiter mass
        v0[0]          = configFile.get<double>("v0x");
        v0[1]         = configFile.get<double>("v0y");
        a            = configFile.get<double>("a");
        nsel_physics = configFile.get<int>("nsel_physics");
        adapt        = configFile.get<bool>("adapt");
        tol          = configFile.get<double>("tol");
        sampling     = configFile.get<int>("sampling");
        nsteps       = configFile.get<int>("nsteps");

        mtot = m[0] + m[1];
        alpha = m[1] / mtot;
        beta  = m[0] / mtot;

        outputFile = new ofstream(configFile.get<string>("output").c_str());
        outputFile->precision(15);

        dt = tFin / nsteps;
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
        printOut(true);

        if (!adapt) {
          cout << "non adaptive" << std::endl;
          while (t < tFin - 0.5 * dt) {
              y = RK4_do_onestep(y, t, dt);
              t += dt;
              printOut(false);
          }
        } 
        else {
          cout << "Adaptive case" << std::endl;
          nsteps = 0;
          double dt_c;
          const int n = 4;
          const double f = 0.9;

          while (t < tFin - 0.5 * dt) {
            dt_c = dt;
            valarray<double> y_full  = RK4_do_onestep(y, t, dt_c);
            valarray<double> y_half1 = RK4_do_onestep(y, t, dt_c / 2.0);
            valarray<double> y_half2 = RK4_do_onestep(y_half1, t + dt_c / 2.0, dt_c / 2.0);
            double d = abs(y_full - y_half2).sum();
            while (d > tol) {
                dt_c = f * dt * pow(tol / d, 1.0 / (n + 1));
                y_full  = RK4_do_onestep(y, t, dt_c);
                y_half1 = RK4_do_onestep(y, t, dt_c / 2.0);
                y_half2 = RK4_do_onestep(y_half1, t + dt_c / 2.0, dt_c / 2.0);
                d = abs(y_full - y_half2).sum();
            }

            y = y_half2;
            t += dt_c;
            dt = dt_c * pow(tol / d, 1.0 / (n + 1)); //je suis pas sure de comprendre si la boucle recommence si on a direct d<epsilon
            dt = std::min(dt, tFin - t);

            ++nsteps;
            printOut(false);
          }
      }
      printOut(true); // final output
    }
};

int main(int argc, char* argv[])
{
  Exercice3 engine(argc, argv);
  engine.run();

  return 0;
}
