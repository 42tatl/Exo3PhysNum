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
  double r0, r1, v0, a;
  double GM=6.674e-11;
  double mtot=0e0;
  double L2x,L2y;
  valarray<double> m = std::valarray<double>(0.e0, 2);
  int N_excit, nsteps;
  int sampling;
  int last;
  int  nsel_physics;
  bool adapt;
  double alpha = 0e0;
  double beta = 0e0;
  double tol= 0e0;
  valarray<double> x0 = std::valarray<double>(0.e0, 4); // Correctly initialized
  valarray<double> x  = std::valarray<double>(0.e0, 4); // Correctly initialized
  ofstream *outputFile;

  void printOut(bool write)
  {
    if((!write && last>=sampling) || (write && last!=1))
    {
      double Energy = compute_energy(x[0],x[1],x[2],x[3]);
      *outputFile << t << " "<< x[0] << " " << x[1] << " "<< x[2] << " " << x[3] << " " \
      << Energy<< " "<< nsteps<< endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }

  std::valarray<double> get_f(double t, const std::valarray<double>& x) {
    std::valarray<double> xdot(0.0, 4);
    //TO DO
    if (nsel_physics == 1) { //two body problem 
      
      double r = sqrt(pow(x[0],2)+ pow(x[1],2));
      xdot[2] = -(GM*m[1]*x[0])/pow(r,3);
      xdot[3] = -(GM*m[1]*x[1])/pow(r,3);
    }
    else if (nsel_physics == 2) { //reduced three body problem
      double rs = sqrt( pow((x[0]+alpha*a ),2)+pow( x[1],2) );
      double rt = sqrt( pow((x[0]-beta*a ),2)+pow( x[1],2) );
      
      double Omega = sqrt(GM*(mtot)/pow(a,3));
      
      xdot[2] = -GM*m[0]/pow(rs, 3)*(x[0] + alpha*a) -GM*m[1]/pow(rt, 3)*(x[0] - beta*a) +pow(Omega,2)*x[0] + 2*Omega*x[3];
      xdot[3] =  -GM*m[0]/pow(rs, 3)*x[1] -GM*m[1]/pow(rt, 3)*x[1] +pow(Omega,2)*x[1] - 2*Omega*x[2];
    }
    else{
        cerr << "No dynamics corresponds to this index" << endl;
        return xdot;
    }

    // Velocities
    xdot[0] = 0;
    xdot[1] = 0;

    return xdot;
  }  


// Function to compute potential energy per mass in R (nsel_physics=1) or in R'(nsel_physics=2)
double get_Epot(double xx, double yy) {
    double V=0;
    if (nsel_physics==1) {
      V=-GM*m[0]/sqrt(pow(xx,2)+pow(yy,2))-GM*m[1]/sqrt(pow(xx,2)+pow(yy,2));
      return V;
    }
    else if (nsel_physics==2){
      double rs = sqrt( pow((xx+alpha*a ),2)+pow( yy,2) );
      double rt = sqrt( pow((xx-beta*a ),2)+pow( yy,2) );  
      double Omega = sqrt(GM*(mtot)/pow(a,3));
      V = -GM*m[0]/rs - GM*m[1]/rt - 0.5*pow(Omega,2)*(pow(xx,2)+pow(yy,2));
      return V;
    }
    else{
      cerr << "No dynamics corresponds to this index" << endl;
      return 0;
    }
}

// Function to compute mechanical energy per mass in R'
double compute_energy(double xx, double yy, double vx, double vy) {
    double T = 0.5*(pow(vx,2)+pow(vy,2));
    if (nsel_physics==1){
      return get_Epot(xx,yy) + T;
    }
    else if (nsel_physics==2){
      return get_Epot(xx,yy) + T;
    }
    else{
      cerr << "No dynamics corresponds to this index" << endl;
      return 0;
    }

    return 0;
}
void initial_condition(void){
  if(nsel_physics==1){
    //TO DO initialize x0
  }
  else{
    //TO DO initialize x0
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
    const double pi=3.1415926535897932384626433832795028841971e0;
    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice3 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    tFin         = configFile.get<double>("tFin");            // t final (overwritten if N_excit >0)
  //  dt       = configFile.get<double>("dt");            // time step (overwritten if nsteps_per >0)
    m[0]         = configFile.get<double>("m1");              // mass of the sun
    m[1]         = configFile.get<double>("m2");              // mass of the earth
    r0           = configFile.get<double>("r0");              // r0
    r1           = configFile.get<double>("r1");              // r1
    L2x          = configFile.get<double>("L2x");              // L2x
    L2y          = configFile.get<double>("L2y");              // L2y
    a            = configFile.get<double>("a");               // demi grand-axe (solei-terre hyp MCU)
    nsel_physics = configFile.get<int>("nsel_physics");       //1) one body problem around mass#2 or 2) one body in rotating reference frame of {1,2}
    adapt        = configFile.get<bool>("adapt");             //if 1=true -> adaptive dt, if 0=false -> fixed dt
    tol          = configFile.get<double>("tol");             //tolerance of the adaptive scheme
    sampling     = configFile.get<int>("sampling");     // number of time steps between two writings on file
    nsteps       = configFile.get<int>("nsteps");        // number of time step per period
    mtot=m[0]+m[1];
    alpha = m[1] / mtot;
    beta = m[0] / mtot;
    //TO DO
    
    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
    //TO DO initialize tFin for nsel_physics=1 and initialize dt for both nsel_physics
    dt=tFin/nsteps;
  }

  ~Exercice3()
  {
    outputFile->close();
    delete outputFile;
  };

  void run()
  {
    t = 0.;
    initial_condition();
    x=x0;
    last = 0;
    printOut(true);
    std::valarray<double> y1;
    std::valarray<double> y2;
    if (adapt==false){
      cout<<"non adaptive"<<endl;
      while(t<tFin-0.5*dt){
        x = RK4_do_onestep(x,t,dt);
        t += dt;
        printOut(false); 
      } 
        
    }
    else{
      //TO DO adaptive case
      cout<<"Adaptive case"<<endl;
      nsteps = 0;
    };
    
    printOut(true); // ecrire le dernier pas de temps
  };

};

int main(int argc, char* argv[])
{
  Exercice3 engine(argc, argv);
  engine.run();

  return 0;
}
