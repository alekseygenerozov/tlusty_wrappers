#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

//Pi
#define PI 3.141592
//Mass of sun in grams
#define Msun 2E33
//Black hole mass in solar
#define  M 1E7*Msun
//Mean molecular weight
#define  mu 0.615
//Mass of proton
#define  mp 1.67E-24
//Boltzmann constant
#define kb 1.38E-16
//Stefan-Boltzmann constant in cgs
#define  sigma 5.67E-5
//Opacity: electron scattering 
#define kappa 0.4
//Newton's constant
#define  G 6.67E-8
//Define speed of light constant
#define c 3E10
//Radius we are considering
//#define  R 400*(2*G*M/pow(c, 2))





//Equation of state: ideal gas plus radiation Pressure
double Pres(double rho, double T);
double PresGas(double rho, double T);
double PresRad( double T);
double alpha(double rho, double T);
double mygamma(double rho, double T);

//Analytic expression for the flux
double F(double u, double nu, double R);

//dlgT/dlgP:
double delrad(double u, const double y[], void* params);
double delad(double u, const double y[], void* params);
double del(double u, const double y[], void* params);

//Acceleration of gravity
double g(double z, double R);

int func(double u, const double y[], double f[], void* params){

  //The dependent variables for the problem
  double z=y[0];
  double rho=y[1];
  double T=y[2];

  double nu=((double *)params)[0];
  double R=((double *)params)[1];



  //Return right hand zide of our ode
  f[0]=(1.0/rho);
  f[2]=(-T/Pres(rho, T))*g(z, R)*del(u, y, params);
  f[1]=(-rho/T)*f[2]-(16*sigma/(3*c))*pow(T, 2)/(kb/(mu*mp))*f[2]-(g(z, R)/(kb*T/(mu*mp)));

  return GSL_SUCCESS;
}

int jac(double u, const double y[], double *dfdy, double dfdt[], void *params){
    return GSL_SUCCESS;
}



//using namespace std;
int main(int argc, char* argv[]){
    //Check to see that number of arguments is appropriate
    if(argc!=4){
      printf("Useage:./a.out Sigma Mdot Radius\n");
      return 1;
    }
    //Getting the surface density Mdot and radius radius from command line args.
    double Sigma0=atof(argv[1]);
    //Second argument is interpreted to be mass accretion rate as fraction of Eddington
    double mdot=atof(argv[2]);
    double R=atof(argv[3])*(2*G*M/pow(c, 2));
    double Mdot=mdot*10*4*PI*G*M/(kappa*c);

    printf("%e\n", Mdot);

    //Loop variable
    int i=0;


    //Get viscosity and total optical depth
    double nu=Mdot/(3*PI*Sigma0);
    double TauTot=0.5*Sigma0*kappa;
    //Get effective temperature and the central temperature 
    double Teff=pow(((9.0/8.0)*nu*Sigma0*G*M/pow(R, 3))/sigma, 0.25);
    double Tc=pow((3.0/8.0)*TauTot, 0.25)*Teff;
    


    //Indicates attempted number of density bisections
    int tries=0;
    //Max number of tries
    int MaxTries=1E4;
    //Steps in the RK4 integrator
    int steps=0;
    //Max number of steps allowed
    int MaxSteps=1E6;

    /*number of dependent variables for our system*/
    int n=3;

    //Brackets for iterative bisection of density
    double rhoc_min_pow=-10;
    double rhoc_max_pow=-1;
    //Central density
    double rhoc_pow= (rhoc_max_pow+rhoc_min_pow)/2.0;
    double rhoc=pow(10, rhoc_pow);
    //Store the central density of the previous step
    double rhoc_prev=rhoc;


    //boundary conditions and stepsize
    double y[]={1, rhoc, Tc};
    double u=rhoc;
    double du=500;
    //Parameters for problem
    double params[]={nu, R};
    //Relative and absolute error tolerances
    double eps_abs=0;
    double eps_rel=1E-6;


    gsl_odeiv2_system sys={func, jac, n, params};
    //gsl_odeiv2_driver* d=gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkck, du, 0, 1E20);
    gsl_odeiv2_evolve* e=gsl_odeiv2_evolve_alloc (n);
    gsl_odeiv2_control* con = gsl_odeiv2_control_standard_new (eps_abs,  eps_rel, 0, 1);
    gsl_odeiv2_step * stepfunc=gsl_odeiv2_step_alloc (gsl_odeiv2_step_rkck, n);


    //Initialize surface density
    double Sigma=0;
    //Tolerance for the convergence of the surface density
    double tol=1E-10;

    //Our output file
    char* bp, tmp;
    size_t size;
    FILE* outfile=open_memstream(&bp, &size);
    FILE* out;
    //Keeping track of the relative errors
    FILE* errors;
    FILE* bisect=fopen("bisect", "w");

 
    //Creating file name based on user input
    fprintf(outfile, "profile"); 
    for(i=1; i<argc; i++){
      //printf("test\n");
      fprintf(outfile, "-%s", argv[i]);
     }
    fflush(outfile);
    fclose(outfile);

  
    printf ("%s\n", bp);
   
    //While we have not converged to the nominal surface density
    while(fabs(Sigma-Sigma0)/Sigma0>tol){
      fprintf(bisect, "%d %lf %.16e\n", tries,  Sigma-Sigma0, rhoc);

      //Reinitialize the solution array
       rhoc_pow=(rhoc_max_pow+rhoc_min_pow)/2.0;
       rhoc_prev=rhoc;
       rhoc=pow(10, rhoc_pow);

       
       //reinitialize column density and the number of steps.
        Sigma=0;
        steps=0;
	

    	
        u=rhoc;
        y[0]=1;
        y[1]=rhoc;
        y[2]=Tc;


        //Increment number of bisections
        tries=tries+1;
	

        if (tries>MaxTries){
            printf("Exceeded max number of tries (%d) for the bisection phase\n", MaxTries);
            break;
        }
        


        //Our output file
        out=fopen(bp, "w");
        errors=fopen("errors", "w");

    	//Inifinite loop with break conditions for integrating our odes
        while (1==1){
            fprintf(out, "%.4e ", u);
            for (i=0; i<n; i++){
    	      //out<<x[i]<<" ";
              fprintf(out, "%.4e ", y[i]);
            }
            fprintf(out, "%.4e ", du);
	    fprintf(out, "%.4e %.4e\n", delad(u, y, params), delrad(u, y, params));


	    


            //printf("%lf\n", u);
            //Go forward one step
            //int status=gsl_odeiv2_driver_apply(d, &u, u+du, y);
            int status=gsl_odeiv2_evolve_apply(e, con, stepfunc, &sys, &u, u+du, &du,  y);

            for (i=0; i<n; i++){
              fprintf(errors, "%.4e %.4e ", ((e->yerr)[i]), eps_rel*(y[i]));
             }
            fprintf(errors, "\n");
            //printf("%lf\n", u);


            if(status != GSL_SUCCESS){
              //printf("test\n");
	      break;
            }

            if (!((y[0]>0) && (y[1]>0) && (y[2]>0))){
                break;
            }

            /*Add to total surface density; Note total surface density is defined as the integral from of the density from negative infinity to infinity*/
            Sigma=2*u;
            if (Sigma>(1+tol)*Sigma0){
              //printf("%lf %lf\n", Sigma, u);
              break;
            }
            
            

    	    //Increment number of steps; if we have exceeded the maximum number of steps allowed then break.
            steps+=1;
            //If we have exceeded our maximum allowed # of steps then break
            // if (steps > MaxSteps)
            //     break;


        }

        fclose(out);
        fclose(errors);
        

        if (Sigma<Sigma0){
    	  rhoc_min_pow=rhoc_pow;
    	}
        else{
    	  //cout<< rhoc_max_pow<<endl;
    	  rhoc_max_pow=rhoc_pow;
    	}
        //printf("%lf\n", rhoc);
    	//cout<<Sigma<<" "<<Sigma0<<" "<<rhoc_min_pow<<" "<<rhoc_max_pow<<" "<<setprecision(10)<<rhoc<<endl;
    }
    printf("%d\n",tries);
    return 0;
}


//our equation of state
double Pres(double rho, double T){
    return (rho*kb*T)/(mu*mp)+(4.0/3.0)*(sigma*pow(T, 4)/c);
}
//Gas pressure and radiation pressure
double PresRad(double T){
    return (4.0/3.0)*(sigma*pow(T, 4)/c);
}
double PresGas(double rho, double T){
    return (rho*kb*T)/(mu*mp);
}
double alpha(double rho, double T){
    double beta=PresRad(T)/PresGas(rho, T);
    return ((1/T)+(4*beta*(1/T)));
}
//Ratio of specific heats
double mygamma(double rho, double T){
     double a=alpha(rho, T);

     return 1 + (  (pow(a, 2)*T)/(rho*(1/PresGas(rho, T))* (3.0/2.0)*  (kb/(mu*mp))  ));
}
//gravitational acceleration
double g(double z, double R){
    return G*M*z/pow(pow(R, 2)+pow(z, 2), 1.5);
}

//the flux (analytically calculated)
double F(double u, double nu, double R){
    return (9.0/4.0)*nu*u*(G*M/pow(R,3));
}


//Adiabatic and radiative temperature gradients.
double delrad(double u, const double y[], void* params){
    double z=y[0];
    double rho=y[1];
    double T=y[2];

    double nu=((double *)params)[0];
    double R=((double *)params)[1];


    return F(u, nu, R)/(16*sigma*pow(T,4)/(  3*kappa*Pres(rho, T)  ) *g(z, R));
}
double delad(double u, const double y[], void* params){
    double z=y[0];
    double rho=y[1];
    double T=y[2];




    //gas pressure and radiation pressure
    double Prad=PresRad(T);
    double Pgas=PresGas(rho, T);

    //quantities that are useful for defining the adiabatic gradient
    double beta=Prad/Pgas;
    double a=alpha(rho, T);
    double g=mygamma(rho, T);

    //Adiabatic temperature gradient
    return ((1 + beta)*((4*beta) + 1))/(((12.*beta) + (1/(g - 1))) +
       pow((4*beta) + 1, 2));
}

//Radiative del
double del(double u, const double y[], void* params){
    //Radiative temperature gradient
    double del_rad= delrad(u, y, params);
    //Adiabatic temperature gradient
    double del_ad = delad(u, y, params);

    //check for convective instability
    if (del_rad>del_ad){
        return del_ad;
    }
    return del_rad;

}



/* //RK4 implementation */
/* void RK4(double* x, int n, double h){ */
/*   //Loop variable */
/*   int i=0; */

/*   //Intermediate slopes and solutions */
/*   double k1[n], k2[n], k3[n], k4[n]; */
/*   double x2[n], x3[n], x4[n]; */

/*   //calculate derivative at starting point */
/*   derivs(x, k1); */
/*   //take half-step using k1 */
/*   for (i=0; i<n; i++){ */
/*     x2[i]=x[i]+(k1[i]*0.5*h); */
/*   } */

/*   //calculate derivative at x2 */
/*   derivs(x2, k2); */
/*   //take a half step using k2 */
/*   for (i=0; i<n; i++){ */
/*     x3[i]=x[i]+(k2[i]*0.5*h); */
/*   } */

/*   //calculate derivative at x3 */
/*   derivs(x3, k3); */
/*   //take a full step using k3 */
/*   for (i=0; i<n; i++){ */
/*     x4[i]=x[i]+(k3[i]*0.5*h); */
/*   } */

/*   //estimate the derivative at x4 */
/*   derivs(x4, k4); */

/*   //Now take a full step */
/*   for (i=0; i<n; i++){ */
/*     x[i]+=(1.0/6.0)*(k1[i]+(2*k2[i])+(2*k3[i])+k4[i])*h; */
/*   } */
/* } */

/* //Our system of ODEs */
/* void derivs(double* x, double* k){ */
/*        double u=x[0]; */
/*        double z=x[1]; */
/*        double rho=x[2]; */
/*        double T=x[3]; */
/*        double nu=x[4]; */
/*        double R=x[5]; */

/*        k[0]=1; */
/*        k[1]=(1.0/rho); */
/*        k[3]=(-T/Pres(rho, T))*g(z, R)*del(x); */
/*        k[2]=(-rho/T)*k[3]-(16*sigma/(3*c))*pow(T, 2)/(kb/(mu*mp))*k[3]-(g(z, R)/(kb*T/(mu*mp))); */
/*        k[4]=0; */
/*        k[5]=0; */
/* } */
