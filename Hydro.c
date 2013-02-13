#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>


#define gamma 1.4
#define n 3
#define safety 0.9

//Returns max value in y
double max(const double y[], int ng){
    int i=0;
    double max=y[0];
    for (i=0; i<ng; i++){
        if (y[i]>max){
            max=y[i];
        }
    }
    return max;
}

//Our equation of state
double pres(const double u[]){
    double rho=u[0];
    double v=u[1]/u[0];
    double E=u[2];


    double e=(E-(0.5*rho*pow(v,2)))/rho;
    return (gamma-1)*rho*e;
}

void flux(const double u[], double F[], int ng){
    //Loop index
    int i=0;
    //Velocity density pressure and energy density
    double v=0;
    double rho=0;
    double p=0;
    double E=0;

    //For each grid point compute fluxes for fluid equation
    for (i=0; i<ng*n; i+=n){
        v=u[i+1]/u[i];
        rho=u[i];
        p=pres(&(u[i]));
        E=u[i+2];



        F[i]=rho*v;
        F[i+1]=(rho*pow(v,2))+p;
        F[i+2]=(E+p)*(v);



    }
}


double lplus(const double u[]){
    double rho=u[0];
    double v=u[1]/u[0];
    double p=pres(u);
    double cs=pow(gamma*p/rho, 0.5);

    return v+cs;
}

double lminus(const double u[]){
    double rho=u[0];
    double v=u[1]/u[0];
    double p=pres(u);
    double cs=pow(gamma*p/rho, 0.5);

    return v-cs;
}

void alpha_plus(const double u[],  double aplus[], int ng){
    //Loop variable
    int i=0;
    double l1, l2;


    for(i=0; i<ng-1; i++){
        l1=lplus(&(u[i]));
        l2=lplus(&(u[i+n]));

        aplus[i]=(l1>l2)?l1:l2;
        aplus[i]=(aplus[i]>0)?aplus[i]:0;
    }
}

void alpha_minus(const double u[], double aminus[], int ng){
    //Loop variable
    int i=0;
    double l1, l2;

    for(i=0; i<ng-1; i++){
        l1=lplus(&(u[i]));
        l2=lplus(&(u[i+n]));

        aminus[i]=(-l1>-l2)?l1:l2;
        aminus[i]=(aminus[i]>0)?aminus[i]:0;
    }
}


//Reimann solver
double hll(const double u[], const double F[], const double aplus[], const double aminus[], int i, int j, int ng){

    //Flux is zero at the boundaries
    if(i<0){
        return 0;
    }
    if (i>=(ng-n)){
        return 0;
    }


    return (aplus[i]*F[i+j]+aminus[i]*F[i+j+n]-(aplus[i]*aminus[i])*(u[i+j+n]-u[i+j]))/(aplus[i]+aminus[i]);
}


void derivs(double t, const double u[], double L[], int ng, double delta_x){
  //Loop variables
  int i=0;
  int j=0;


  double F[ng*n];
  //Calculating the flux at each of the grid points
  flux(u, F, ng);


  double aplus[ng-1], aminus[ng-1];
  alpha_plus(u,  aplus, ng);
  alpha_minus(u, aminus, ng);



  //For each of the grid points
  for (i=0; i<ng; i+=n){
      for (j=0; j<n; j++){
          double Fi_down=hll(u, F, aplus, aminus, i-n, j, ng);
          double Fi_up=hll(u, F, aplus, aminus, i, j, ng);
          //Time derivative
          L[i+j]=-(Fi_up-Fi_down)/(delta_x);
      }
  }



}

/*Forward euler scheme: given current solution vector and time it computes vector of derivatives and then
evolves the solution vector forward by one time step*/
void Euler(double t, double u[], double L[], double delta_t,  int ng, double delta_x, void(* myderivs)(double,  const double *, double *, int,  double)){
    //Loop variable
    int i=0;

    //Calculate derivatives
    myderivs(t, u, L, ng, delta_x);

    for (i=0; i<ng*n; i++){
        u[i]+=delta_t*L[i];
    }
}

//Return the Courant-Friedrick-Levy time step multiplied by some safety factor
double cfl(const double aplus[], const double aminus[], int ng, double delta_x){
    //Maximum of alpha pluses and minuses
    double m1=max(aplus, ng-1);
    double m2=max(aminus, ng-1);
    //Maximum of the alpha pluses and minuses
    double m=(m1>m2)?m1:m2;


    return safety*(delta_x/m);
}


int main(int argc, char* argv[]){

    int ng=2;
    double u[]={1,1,1, 0.5, 0, 0};
    double F[ng*n];
    flux(u, F, ng);

    double aplus[ng-1], aminus[ng-1];

    alpha_plus(u, aplus, ng);
    alpha_minus(u, aminus, ng);
    printf("%lf %lf\n", aminus[i], aplus[i]);


//    FILE* out=fopen("hydro.txt", "w");
//
//    //Loop variable
//    int i=0;
//    //Number of grid points
//    int ng=100;
//    //x boundaries of our numerical domain
//    double xmin=0;
//    double xmax=1;
//    //size of each of our cells
//    double delta_x=(xmax-xmin)/ng;
//
//
//
//    //Initial condition
//    double u[ng*n];
//    double rho_l=1;
//    double v_l=0;
//    double p_l=1;
//    double rho_r=0.1;
//    double v_r=0;
//    double p_r=0.125;
//
//    printf("%d\n", (ng*n)/2);
//
//    for(i=0; i<((ng*n)/2); i+=n){
//        u[i]=rho_l;
//        u[i+1]=rho_l*v_l;
//        u[i+2]=p_l/(gamma-1)+(rho_l*pow(v_l,2));
//
//    }
//    printf("%d\n", i);
//    for (i=((ng*n)/2); i<(ng*n); i+=n){
//        u[i]=rho_r;
//        u[i+1]=rho_r*v_r;
//        u[i+2]=p_r/(gamma-1)+(rho_r*pow(v_r,2));
//    }
//
//    //Alpha's
//    double aplus[ng-1], aminus[ng-1];
//    //Time derivatives
//    double L[ng*n];
//    //Time step we are going to take
//    double delta_t=0;
//    //Keeping track of time that went by
//    double t=0;
//    //Keep track of number of iterations
//    int iterations=0;
//    while((t<1)&&(iterations<1)){
//        //Calculate alpha at each grid point
//        alpha_minus(u, aminus, ng);
//        alpha_plus(u, aplus, ng);
//        //Calculate time step
//        delta_t=cfl(aplus, aminus, ng, delta_x);
//
//
//        for (i=0; i<ng-1; i++){
//            printf("%lf %lf\n", aminus[i], aplus[i]);
//        }
//
//        for (i=0; i<ng*n; i+=n){
//            fprintf(out, "%lf %lf %lf\n", t, delta_x*i, u[i]);
//        }
//
//        //Evolve forward one time step
//        Euler(t, u, L, delta_t, ng, delta_x, &derivs);
//        t+=delta_t;
//        iterations++;
//
//    }

    return 0;
}


