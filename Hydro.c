#include <math.h>
#include <stdio.h>
#include <stdlib.h>



#define gamma 1.4
#define n 3
#define safety 0.9
#define theta 1.5


struct cell{
  //Values of flux, conservative variables, and primitive variables for a cell
  double* F;
  double* u;
  double* v;

  double u0[n];

  //Size of cell
  double delta_x;

  //left and right biased values for each interface
  double F_ll[n], F_lr[n], F_rl[n], F_rr[n];
  double u_ll[n], u_lr[n], u_rl[n], u_rr[n];
  double v_ll[n], v_lr[n], v_rl[n], v_rr[n];
  //fluxes at interfaces
  double Fi_l[n], Fi_r[n];
  //Time derivative estimated from the fluxes at interfaces.
  double L[n];

  

  double aplus_r, aplus_l, aminus_r, aminus_l;



};

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

//Returns max value in y
double min(const double y[], int ng){
    int i=0;
    double min=y[0];
    for (i=0; i<ng; i++){
        if (y[i]<min){
            min=y[i];
        }
    }
    return min;
}
//Returns the sign of x
double sgn(double x){
    if (x<0){
        return -1;
    }
//    else if(x==0){
//        return 0;
//    }
    return 1;
}
//Minmod function
double minmod(double x, double y, double z){
    double tmp[3]={fabs(x), fabs(y), fabs(z)};


    return (1./4.)*fabs(sgn(x)+sgn(y))*(sgn(x)+sgn(z))*min(tmp, 3);
}



//void plm(const double c[], double c1, char x){
//    if (x=='l'){
//            c1=c[1]+minmod(theta*(c[1]-c[0]), 0.5*(c[2]-c[0]), theta*(c[2]-c[1]));
//    }
//
//}

//Computes specific energy from primitive variables given our equation of state
double en(const double v[]){
  double rho=v[0];
  double vel=v[1];
  double p=v[2];

  return p/(rho*(gamma-1));

}

//Converts conservative variables to primitive variables
void to_prim(const double u[], double v[]){
  double rho=u[0];
  v[0]=rho;


  double vel=u[1]/rho;
  v[1]=vel;

  double E=u[2];
  double pres=(E-(1./2.)*rho*pow(vel,2))*(gamma-1);
  v[2]=pres;
 }

//Function that converts primitive variables to conservative variables
void to_cons(const double v[], double u[]){
  double rho=v[0];
  double vel=v[1];
  double p=v[2];
  //specific internal energy
  double e=en(v);



  u[0]=rho;
  u[1]=rho*vel;
  u[2]=(1./2.)*rho*pow(vel,2)+(rho*e);

}

//Calculates the flux given a vecor u at the and stores it in vector F
void to_flux(const double v[], double F[]){
    double u[n];
    to_cons(v, u);

    double rho=v[0];
    double vel=v[1];
    double p=v[2];
    double E=u[2];

    F[0]=rho*vel;
    F[1]=(rho*pow(vel,2))+p;
    F[2]=(E+p)*(vel);
}


double lplus(const double v[]){
    double rho=v[0];
    double vel=v[1];
    double p=v[2];
    double cs=pow(gamma*p/rho, 0.5);

    return vel+cs;
}

double lminus(const double v[]){
    double rho=v[0];
    double vel=v[1];
    double p=v[2];
    double cs=pow(gamma*p/rho, 0.5);

    return vel-cs;
}

double alpha_plus(const double vl[], const double vr[]){
    //Lambda's
    double l1, l2;
    //Calculating the values of lambdas
    l1=lplus(vr);
    l2=lplus(vl);

    double tmp[]={l1, l2, 0};
    return max(tmp, n);
}

double alpha_minus(const double vl[], const double vr[]){
    //Lambda's
    double l1, l2;
    //Calculating the values of lambdas
    l1=lminus(vr);
    l2=lminus(vl);

    //printf("%lf %lf test\n",l1, l2);

    double tmp[]={-l1, -l2, 0};
    return max(tmp, n);
}


//Reimann solver
void hll(double F_i[], const double u_l[], const double u_r[], const double F_l[], const double F_r[], const double aplus, const double aminus){
    int i=0;
    for(i=0; i<n; i++){
      
      F_i[i]=(aplus*F_l[i]+aminus*F_r[i]-(aplus*aminus)*(u_r[i]-u_l[i]))/(aplus+aminus);
      
    }
}
//plm functions for calculating left and right biased values
double plm_l(struct cell* grid, int i, int j){
  return grid[i].v[j]+0.5*minmod(theta*(grid[i].v[j]-grid[i-1].v[j]), 0.5*(grid[i+1].v[j]-grid[i-1].v[j]), theta*(grid[i+1].v[j]-grid[i].v[j]));
}


double plm_r(struct cell* grid, int i, int j){
  return grid[i+1].v[j]-0.5*minmod(theta*(grid[i+1].v[j]-grid[i].v[j]), 0.5*(grid[i+2].v[j]-grid[i].v[j]), theta*(grid[i+2].v[j]-grid[i+1].v[j]));
}



double cfl(struct cell* grid, int ng, double delta_x){
    //Loop variable
    int i=0;
    double aplus[ng+1], aminus[ng+1];

    for (i=0; i<ng; i++){
        aplus[i]=grid[i].aplus_l;
        aminus[i]=grid[i].aminus_l;
    }
    aplus[ng]=grid[ng-1].aplus_r;
    aminus[ng]=grid[ng-1].aminus_r;

    double m1=max(aplus, ng+1);
    double m2=max(aminus, ng+1);
    double m=(m1>m2)?m1:m2;

    return safety*(delta_x/m);

}

//Given grid and its size this function computes the fluxes at all of the interfaces
void interfaces(struct cell* grid, int ng){
  int i=0;
  int j=0;

  int tmp1, tmp2;

  for(i=0; i<ng; i++){
    for (j=0; j<n; j++){
        //Computing left and right biased values for primitive, conservative and flux variables
        grid[i].v_ll[j]=plm_l(grid, i-1, j);
        grid[i].v_rl[j]=plm_l(grid, i, j);
        grid[i].v_lr[j]=plm_r(grid, i-1, j);
        grid[i].v_rr[j]=plm_r(grid, i, j);
        to_cons(grid[i].v_ll, grid[i].u_ll);
        to_cons(grid[i].v_rl, grid[i].u_rl);
        to_cons(grid[i].v_lr, grid[i].u_lr);
        to_cons(grid[i].v_rr, grid[i].u_rr);
        //Caculating left and right biased values of fluxes
        to_flux(grid[i].v_ll, grid[i].F_ll);
        to_flux(grid[i].v_rl, grid[i].F_rl);
        to_flux(grid[i].v_lr, grid[i].F_lr);
        to_flux(grid[i].v_rr, grid[i].F_rr);



        grid[i].u_ll[j]=grid[i-1].u[j];
        grid[i].u_lr[j]=grid[i].u[j];
        grid[i].u_rl[j]=grid[i].u[j];
        grid[i].u_rr[j]=grid[i+1].u[j];

        grid[i].F_ll[j]=grid[i-1].F[j];
        grid[i].F_lr[j]=grid[i].F[j];
        grid[i].F_rl[j]=grid[i].F[j];
        grid[i].F_rr[j]=grid[i+1].F[j];

    }


    //Computing alpha's for each grid cell
    grid[i].aplus_l=alpha_plus(grid[i].v_ll, grid[i].v_lr);
    grid[i].aplus_r=alpha_plus(grid[i].v_rl, grid[i].v_rr);
    grid[i].aminus_l=alpha_minus(grid[i].v_ll, grid[i].v_lr);
    grid[i].aminus_r=alpha_minus(grid[i].v_rl, grid[i].v_rr);


    //Call hll function to compute the fluxes at interfaces
    hll(grid[i].Fi_l, grid[i].u_ll, grid[i].u_lr, grid[i].F_ll, grid[i].F_lr, grid[i].aplus_l, grid[i].aminus_l);
    hll(grid[i].Fi_r, grid[i].u_rl, grid[i].u_rr, grid[i].F_rl, grid[i].F_rr, grid[i].aplus_r, grid[i].aminus_r);
    //Calculate the time  derivative for conservative variables using the interface fluxes calculated above 
    for (j=0; j<n; j++){
      grid[i].L[j]= (-grid[i].Fi_r[j]+grid[i].Fi_l[j])/grid[i].delta_x;
    }


    //printf("%lf %lf\n", grid[i].aminus_l ,grid[i].Fi_l[0]);

  }
}

//Shu Osher implementation
void ShuOsher(struct cell* grid, int ng,  double delta_t){
    //declare loop variables
    int i,j=0;


    //Take first step and store initial values of conservative variable in u0
    interfaces(grid, ng);
    for (i=0; i<ng; i++){
        for (j=0; j<n; j++){
            grid[i].u0[j]=grid[i].u[j];
            grid[i].u[j]+=delta_t*grid[i].L[j];
        }
        //Update fluxes and primitive variables in each cell
        to_prim(grid[i].u, grid[i].v);
        to_flux(grid[i].v, grid[i].F);
    }
    //Update the grid after first step note that the initial values of the conservative variables are stored in member u0 for each cell.
    interfaces(grid, ng);
    for (i=0; i<ng; i++){
        for (j=0; j<n; j++){
            grid[i].u[j]=(3./4.)*grid[i].u0[j]+(1./4.)*grid[i].u[j]+(1./4.)*delta_t*grid[i].L[j];
        }
        //Update fluxes and primitive variables in each cell
        to_prim(grid[i].u, grid[i].v);
        to_flux(grid[i].v, grid[i].F);
    }

    //Now take the final step
    interfaces(grid, ng);
    for (i=0; i<ng; i++){
        for (j=0; j<n; j++){
            grid[i].u[j]=(1./3.)*grid[i].u0[j]+(2./3.)*grid[i].u[j]+(2./3.)*delta_t*grid[i].L[j];
        }
        //Update fluxes and primitive variables in each cell
        to_prim(grid[i].u, grid[i].v);
        to_flux(grid[i].v, grid[i].F);
    }



}

//Evolves our grid of values forward one time step delta_t
void evolve (struct cell* grid,  int ng, double delta_t){
  //loop variable
  int i,j=0;

  ShuOsher(grid, ng, delta_t);

//  for (i=0; i<ng; i++){
//    for (j=0; j<n; j++){
//      //Update the conservative variables for now this is just simple Euler form.
//      grid[i].u[j]+=delta_t*grid[i].L[j];

//    }
    //Update fluxes and primitive variables in each cell
//    to_prim(grid[i].u, grid[i].v);
//    to_flux(grid[i].v, grid[i].F);
//}

}



int main(int argc, char* argv[]){
  //File to which we write output
  FILE* out=fopen("hydro.txt", "w");

  //Loop variable
  int i=0;
  //Number of grid points
  int ng=200;
  //x boundaries of our numerical domain
  double xmin=0;
  double xmax=1;
  //size of each of our cells
  double delta_x=(xmax-xmin)/ng;
  //Allocating our grid
  struct cell* grid=malloc((ng+4)*sizeof(struct cell));
  struct cell* grid2=&(grid[2]);



  //Initial condition
  double rho_l=1;
  double v_l=0;
  double p_l=1;
  double rho_r=0.1;
  double v_r=0;
  double p_r=0.125;
  //Initializing grid using our initial condition
  for(i=0; i<((ng+4)/2); i++){
    grid[i].v=malloc(n*sizeof(double));
    grid[i].u=malloc(n*sizeof(double));
    grid[i].F=malloc(n*sizeof(double));

    grid[i].v[0]=rho_l;
    grid[i].v[1]=v_l;
    grid[i].v[2]=p_l;

    to_cons(grid[i].v, grid[i].u);
    to_flux(grid[i].v, grid[i].F);

    grid[i].delta_x=delta_x;
  }

  for (i=((ng+4)/2); i<(ng+4); i++){
    grid[i].v=malloc(n*sizeof(double));
    grid[i].u=malloc(n*sizeof(double));
    grid[i].F=malloc(n*sizeof(double));

    grid[i].v[0]=rho_r;
    grid[i].v[1]=v_r;
    grid[i].v[2]=p_r;

    to_cons(grid[i].v, grid[i].u);
    to_flux(grid[i].v, grid[i].F);

    grid[i].delta_x=delta_x;
  }
  //The following lines should make bdry cells which behave like mirrors
  grid[1].v=grid[2].v;
  grid[0].v=grid[2].v;
  grid[ng].v=grid[ng-1].v;
  grid[ng+1].v=grid[ng-1].v;

  grid[1].u=grid[2].u;
  grid[0].u=grid[2].u;
  grid[ng].u=grid[ng-1].u;
  grid[ng+1].u=grid[ng-1].u;

  grid[1].F=grid[2].F;
  grid[0].F=grid[2].F;
  grid[ng].F=grid[ng-1].F;
  grid[ng+1].F=grid[ng-1].F;


  

  //Calculate fluxes at interfaces
  interfaces(grid2,ng);

  printf("test\n");

  //Calculate cfl time-step
  double delta_t=cfl(grid2, ng, delta_x);
  //Keeping track of time that went by
  double t=0;
  //Maximum time to which we would like to integrate
  double tmax=0.2;
  //Keep track of number of iterations
  int iterations=0;


  while((t<tmax)){
    //Calculate time step
    delta_t=cfl(grid2, ng, delta_x);
    //If cfl time-step is bigger than the  
    if (delta_t>(tmax-t)){
      delta_t=tmax-t;
    }

    //Calculate fluxes at interfaces
    //interfaces(grid2,ng);

    //Evolve forward one-time step
    evolve(grid2, ng, delta_t);
    t+=delta_t;
    iterations++;



  }

 //Write result to output file.
      for (i=0; i<ng; i++){
          fprintf(out, "%lf %lf %lf\n", t, delta_x*(i+(1./2.)), (grid2[i]).v[0]);
      }

  return 0;
}


