#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>


#define gamma 1.4
#define n 3
#define safety 0.9
#define theta 1.5


struct cell{
    //Values of flux, conservative variables, and primitive variables for a cell
    double F[n];
    double u[n];
    double v[n];

    //Size of cell
    double delta_x;

    //left and right biased values for each interface
    double F_ll[n], F_lr[n], F_rl[n], F_rr[n];
    double u_ll[n], u_lr[n], u_rl[n], u_rr[n];
    double v_ll[n], v_lr[n], v_rl[n], v_rr[n];

    double Fi_l[n], Fi_r[n];

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
    return 1;
}
//Minmod function
double minmod(double x, double y, double z){
    double tmp[3]={fabs(x), fabs(y), fabs(z)};


    return (1./4.)*abs(sgn(x)+sgn(y))*(sgn(x)+sgn(z))*min(tmp, 3);

}
//converts to primitive variables
void to_prim(const double u[], double v[]){

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

    //printf("%lf %lf %lf %lf\n", rho, vel, p, cs);

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


//void derivs(double t, const double u[], double L[], int ng, double delta_x){
//
//    //Loop variables
//    int i=0;
//    int j=0;
//
//
//    double F[ng*n];
//    //Calculating the flux at each of the grid points
//    flux(u, F, ng);
//
//    //Flux and u with 2 ghost cells appended
//    double F_ghost[(ng+4)*n];
//    double u_ghost[(ng+4)*n];
//    //Filling in our second flux vector with ghost cells
//    for(i=0; i<(ng+4)*n; i++){
//        //Boundary cells
//        if (i<2*n){
//            F_ghost[i]=F[i%3];
//            u_ghost[i]=u[i%3];
//        }
//        else if (i>((ng+2)*n)-1){
//            F_ghost[i]=F[(ng*n)-3+(i%3)];
//            u_ghost[i]=u[(ng*n)-3+(i%3)];
//        }
//        //Middle cells
//        else{
//            F_ghost[i]=F[i-(2*n)];
//            u_ghost[i]=u[i-(2*n)];
//        }
//    }
//
//
//
//
//
//
//
//    double aplus[ng+3], aminus[ng+3];
//    alpha_plus(u_ghost,  aplus, ng+4);
//    alpha_minus(u_ghost, aminus, ng+4);
//
//    //Loop variable
//    int k=0;
//
//    //For each of the grid points. N.B. The index i represents the index of the grid point.
//    for (i=0; i<ng; i++){
//        //For each of the three variables.
//        for (j=0; j<n; j++){
//            double Fi_down=hll(&(u_ghost[2*n]), &(F_ghost[2*n]), &(aplus[2]), &(aminus[2]), i-1, j, ng);
//            double Fi_up=hll(&(u_ghost[2*n]), &(F_ghost[2*n]), &(aplus[2]), &(aminus[2]), i, j, ng);
//            //Time derivative
//            L[(n*i)+j]=-(Fi_up-Fi_down)/(delta_x);
//            //printf("%d %lf %lf %lf\n", i+j, L[(n*i)+j], Fi_up, Fi_down);
//
//
//
//        }
//    }
//    //printf("\n");
//
//
//}

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

void ShuOsher(double t, double u[], double L[], double delta_t,  int ng, double delta_x, void(* myderivs)(double,  const double *, double *, int,  double)){
    //Loop variable
    int i=0;

    //Store solution vectors at intermediate steps
    double u1[ng*n], u2[ng*n];

    //Calculate derivatives
    myderivs(t, u, L, ng, delta_x);
    for (i=0; i<ng*n; i++){
        u1[i]=u[i]+delta_t*L[i];
    }
    /*Calculate derivatives with the new u
      Note that the explicit value of time
      passed to the derivatives function does not matter.
      Formally, I should probably pass t + delta_t or something*/
    myderivs(t, u1, L, ng, delta_x);
    for (i=0; i<(ng*n); i++){
        u2[i]=(3./4.)*u[i]+(1./4.)*u1[i]+(1./4.)*delta_t*L[i];
    }
    /*Calculate derivatives at u2. Note once again that the explicit
    value of t does not matter.*/
    myderivs(t, u2, L, ng, delta_x);
    for (i=0; i<(ng*n); i++){
        u[i]=(1./3.)*u[i]+(2./3.)*u2[i]+(2./3.)*delta_t*L[i];
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

//Given grid and its size this function computes the fluxes at all of the interfaces
void interfaces(struct cell* grid, int ng){
    int i=0;
    int j=0;

    for(i=1; i<ng-1; i++){
        for (j=0; j<n; j++){
            //printf("%lf\n", grid[i].v[j]);

            grid[i].v_ll[j]=grid[i-1].v[j];
            grid[i].v_lr[j]=grid[i].v[j];
            grid[i].v_rl[j]=grid[i].v[j];
            grid[i].v_rr[j]=grid[i+1].v[j];

            //printf("%lf %lf %lf\n", grid[i-1].v[j], grid[i].v_lr[j], grid[i].v_ll[j]);

            grid[i].u_ll[j]=grid[i-1].u[j];
            grid[i].u_lr[j]=grid[i].u[j];
            grid[i].u_rl[j]=grid[i].u[j];
            grid[i].u_rr[j]=grid[i+1].u[j];

            grid[i].F_ll[j]=grid[i-1].F[j];
            grid[i].F_lr[j]=grid[i].F[j];
            grid[i].F_rl[j]=grid[i].F[j];
            grid[i].F_rr[j]=grid[i+1].F[j];



        }



        grid[i].aplus_l=alpha_plus(grid[i].v_ll, grid[i].v_lr);
        grid[i].aplus_r=alpha_plus(grid[i].v_rl, grid[i].v_rr);
        grid[i].aminus_l=alpha_minus(grid[i].v_ll, grid[i].v_lr);
        grid[i].aminus_r=alpha_minus(grid[i].v_rl, grid[i].v_rr);



        hll(grid[i].Fi_l, grid[i].u_ll, grid[i].u_lr, grid[i].F_ll, grid[i].F_lr, grid[i].aplus_l, grid[i].aminus_l);
        hll(grid[i].Fi_r, grid[i].u_rl, grid[i].u_rr, grid[i].F_rl, grid[i].F_rr, grid[i].aplus_r, grid[i].aminus_r);


        //printf("%lf %lf\n", grid[i].aminus_l ,grid[i].Fi_l[0]);

    }


    for (j=0; j<n; j++){
        grid[0].Fi_r[j]=grid[1].Fi_l[j];
        grid[0].Fi_l[j]=grid[0].F[j];
    }

    for (j=0; j<n; j++){
        grid[ng-1].Fi_l[j]=grid[ng-2].Fi_r[j];
        grid[0].Fi_r[j]=grid[ng-1].F[j];
    }



}


int main(int argc, char* argv[]){
    FILE* out=fopen("hydro.txt", "w");

    //Loop variable
    int i=0;
    //Number of grid points
    int ng=10;
    //x boundaries of our numerical domain
    double xmin=0;
    double xmax=5;
    //size of each of our cells
    double delta_x=(xmax-xmin)/ng;

    //Allocating our grid
    struct cell* grid=malloc(ng*sizeof(struct cell));



    //Initial condition
    double rho_l=10;
    double v_l=0;
    double p_l=100;
    double rho_r=1;
    double v_r=0;
    double p_r=1;
    //Initializing grid using our initial condition
    for(i=0; i<((ng)/2); i++){
        grid[i].v[0]=rho_l;
        grid[i].v[1]=v_l;
        grid[i].v[2]=p_l;

        to_cons(grid[i].v, grid[i].u);
        to_flux(grid[i].v, grid[i].F);


    }
    printf("%d\n", i);
    for (i=(ng/2); i<(ng); i++){
        grid[i].v[0]=rho_r;
        grid[i].v[1]=v_r;
        grid[i].v[2]=p_r;

        to_cons(grid[i].v, grid[i].u);
        to_flux(grid[i].v, grid[i].F);
    }
    printf("test\n");
    interfaces(grid,ng);

    for (i=0; i<ng; i++){
        printf("%lf\n", grid[i].Fi_l[1]);
    }


//    //Alpha's
//    double aplus[ng-1], aminus[ng-1];
//    //Time derivatives
//    double L[ng*n];
//    //Time step we are going to take
//    double delta_t=0;
//    //Keeping track of time that went by
//    double t=0;
//    //Maximum time to which we would like to integrate
//    double tmax=0.4;
//    //Keep track of number of iterations
//    int iterations=0;
//
//    //Write the initial conditionto output file
//    for (i=0; i<ng*n; i+=n){
//        fprintf(out, "%lf %lf %lf\n", t, delta_x*(i/n), u[i]);
//    }
//    while((t<tmax)){
//        //Calculate alpha at each grid point
//        alpha_minus(u, aminus, ng);
//        alpha_plus(u, aplus, ng);
//        //Calculate time step
//        delta_t=cfl(aplus, aminus, ng, delta_x);
//        if (delta_t>(tmax-t)){
//            delta_t=tmax-t;
//        }
//
//
//
//
//        //Evolve forward one time step
//        Euler(t, u, L, delta_t, ng, delta_x, &derivs);
//        t+=delta_t;
//        iterations++;
//
//        //Write result to output file.
//        for (i=0; i<ng*n; i+=n){
//            fprintf(out, "%lf %lf %lf\n", t, delta_x*(i/n), u[i]);
//        }
//
//    }

    return 0;
}


