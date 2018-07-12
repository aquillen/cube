/**
 * resolved mass spring model
 * using the leap frog integrator. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "rebound.h"
#include "tools.h"
#include "output.h"
#include "spring.h"


int NS;  // number of springs
struct spring* springs;  // springs structure
void reb_springs();// to pass springs to display
double surfdist;       // for identifying surface particles
double Pressure; // vertical pressure
void table_top();
void top_push();


double gamma_fac; // for adjustment of gamma of all springs
double t_damp;    // end faster damping, relaxation
double t_print;   // for printouts

char froot[30];   // output files

void heartbeat(struct reb_simulation* const r);

void additional_forces(struct reb_simulation* r){
   zero_accel(r);
   spring_forces(r); // spring forces
   top_push(r); // top force
   table_top(r); // bottom force 
}

int main(int argc, char* argv[]){
	struct reb_simulation* const r = reb_create_simulation();
	// Setup constants
	r->integrator	= REB_INTEGRATOR_LEAPFROG;
	r->gravity	= REB_GRAVITY_NONE;
	r->boundary	= REB_BOUNDARY_NONE;
	r->G 		= 0;		
        r->additional_forces = additional_forces;  // setup callback function for additional forces
        double mcube = 1.0;          // total mass of all particles
        double rcube = 1.0;          // length of cube 
        double tmax = 0.0;           // if 0 integrate forever

// things to set! ////////////////////// could be read in with parameter file
        double dt;             // timestep
        double b_distance;     // mininum interparticle spacing
        double mush_fac1;       // mush_fac*b_distance is maximum spring length
        double mush_fac2;       // 
        double mush_fac3;       // 
        double ggamma  =0.01;  // damping parm
        double ks1  = 0.0;      // spring constant
        double ks2  = 0.0;      // spring constants
        double ks3  = 0.0;      // spring constants

    if (argc ==1){
        strcpy(froot,"t1");   // to make output files
	dt	   = 1e-3;    // Timestep
        b_distance = 0.15;    // for creating random sphere, min separation between particles
        mush_fac1    = 2.3;    // ratio of smallest spring distance to minimum interparticle dist
        mush_fac2    = 2.3;    // ratio of smallest spring distance to minimum interparticle dist
        mush_fac3    = 2.3;    // ratio of smallest spring distance to minimum interparticle dist
        // spring damping
        ggamma   = 1.0;    // initial factor for initial damping value for springs
        t_damp   = 1.0;    // gamma to final values for all springs at this time
        ks1      = 8e-2;   // spring constant
        ks2      = 8e-2;   // spring constant
        ks3      = 8e-2;   // spring constant
        t_print =  100000.0;  // printouts 
        surfdist=0.1;         // for identifying surface particles
        Pressure=0.0;         // apply presure
     }
     else{
        FILE *fpi; // read in a parameter file
        fpi = fopen(argv[1],"r");
        char line[300];
        fgets(line,300,fpi);  sscanf(line,"%s",froot);
        fgets(line,300,fpi);  sscanf(line,"%lf",&dt);
        fgets(line,300,fpi);  sscanf(line,"%lf",&tmax);
        fgets(line,300,fpi);  sscanf(line,"%lf",&t_print);
        fgets(line,300,fpi);  sscanf(line,"%lf",&b_distance);
        fgets(line,300,fpi);  sscanf(line,"%lf",&mush_fac1);
        fgets(line,300,fpi);  sscanf(line,"%lf",&mush_fac2);
        fgets(line,300,fpi);  sscanf(line,"%lf",&mush_fac3);
        fgets(line,300,fpi);  sscanf(line,"%lf",&ks1);
        fgets(line,300,fpi);  sscanf(line,"%lf",&ks2);
        fgets(line,300,fpi);  sscanf(line,"%lf",&ks3);
        fgets(line,300,fpi);  sscanf(line,"%lf",&ggamma);
        fgets(line,300,fpi);  sscanf(line,"%lf",&gamma_fac);
        fgets(line,300,fpi);  sscanf(line,"%lf",&t_damp);
        fgets(line,300,fpi);  sscanf(line,"%lf",&surfdist);
        fgets(line,300,fpi);  sscanf(line,"%lf",&Pressure);

     printf("parm file read in\n");

     }
     

/// end of things to set /////////////////////////

        r->dt=dt;            // set integration timestep
	const double boxsize = 1.1*rcube;    // display window
	reb_configure_box(r,boxsize,1,1,1);
	r->softening      = b_distance/100.0;	// Gravitational softening length
// viewer +x to right, +y to up, z back and forth along line of sight


   struct spring spring_mush; // spring parameters for mush
   // properties of springs
   spring_mush.gamma          = ggamma; // damping coefficient
   spring_mush.ks             = ks1; // spring constant
   spring_mush.k_heat         = 1.0;    // heat transport coefficient
   double mush_distance1=b_distance*mush_fac1; 
   double mush_distance2=b_distance*mush_fac2; 
   double mush_distance3=b_distance*mush_fac3; 
       // distance for connecting and reconnecting springs

   FILE *fpr;
   char fname[200];
   sprintf(fname,"%s_run.txt",froot); // for simulation info
   fpr = fopen(fname,"w");

   NS=0; // start with no springs 


   // create particle distribution

   rand_rectangle(r,b_distance,rcube,rcube, rcube,mcube);
   
   // make springs, all pairs connected within interparticle distance mush_distance
   connect_springs_dist(r,mush_distance1, 0, r->N, spring_mush);
   spring_mush.ks = ks2; // add longer springs!
   connect_springs_dist(r,mush_distance2, 0, r->N, spring_mush);
   spring_mush.ks = ks3; // add longer springs!
   connect_springs_dist(r,mush_distance3, 0, r->N, spring_mush);
   
   reb_springs(r); // pass spring index list to display
   set_gamma_fac(1.0/gamma_fac);  // start with enhanced damping

   r->heartbeat = heartbeat;

   if (tmax ==0.0) // start the integration!!!!
      reb_integrate(r, INFINITY);
   else
      reb_integrate(r, tmax);
}


void heartbeat(struct reb_simulation* const r){
        static int index = 0;
	if (reb_output_check(r,10.0*r->dt)){
		reb_output_timing(r,0); // print time of simulation run on screen
	}
        if (fabs(r->t - t_damp) < 0.9*r->dt) set_gamma_fac(gamma_fac); 
            // damp initial bounce only , end damping
            // reset gamma only at t near t_damp
	
        // stuff to do every timestep
        // nothing!

        if (reb_output_check(r,t_print)) {
            write_particles(r,froot,index); //   output particle positions
            index++;
        }


}

// make a spring index list to pass to viewer
void reb_springs(struct reb_simulation* const r){
   r->NS = NS;
   r->springs_ii = malloc(NS*sizeof(int));
   r->springs_jj = malloc(NS*sizeof(int));
   for(int i=0;i<NS;i++){
     r->springs_ii[i] = springs[i].i;
     r->springs_jj[i] = springs[i].j;
   }
}



// I need to mark top surface in the first call as it can move!
void top_push(struct reb_simulation* const r){
    static int *surflist; // list of top surface  particles!
    static int Nsurf=0;  // number of particles on top surface
    static int first=0;
    if (first==0){ // first time the routine is called
       first=1;
       surflist = malloc(sizeof(int)*r->N);
       for(int i=0;i<r->N;i++){
           if (r->particles[i].y > 0.5 - surfdist){
              surflist[Nsurf] = i;
              Nsurf++;
           }
       }
       printf("Nsurf top =%d\n",Nsurf);
    }
    for(int j=0;j<Nsurf;j++){
       int i = surflist[j];
       r->particles[i].ay -= Pressure/Nsurf / r->particles[i].m;
    }
}

void table_top(struct reb_simulation* const r) {
    static int *surflist; // list of bottom surface  particles!
    static int Nsurf=0;  // number of particles on bottom surface
    static int first=0;
    if (first==0){ // first time the routine is called
       first=1;
       surflist = malloc(sizeof(int)*r->N);
       for(int i=0;i<r->N;i++){
           if (r->particles[i].y < -0.5 + surfdist){
              surflist[Nsurf] = i;
              Nsurf++;
           }
       }
       printf("Nsurf bottom=%d\n",Nsurf);
    }
    for(int j=0;j<Nsurf;j++){
       int i = surflist[j];
       r->particles[i].ay = 0.0;
    }
}



