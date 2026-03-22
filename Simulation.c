/*This is bohr's model of hydrogen atom
 * file name is simulation.c*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "sim_essen/particle/particle.h"
#include "sim_essen/acc_cal/acc.h"
#include "sim_essen/verlet/integrator.h"

#define p_mass 1.67262193e-27
#define e_mass 9.1093837e-31
#define p_charge 1.60217663e-19
#define e_charge -1.60217663e-19
#define k 8.98755e9
#define c 299792458
#define h 6.62607015e-34

struct photon
         {
  double x;
  double E;
  double lamda;
};


int main()
{
  struct particle e1,p;
  struct photon pht;

  srand(time(NULL));
  double decay_const = 1e13;
  int jumped = 0;
  int photon = 0;


  printf("Enter state of electron\n");
  scanf("%d",&e1.state);

  e1.x_position = 0.53e-10*pow(e1.state,2);
  e1.y_position = 0;
  p.x_position = 0;
  p.y_position =0;
  
  double potential_energy,kinetic_energy,total_energy;
  double r ;
  double angle ;
  double dt = 1e-19;
  double Time = 2e-15*pow(e1.state,3);
  long long steps = (long long)(Time/dt);

  e1.vx = 0;
  e1.vy = 2.18e6/e1.state;
  p.vx = 0;
  p.vy = -(e_mass/p_mass)*2.18e6/e1.state;

  FILE *fp = fopen("orbit.csv", "w");
  if (fp == NULL) {
      printf("Error opening file!\n");
      return 1;
  }

  fprintf(fp, "step,ex,ey,px,py,te,ke,pe\n");

  FILE *phton = fopen("photon.csv","w");
  if (phton == NULL){
    printf("Error opening file\n");
    return 1;
  }
  fprintf(phton,"x,E,wl\n");

  acc(&p,&e1,e_charge,p_charge);
   

  for (long long i = 0; i < steps; i++)
  {
    double probability = 1 - exp( -decay_const * dt);
    double roll = (double)rand()/RAND_MAX;

    integratorfh(&p,&e1,dt);


   r = hypot(e1.x_position - p.x_position,e1.y_position - p.y_position);
   acc(&p,&e1,e_charge,p_charge);
   
   
    integratorlh(&p, &e1, dt);

   potential_energy = (k * p_charge * e_charge) / r;
   kinetic_energy = (0.5 * e_mass*e1.vx*e1.vx + 0.5*e_mass*e1.vy*e1.vy) + (0.5 * p_mass*p.vx*p.vx + 0.5*p_mass*p.vy*p.vy);

   total_energy = kinetic_energy + potential_energy;

   fprintf(fp, "%lld,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e\n", i, e1.x_position, e1.y_position, p.x_position, p.y_position,total_energy,kinetic_energy,potential_energy);

   if (i % 500 == 0) {
        printf("%lld | %.3e | %.2f eV\n", i, r, total_energy / 1.602e-19);
    }

    if (!jumped && roll < probability && e1.state >1)
    {
      jumped = 1;
      int new_state = e1.state -1;

      double Ei = -13.6/(e1.state * e1.state);
      double Ef = -13.6/(new_state * new_state);

      double photon_energy = Ei - Ef;

      printf("Transition %d -> %d | Photon = %.2f eV \n",e1.state,new_state,photon_energy);


      r = 0.53e-10 * new_state*new_state;
      double r_old = hypot(e1.x_position - p.x_position,e1.y_position - p.y_position);
      double v_old = hypot(e1.vx,e1.vy);
      double v_new = 2.18e6/new_state;
      double dx = e1.x_position - p.x_position;
      double dy = e1.y_position - p.y_position;

      e1.x_position = p.x_position + (r)/r_old*dx;
      e1.y_position = p.y_position + (r)/r_old*dy;
      e1.vx *= v_new / v_old;
      e1.vy *= v_new/v_old;

      acc(&p,&e1,e_charge,p_charge);
      photon = 1;
      pht.x = e1.x_position;
      pht.E = photon_energy;
      pht.lamda = ((h*c)/(photon_energy*1.602e-19));
      
      e1.state = new_state;


    }

    if (photon)
    {
      fprintf(phton,"%lf,%lf,%lf \n",pht.x,pht.E,pht.lamda);
      if(pht.x >= 0 )
      {
        pht.x += 0.01*c*dt;
      }
      else {
        pht.x -= 0.01*c*dt;
      }
    }


  }
  
  fclose(fp);
  fclose(phton);
  printf("Simulation complete. Data saved to 'orbit.csv'.\n");
  return 0;
}
