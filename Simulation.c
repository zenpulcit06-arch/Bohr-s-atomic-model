/*This is bohr's model of hydrogen atom
 * file name is simulation.c*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "sim_essen/particle/particle.h"

#define p_mass 1.67262193e-27
#define e_mass 9.1093837e-31
#define p_charge 1.60217663e-19
#define e_charge -1.60217663e-19
#define k 8.98755e9


int main()
{
  struct particle e1,p;

  srand(time(NULL));
  double decay_const = 1e13;
  int jumped = 0;


  printf("Enter state of electron\n");
  scanf("%d",&e1.state);

  e1.x_position = 0.53e-10*pow(e1.state,2);
  e1.y_position = 0;
  p.x_position = 0;
  p.y_position =0;
  
  double potential_energy,kinetic_energy,total_energy;
  double r ;
  double angle ;
  double dt = 1e-18;
  double Time = 2e-15*pow(e1.state,3);
  long long steps = (long long)(Time/dt);
  double force;

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
  r = hypot(e1.x_position - p.x_position,e1.y_position-p.y_position);

  force = (k*p_charge*e_charge)/pow(r,2);

  e1.ax = (force/e_mass)*((e1.x_position - p.x_position)/r);
  e1.ay = (force/e_mass)*((e1.y_position - p.y_position)/r);

  p.ax = (force/p_mass)*((p.x_position - e1.x_position)/r);
  p.ay = (force/p_mass)*((p.y_position - e1.y_position)/r);
   

  for (long long i = 0; i < steps; i++)
  {
    double probability = 1 - exp( decay_const * dt);
    double roll = (double)rand()/RAND_MAX;

    e1.x_position += e1.vx*dt + 0.5*e1.ax*dt*dt;
    e1.y_position += e1.vy*dt + 0.5*e1.ay*dt*dt;

    p.y_position += p.vy*dt + 0.5*p.ay*dt*dt;
    p.x_position += p.vx*dt + 0.5*p.ax*dt*dt;

    double e_ax_old = e1.ax; double e_ay_old = e1.ay;
    double p_ax_old = p.ax; double p_ay_old = p.ay;


   r = hypot(e1.x_position - p.x_position,e1.y_position - p.y_position);
   force = (k*p_charge*e_charge)/pow(r,2);

   e1.ax = (force/e_mass)*((e1.x_position - p.x_position)/r);
   e1.ay = (force/e_mass)*((e1.y_position - p.y_position)/r);

   p.ax = (force/p_mass)*((p.x_position - e1.x_position)/r);
   p.ay = (force/p_mass)*((p.y_position - e1.y_position)/r);
   
   e1.vx += 0.5 * (e_ax_old + e1.ax) * dt;
    e1.vy += 0.5 * (e_ay_old + e1.ay) * dt;
    p.vx  += 0.5 * (p_ax_old + p.ax) * dt;
    p.vy  += 0.5 * (p_ay_old + p.ay) * dt;


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

      e1.state = new_state;

      r = 0.53e-10 * new_state*new_state;
      e1.x_position = r;
      e1.y_position = 0;
      e1.vx = 0;
      e1.vy = 2.18e6/new_state;

      r = hypot(e1.x_position - p.x_position, e1.y_position - p.y_position);
      force = (k * p_charge * e_charge) / (r * r);

      e1.ax = (force / e_mass) * ((e1.x_position - p.x_position) / r);
      e1.ay = (force / e_mass) * ((e1.y_position - p.y_position) / r);

      p.ax = (force / p_mass) * ((p.x_position - e1.x_position) / r);
      p.ay = (force / p_mass) * ((p.y_position - e1.y_position) / r);

    }

  }
  
  fclose(fp);
  printf("Simulation complete. Data saved to 'orbit.csv'.\n");
  return 0;
}
