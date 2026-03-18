/*This is bohr's model of hydrogen atom
 * file name is simulation.c*/

#include <stdio.h>
#include <math.h>

#define p_mass 1.67262193e-27
#define e_mass 9.1093837e-31
#define p_charge 1.60217663e-19
#define e_charge -1.60217663e-19
#define k 8.98755e9

struct proton
{
  double x_position;
  double y_position;
  double vx;
  double vy;
  double ax;
  double ay;
};


struct electron {
  double x_position;
  double y_position;
  double vx;
  double vy;
  double ax;
  double ay;
  int state;
};

int main()
{
  struct electron e1;
  struct proton p;


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

  }
  
  fclose(fp);
  printf("Simulation complete. Data saved to 'orbit.csv'.\n");
  return 0;
}
