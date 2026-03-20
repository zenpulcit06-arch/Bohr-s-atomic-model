#include <math.h>
#include "acc.h"
#include "../particle/particle.h"
#define e_mass 9.1093837e-31
#define p_mass 1.67262193e-27

void acc(struct particle *p, struct particle *e,double e_ch, double p_ch )
{ 
  double r = hypot(e->x_position - p->x_position, e->y_position - p->y_position);
  double force =8.98755e9 * e_ch*p_ch/pow(r,2);

  
  e->ax = (force/e_mass)*((e -> x_position - p->x_position)/r);
  e->ay = (force/e_mass)*((e->y_position - p->y_position)/r);

  p->ax = (force/p_mass)*((p->x_position - e->x_position)/r);
  p->ay = (force/p_mass)*((p->y_position - e->y_position)/r);
}
