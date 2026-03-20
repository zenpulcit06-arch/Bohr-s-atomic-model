#include "integrator.h"
#include "../particle/particle.h"

void integratorfh(struct particle *p, struct particle *e, double dt)
{
  e->x_position += e->vx*dt + 0.5*e->ax*dt*dt;
  e->y_position += e->vy*dt + 0.5*e->ay*dt*dt;


  p->x_position += p->vx*dt + 0.5*p->ax*dt*dt;
  p->y_position += p->vy*dt + 0.5*p->ay*dt*dt;

    e->vx += 0.5 * e->ax * dt;
    e->vy += 0.5 * e->ay * dt;
    p->vx += 0.5 * p->ax * dt;
    p->vy += 0.5 * p->ay * dt;

}

void integratorlh(struct particle *p, struct particle *e, double dt)
{
   e->vx += 0.5 * (e->ax) * dt;
   e->vy += 0.5 * ( e->ay) * dt;
   p->vx += 0.5 * (p->ax) * dt;
   p->vy += 0.5 * ( p->ay) * dt;

}
