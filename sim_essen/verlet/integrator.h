#ifndef INTEGRATOR
#define INTEGRATOR
#include "../particle/particle.h"

void integratorfh(struct particle *p, struct particle *e, double dt);
void integratorlh(struct particle *p, struct particle *e, double dt);

#endif // !INTEGRATOR

