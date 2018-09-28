/*
 * particle_types.h
 *
 *  Created on: Apr 9, 2016
 *      Author: snytav
 */

#ifndef PARTICLE_TYPES_H_
#define PARTICLE_TYPES_H_

typedef enum {ION, PLASMA_ELECTRON, BEAM_ELECTRON} particle_sorts;

int particle_sort_name(char *name,particle_sorts sort);

typedef enum {DENSITY,KINETIC_ENERGY,VELOCITY,CURRENT
} diagnostics_operations;



#endif /* PARTICLE_TYPES_H_ */
