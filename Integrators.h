/*
 * Integrators.h
 *
 *  Created on: Oct 15, 2015
 *      Author: thomas
 */

#ifndef INTEGRATORS_H_
#define INTEGRATORS_H_

#include "nbsys.h"

class Integrator {
public:
	nbsys *nbs;		// N-body system.
	double dt;
	double *rx_next, *ry_next, *rz_next, *vx_next, *vy_next, *vz_next; 	// Position, velocity and acceleration components.

	Integrator(std::string filename);
	virtual ~Integrator();

	virtual void start() = 0;			// Start the integrator. E.g. jumpstart for Leap Frog.
	virtual void stop() = 0;			// Stop the integrator.	E.g. stop for Leap Frog.
	virtual void iterate() = 0;		// Iterate once.
};

class Euler: public Integrator {
public:

	Euler(std::string filename, double dt);
	virtual ~Euler();

	void start();
	void stop();
	void iterate();
	void iterate_ind(int i);
};

class Euler_Comer: public Integrator {
public:

	Euler_Comer(std::string filename, double dt);
	virtual ~Euler_Comer();

	void start();
	void stop();
	void iterate();
	void iterate_ind(int i);
};

class Mittelung: public Integrator {
public:

	Mittelung(std::string filename, double dt);
	virtual ~Mittelung();

	void start();
	void stop();
	void iterate();
	void iterate_ind(int i);
};

class LeapFrog: public Integrator {
public:

	LeapFrog(std::string filename, double dt);
	virtual ~LeapFrog();

	void start();
	void stop();
	void iterate();
	void iterate_ind(int i);
};

class Verlet: public Integrator {
public:

	double *rx_last, *ry_last, *rz_last; 	// Position, velocity and acceleration components.

	Verlet(std::string filename, double dt);
	virtual ~Verlet();

	void start();
	void stop();
	void iterate();
	void iterate_ind(int i);
};

class HermitePC: public Integrator {
public:

	double *rx_p, *ry_p, *rz_p, *vx_p, *vy_p, *vz_p, *ax_p, *ay_p, *az_p, *adotx_p, *adoty_p, *adotz_p; 	// Position, velocity and acceleration components.
	double *a2x, *a2y, *a2z, *a3x, *a3y, *a3z; 	// Higher time derivatives of a.

	HermitePC(std::string filename, double dt);
	virtual ~HermitePC();

	void start();
	void stop();
	void iterate();
	void iterate_ind(int i);
	void iterate_correction_ind(int i);
	void calc_acc_derivatives(int i);
};

class HermiteIterated: public Integrator {
public:
	int iteration_depth;
	double *rx_p, *ry_p, *rz_p, *vx_p, *vy_p, *vz_p, *ax_p, *ay_p, *az_p, *adotx_p, *adoty_p, *adotz_p; 	// Position, velocity and acceleration components.
	double *a2x, *a2y, *a2z, *a3x, *a3y, *a3z; 	// Higher time derivatives of a.

	HermiteIterated(std::string filename, double dt, int iteration_depth);
	virtual ~HermiteIterated();

	void start();
	void stop();
	void iterate();
	void iterate_ind(int i);
	void iterate_correction_ind(int i);
	void calc_acc_derivatives(int i);
};

class RungeKutta: public Integrator {
public:

	// Variables for temporary values.
	double *rx_t, *ry_t, *rz_t, *vx_t, *vy_t, *vz_t, *ax_t, *ay_t, *az_t;
	double *rx_aux, *ry_aux, *rz_aux;

	RungeKutta(std::string filename, double dt);
	virtual ~RungeKutta();

	void start();
	void stop();
	void iterate();
	void iterate_ind(int i);
};


#endif /* INTEGRATORS_H_ */
