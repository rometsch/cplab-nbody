/*
 * nbodysys.h
 *
 *  Created on: Oct 12, 2015
 *      Author: thomas
 */

#ifndef NBSYS_
#define NBSYS_

#include <iostream>
#include <fstream>
#include <cmath>


class nbsys {
public:
	int N;		// Number of particles.
	double t_max;	// Maximal integration time.
	double eta;		// Time step parameter.
	double *rx, *ry, *rz, *vx, *vy, *vz, *ax, *ay, *az, *adotx, *adoty, *adotz;	// Position, velocity and acceleration components.
	double *m;		// Mass of particles.

	double Rx, Ry, Rz;		// Center of mass (COM) coordinate.
	double Pcom, Pcomx, Pcomy, Pcomz;	// Center of mass momentum
	double P, Px, Py, Pz;	// Total momentum with reference to COM system.
	double J, Jx, Jy, Jz;	// Total angular momentum.


	double E;		// Total energy.

	// Quantities for 2-body systems.
	double j, jx, jy, jz;	// Angular momentum.
	double e, ex, ey, ez;	// Runge-Lenz vector.
	double ae;		// Major axis of ellipse.

	nbsys(std::string filename);
	virtual ~nbsys();

	void calc_energy();			// Calculate the total energy.
	void calc_com();			// Calculate the center of mass, its velocity.
	void calc_momentum();		// Calculate the total momentum.
	void calc_angular_momentum();

	void calc_acc();	// Calculate the acceleration for all particles.
	void calc_dot_acc();	// Calculate the time derivative oft the acceleration for all particles.
	void acc(int i);	// Calculate the acceleration for particle i.
	void dot_acc(int i);	// Calculate time derivative of acceleration for particle i.

	void load_data(std::string filename);		// Load initial data from file.
	void transform_to_com();	// Transform the coordinates and velocities to
	void print_system();		// Print r,v,a for all particles.

};

#endif /* NBSYS_ */
