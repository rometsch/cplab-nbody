/*
 * Logger.h
 *
 *  Created on: Oct 15, 2015
 *      Author: thomas
 */

#ifndef LOGGER_H_
#define LOGGER_H_


#include "Integrators.h"
#include <vector>



class Logger {
public:
	Integrator *integrator;
	nbsys *nbs;

	// Vector (1,..,N) of vectors to log each particles component.
	std::vector<std::vector<double> > X, Y, Z;

	// Vectors to to system variables (energy, momentum,...)
	std::vector<double> E,
	Rx, Ry, Rz,		// Center of mass (COM) coordinate.
	Pcom, Pcomx, Pcomy, Pcomz,	// Center of mass momentum
	P, Px, Py, Pz,	// Total momentum with reference to COM system.
	J, Jx, Jy, Jz;	// Total angular momentum.

	// Quantities for 2-body systems.
	std::vector<double> j, jx, jy, jz,	// Angular momentum.
	e, ex, ey, ez,	// Runge-Lenz vector.
	ae;		// Major axis of ellipse.

	// Time keepers.
	std::vector<double> simulation_time, real_time;
	double sim_time;

	Logger(Integrator *integrator);
	virtual ~Logger();

	void iterate(int Niter, int step);
	void snapshot();
	void plot_trajectories();
};

#endif /* LOGGER_H_ */
