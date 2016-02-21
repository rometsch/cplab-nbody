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
#include <sstream>

#include <iostream>


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

	void iterate(double time, int step);		// Iterate Niter times and save every step-th system status.
	void iterate(int step);		// Iterate until time of the system reached tmax. Time means the time of the simulated system rather than real time. Save every step-th system status.
	void snapshot();
	void plot_trajectories();
	void plot_system();
	void plot_drift(std::vector<double> &var, std::string var_name);
	void plot_energy_drift();
	void plot_lenz_drift();
	void output_drift(std::ostream &out, std::vector<double> &var, std::string var_name);
	void output_energy_drift(std::ostream &out);
	void output_lenz_drift(std::ostream &out);
	void output_drift_two_body(std::ostream &out);


	void calc_delta(std::vector<double> &var, std::vector<double> &delta, std::vector<double> &t);
};

#endif /* LOGGER_H_ */
