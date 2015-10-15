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

	std::vector<std::vector<double> > X;
	std::vector<std::vector<double> > Y;
	std::vector<std::vector<double> > Z;

	Logger(Integrator *integrator);
	virtual ~Logger();

	void iterate(int Niter, int step);
	void plot_trajectories();
};

#endif /* LOGGER_H_ */
