/*
 * Logger.cpp
 *
 *  Created on: Oct 15, 2015
 *      Author: thomas
 */

#include "Logger.h"
#include "gnuplot_i.hpp"

Logger::Logger(Integrator *integrator) {
	this->integrator = integrator;
	this->nbs = integrator->nbs;
	// Initialize log vectors.
	for (int i=0; i<this->nbs->N; i++) {
		X.push_back(std::vector<double> ());
		Y.push_back(std::vector<double> ());
		Z.push_back(std::vector<double> ());
	}
}

Logger::~Logger() {
	// TODO Auto-generated destructor stub
}

void Logger::iterate(int Niter, int step) {
	// Iterate Niter times and keep every step-th system status.
	for (int i=0; i<Niter; i++) {
		integrator->iterate();
		this->nbs->calc_energy();
		std::cout << this->nbs->E << std::endl;
		// Store every step-th set of positions.
		if (i%step == 0) {
			for (int j=0; j<nbs->N; j++) {
				X[j].push_back(this->nbs->rx[j]);
				Y[j].push_back(this->nbs->ry[j]);
				Z[j].push_back(this->nbs->rz[j]);
			}
		}
	}
}

void Logger::plot_trajectories() {
	// Plot the stored trajectories.
	Gnuplot gp;
	gp << "set size ratio -1\n";
	for (int i=0; i<this->nbs->N; i++) {
		std::stringstream title;
		title << "Particle " << i << " with m = " << this->nbs->m[i];
		gp.plot_xyz(X[i],Y[i],Z[i], title.str());
	}
}
