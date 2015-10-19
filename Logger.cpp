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
			this->snapshot();
		}
	}
}

void Logger::snapshot() {
	// Store the state of the system.
	for (int j=0; j<nbs->N; j++) {
		this->X[j].push_back(this->nbs->rx[j]);
		this->Y[j].push_back(this->nbs->ry[j]);
		this->Z[j].push_back(this->nbs->rz[j]);
	}
	// Vectors to to system variables (energy, momentum,...)
	this->E.push_back(this->nbs->E);
	this->Rx.push_back(this->nbs->Rx);
	this->Ry.push_back(this->nbs->Ry);
	this->Rz.push_back(this->nbs->Rz);
	this->Pcom.push_back(this->nbs->Pcom);
	this->Pcomx.push_back(this->nbs->Pcomx);
	this->Pcomy.push_back(this->nbs->Pcomy);
	this->Pcomz.push_back(this->nbs->Pcomz);
	this->P.push_back(this->nbs->P);
	this->Px.push_back(this->nbs->Px);
	this->Py.push_back(this->nbs->Py);
	this->Pz.push_back(this->nbs->Pz);
	this->J.push_back(this->nbs->J);
	this->Jx.push_back(this->nbs->Jx);
	this->Jy.push_back(this->nbs->Jy);
	this->Jz.push_back(this->nbs->Jz);

	this->j.push_back(this->nbs->j);
	this->jx.push_back(this->nbs->jx);
	this->jy.push_back(this->nbs->jy);
	this->jz.push_back(this->nbs->jz);
	this->e.push_back(this->nbs->e);
	this->ex.push_back(this->nbs->ex);
	this->ey.push_back(this->nbs->ey);
	this->ez.push_back(this->nbs->ez);
	this->ae.push_back(this->nbs->ae);

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
