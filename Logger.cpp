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
	this->sim_time = 0;
}

Logger::~Logger() {
	// TODO Auto-generated destructor stub
}

void Logger::iterate(double time,int step) {
	// Iterate time (in units of t_max) and keep every step-th system status.
	integrator->start();
	int cnt = 0;
	while (this->sim_time < time*this->nbs->t_max) {
		integrator->iterate();
		this->sim_time += this->integrator->dt;
		// Store every step-th set of positions.
		if (++cnt%step == 0) {
			this->nbs->calc_energy();
			this->nbs->calc_two_body_vars();
			this->snapshot();
		}
	}
	integrator->stop();
}

void Logger::iterate(int step) {
	// Iterate Niter times and keep every step-th system status.
	integrator->start();
	int cnt = 0;
	while (this->sim_time < this->nbs->t_max) {
		integrator->iterate();
		this->sim_time += this->integrator->dt;
		this->nbs->calc_energy();
		this->nbs->calc_two_body_vars();
		// Store every step-th set of positions.
		if (++cnt%step == 0) {
			this->snapshot();
		}
	}
	integrator->stop();
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

	// Time
	this->simulation_time.push_back(sim_time);
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

void Logger::plot_system() {
	// Plot the stored trajectories.
	// First Copy data to vectors.
	std::vector<double> x,y,z;
	for (int i=0; i<this->nbs->N; i++) {
		x.push_back(this->nbs->rx[i]);
		y.push_back(this->nbs->ry[i]);
		z.push_back(this->nbs->rz[i]);
	}

	Gnuplot gp;
	gp << "set size ratio -1\n";
	std::stringstream title;
	title << this->nbs->N << "-body system" << std::endl;
	gp.plot_xyz(x, y, z, title.str());
}

void Logger::plot_drift(std::vector<double> &var, std::string var_name) {
	// Plot the difference of var at a given time to var at start as a function of time.
	std::vector<double> delta, t;
	this->calc_delta(var, delta, t);
//	std::cout << this->simulation_time.size() << std::endl;
//	std::cout << delta.size() << " , " << t.size() << std::endl;

	Gnuplot gp;
	std::stringstream title;
	title << var_name << " drift" << std::endl;
	gp << "set logscale y 10 \n";
	gp << "set format y \"%s*10^{%S}\"\n";
	gp.plot_xy(t, delta, title.str());
	gp << "pause mouse any \"Any key or button will terminate\"\n";
}

void Logger::plot_energy_drift() {
	// Plot drift of energy with time.
	this->plot_drift(this->E,"energy");
}

void Logger::plot_lenz_drift() {
	this->plot_drift(this->e, "lenz");
}

void Logger::output_drift(std::ostream &out, std::vector<double> &var, std::string var_name) {
	// Plot the difference of var at a given time to var at start as a function of time.
	std::vector<double> delta, t;
	this->calc_delta(var, delta, t);
	out << "# time\tdelta_" << var_name << std::endl;
	for (int i=0; i<t.size(); i++) {
		out << t[i] << "\t" << delta[i] << std::endl;
	}
}

void Logger::output_energy_drift(std::ostream &out) {
	this->output_drift(out,this->E, "energy");
}
void Logger::output_lenz_drift(std::ostream &out) {
	this->output_drift(out,this->e, "lenz");
}

void Logger::output_drift_two_body(std::ostream &out) {
	// Plot the difference of var at a given time to var at start as a function of time.
	std::vector<double> delta_E,delta_e, delta_ae, t, dummy;
	this->calc_delta(this->E, delta_E, t);
	this->calc_delta(this->e, delta_e, dummy);
	this->calc_delta(this->ae, delta_ae, dummy);
	out << "# time\tdelta_E\tdelta_e\tdelta_ae" << std::endl;
	for (int i=0; i<t.size(); i++) {
		out << t[i] << "\t" << delta_E[i] << "\t" << delta_e[i] << "\t" << delta_ae[i]  << std::endl;
	}
}

void Logger::calc_delta(std::vector<double> &var, std::vector<double> &delta, std::vector<double> &t) {
	// Calculate the difference between var[0] and all entries and save the simulation time for these entries in t.
	// Auxiliary function for plots and data save routines.
	double var0 = var[0];
	for (int i=1; i<this->simulation_time.size(); i++) {
		delta.push_back(fabs(var[i]-var0));
		t.push_back(this->simulation_time[i]);
	}
}
