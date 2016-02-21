/*
 * Integrators.cpp
 *
 *  Created on: Oct 15, 2015
 *      Author: thomas
 */

#include "Integrators.h"

Integrator::Integrator(std::string filename) {
	// Create new N-body system.
	this->nbs = new nbsys::nbsys(filename);
	this->dt = 1;

	// Allocate new arrays.
	// Position.
	this->rx_next = new double[this->nbs->N];
	this->ry_next = new double[this->nbs->N];
	this->rz_next = new double[this->nbs->N];
	// Velocity.`
	this->vx_next = new double[this->nbs->N];
	this->vy_next = new double[this->nbs->N];
	this->vz_next = new double[this->nbs->N];

	for (int i=0; i<this->nbs->N; i++) {
		this->rx_next[i] = 0;
		this->ry_next[i] = 0;
		this->rz_next[i] = 0;
		this->vx_next[i] = 0;
		this->vy_next[i] = 0;
		this->vz_next[i] = 0;
	}
}

Integrator::~Integrator() {
	// Delete all dynamically generated objects.
	delete this->nbs;
}

void Integrator::set_dt(double dt) {
	// Set dt value.
	this->dt = dt;
}

Euler::Euler(std::string filename, double dt)
	: Integrator(filename) {
	// Create new Integrator.
	this->dt = dt;
	this->name = "Euler";
}

Euler::~Euler() {

}

void Euler::start() {}	//	Nothing to do.
void Euler::stop() {}	//	Nothing to do.

void Euler::iterate_ind(int i) {
	// Iterate velocity.
	this->vx_next[i] = this->nbs->vx[i] + this->nbs->ax[i]*this->dt;
	this->vy_next[i] = this->nbs->vy[i] + this->nbs->ay[i]*this->dt;
	this->vz_next[i] = this->nbs->vz[i] + this->nbs->az[i]*this->dt;
	// Iterate position.
	this->rx_next[i] = this->nbs->rx[i] + this->nbs->vx[i]*this->dt;
	this->ry_next[i] = this->nbs->ry[i] + this->nbs->vy[i]*this->dt;
	this->rz_next[i] = this->nbs->rz[i] + this->nbs->vz[i]*this->dt;
}

void Euler::iterate() {
	// Iterate all particles.
	this->nbs->calc_acc();
	for (int i=0; i<this->nbs->N; i++) {
		this->iterate_ind(i);
	}
	// Store the new positions and velocities in the system.
	for (int i=0; i<this->nbs->N; i++) {
		this->nbs->rx[i] = this->rx_next[i];
		this->nbs->ry[i] = this->ry_next[i];
		this->nbs->rz[i] = this->rz_next[i];
		this->nbs->vx[i] = this->vx_next[i];
		this->nbs->vy[i] = this->vy_next[i];
		this->nbs->vz[i] = this->vz_next[i];
	}
}

EulerComer::EulerComer(std::string filename, double dt)
	: Integrator(filename) {
	// Create new Integrator.
	this->dt = dt;
	this->name = "EulerComer";
}

EulerComer::~EulerComer() {

}

void EulerComer::start() {}	//	Nothing to do.
void EulerComer::stop() {}	//	Nothing to do.


void EulerComer::iterate_ind(int i) {
	// Iterate velocity.
	this->vx_next[i] = this->nbs->vx[i] + this->nbs->ax[i]*this->dt;
	this->vy_next[i] = this->nbs->vy[i] + this->nbs->ay[i]*this->dt;
	this->vz_next[i] = this->nbs->vz[i] + this->nbs->az[i]*this->dt;
	// Iterate position.
	this->rx_next[i] = this->nbs->rx[i] + this->vx_next[i]*this->dt;
	this->ry_next[i] = this->nbs->ry[i] + this->vy_next[i]*this->dt;
	this->rz_next[i] = this->nbs->rz[i] + this->vz_next[i]*this->dt;
}

void EulerComer::iterate() {
	// Iterate all particles.
	this->nbs->calc_acc();
	for (int i=0; i<this->nbs->N; i++) {
		this->iterate_ind(i);
	}
	// Store the new positions and velocities in the system.
	for (int i=0; i<this->nbs->N; i++) {
		this->nbs->rx[i] = this->rx_next[i];
		this->nbs->ry[i] = this->ry_next[i];
		this->nbs->rz[i] = this->rz_next[i];
		this->nbs->vx[i] = this->vx_next[i];
		this->nbs->vy[i] = this->vy_next[i];
		this->nbs->vz[i] = this->vz_next[i];
	}
}


Mittelung::Mittelung(std::string filename, double dt)
	: Integrator(filename) {
	// Create new Integrator.
	this->dt = dt;
	this->name = "Mittelung";
}

Mittelung::~Mittelung() {

}

void Mittelung::start() {}	//	Nothing to do.
void Mittelung::stop() {}	//	Nothing to do.



void Mittelung::iterate_ind(int i) {
	// Iterate velocity.
	this->vx_next[i] = this->nbs->vx[i] + this->nbs->ax[i]*this->dt;
	this->vy_next[i] = this->nbs->vy[i] + this->nbs->ay[i]*this->dt;
	this->vz_next[i] = this->nbs->vz[i] + this->nbs->az[i]*this->dt;
	// Iterate position.
	this->rx_next[i] = this->nbs->rx[i] + 0.5*(this->vx_next[i]+this->nbs->vx[i])*this->dt;
	this->ry_next[i] = this->nbs->ry[i] + 0.5*(this->vy_next[i]+this->nbs->vy[i])*this->dt;
	this->rz_next[i] = this->nbs->rz[i] + 0.5*(this->vz_next[i]+this->nbs->vz[i])*this->dt;
}

void Mittelung::iterate() {
	// Iterate all particles.
	this->nbs->calc_acc();
	for (int i=0; i<this->nbs->N; i++) {
		this->iterate_ind(i);
	}
	// Store the new positions and velocities in the system.
	for (int i=0; i<this->nbs->N; i++) {
		this->nbs->rx[i] = this->rx_next[i];
		this->nbs->ry[i] = this->ry_next[i];
		this->nbs->rz[i] = this->rz_next[i];
		this->nbs->vx[i] = this->vx_next[i];
		this->nbs->vy[i] = this->vy_next[i];
		this->nbs->vz[i] = this->vz_next[i];
	}
}


LeapFrog::LeapFrog(std::string filename, double dt)
	: Integrator(filename) {
	// Create new Integrator.
	this->dt = dt;
	this->name = "LeapFrog";
}

LeapFrog::~LeapFrog() {

}

void LeapFrog::iterate_ind(int i) {
	// Jump half time.

	// Iterate velocity.
	this->vx_next[i] = this->nbs->vx[i] + this->nbs->ax[i]*this->dt;
	this->vy_next[i] = this->nbs->vy[i] + this->nbs->ay[i]*this->dt;
	this->vz_next[i] = this->nbs->vz[i] + this->nbs->az[i]*this->dt;
	// Iterate position (will be net half time positions).
	this->rx_next[i] = this->nbs->rx[i] + this->vx_next[i]*this->dt;
	this->ry_next[i] = this->nbs->ry[i] + this->vy_next[i]*this->dt;
	this->rz_next[i] = this->nbs->rz[i] + this->vz_next[i]*this->dt;
}

void LeapFrog::start() {
	// Calculate the positions after half the timestep to get the firt half time values.
	for (int i=0; i<this->nbs->N; i++) {
		this->nbs->rx[i] += 0.5*this->nbs->vx[i]*this->dt;
		this->nbs->ry[i] += 0.5*this->nbs->vy[i]*this->dt;
		this->nbs->rz[i] += 0.5*this->nbs->vz[i]*this->dt;
	}
}

void LeapFrog::stop() {
	// Last step in position calculation. Is actually the same like jump start iteration.
	this->start();
}

void LeapFrog::iterate() {

	// Calculate accelerations.
	this->nbs->calc_acc();

	// Iterate all particles.
	for (int i=0; i<this->nbs->N; i++) {
		this->iterate_ind(i);
	}
	// Store the new positions and velocities in the system.
	for (int i=0; i<this->nbs->N; i++) {
		this->nbs->rx[i] = this->rx_next[i];
		this->nbs->ry[i] = this->ry_next[i];
		this->nbs->rz[i] = this->rz_next[i];
		this->nbs->vx[i] = this->vx_next[i];
		this->nbs->vy[i] = this->vy_next[i];
		this->nbs->vz[i] = this->vz_next[i];
	}
}


Verlet::Verlet(std::string filename, double dt)
	: Integrator(filename) {
	// Create new Integrator.
	this->name = "Verlet";
	this->dt = dt;
	this->rx_last = new double[this->nbs->N];
	this->ry_last = new double[this->nbs->N];
	this->rz_last = new double[this->nbs->N];
	for (int i=0; i<this->nbs->N; i++) {
		this->rx_last[i] = 0;
		this->ry_last[i] = 0;
		this->rz_last[i] = 0;
	}
}

Verlet::~Verlet() {
	 delete[] rx_last;
	 delete[] ry_last;
	 delete[] rz_last;
}

void Verlet::iterate_ind(int i) {
	// Iterate position.
	this->rx_next[i] = 2*this->nbs->rx[i] - this->rx_last[i] + this->nbs->ax[i]*this->dt*this->dt;
	this->ry_next[i] = 2*this->nbs->ry[i] - this->ry_last[i] + this->nbs->ay[i]*this->dt*this->dt;
	this->rz_next[i] = 2*this->nbs->rz[i] - this->rz_last[i] + this->nbs->az[i]*this->dt*this->dt;

	// Calculate velocities. Note that velocity is always one time step behind.
	// Velocities for last timestep will be calculated in stop();
	this->vx_next[i] = (this->rx_next[i] - this->rx_last[i])/2/this->dt;
	this->vy_next[i] = (this->ry_next[i] - this->ry_last[i])/2/this->dt;
	this->vz_next[i] = (this->rz_next[i] - this->rz_last[i])/2/this->dt;
}


void Verlet::start() {
	// Calculate the positions half the timestep ago to get the firt half time values.
	for (int i=0; i<this->nbs->N; i++) {
		this->rx_last[i] = this->nbs->rx[i] - this->nbs->vx[i]*this->dt + 0.5*this->nbs->ax[i]*this->dt*this->dt;
		this->ry_last[i] = this->nbs->ry[i] - this->nbs->vy[i]*this->dt + 0.5*this->nbs->ay[i]*this->dt*this->dt;
		this->rz_last[i] = this->nbs->rz[i] - this->nbs->vz[i]*this->dt + 0.5*this->nbs->az[i]*this->dt*this->dt;
	}
}

void Verlet::stop() {
	// Calculate the velocities for the last step.
	// Just iterate all particles but only save the velocities.
	this->nbs->calc_acc();
	for (int i=0; i<this->nbs->N; i++) {
		this->iterate_ind(i);
	}
	for (int i=0; i<this->nbs->N; i++) {
		this->nbs->vx[i] = this->vx_next[i];
		this->nbs->vy[i] = this->vy_next[i];
		this->nbs->vz[i] = this->vz_next[i];
	}

}

void Verlet::iterate() {

	// Calculate accelerations.
	this->nbs->calc_acc();

	// Iterate all particles.
	for (int i=0; i<this->nbs->N; i++) {
		this->iterate_ind(i);
	}
	// Keep the current positions as last positions for next iteration.
	for (int i=0; i<this->nbs->N; i++) {
		this->rx_last[i] = this->nbs->rx[i];
		this->ry_last[i] = this->nbs->ry[i];
		this->rz_last[i] = this->nbs->rz[i];
	}

	// Store the new positions and velocities in the system.
	for (int i=0; i<this->nbs->N; i++) {
		this->nbs->rx[i] = this->rx_next[i];
		this->nbs->ry[i] = this->ry_next[i];
		this->nbs->rz[i] = this->rz_next[i];
		this->nbs->vx[i] = this->vx_next[i];
		this->nbs->vy[i] = this->vy_next[i];
		this->nbs->vz[i] = this->vz_next[i];
	}
}



HermitePC::HermitePC(std::string filename, double dt)
	: Integrator(filename) {
	// Create new Integrator.
	this->name = "HermitePC";
	this->dt = dt;
	this->rx_p = new double[this->nbs->N];
	this->ry_p = new double[this->nbs->N];
	this->rz_p = new double[this->nbs->N];
	this->vx_p = new double[this->nbs->N];
	this->vy_p = new double[this->nbs->N];
	this->vz_p = new double[this->nbs->N];
	this->ax_p = new double[this->nbs->N];
	this->ay_p = new double[this->nbs->N];
	this->az_p = new double[this->nbs->N];
	this->adotx_p = new double[this->nbs->N];
	this->adoty_p = new double[this->nbs->N];
	this->adotz_p = new double[this->nbs->N];
	this->a2x = new double[this->nbs->N];
	this->a2y = new double[this->nbs->N];
	this->a2z = new double[this->nbs->N];
	this->a3x = new double[this->nbs->N];
	this->a3y = new double[this->nbs->N];
	this->a3z = new double[this->nbs->N];

	for (int i=0; i<this->nbs->N; i++) {
		this->rx_p[i] = 0;
		this->ry_p[i] = 0;
		this->rz_p[i] = 0;
		this->vx_p[i] = 0;
		this->vy_p[i] = 0;
		this->vz_p[i] = 0;
		this->ax_p[i] = 0;
		this->ay_p[i] = 0;
		this->az_p[i] = 0;
		this->adotx_p[i] = 0;
		this->adoty_p[i] = 0;
		this->adotz_p[i] = 0;
		this->a2x[i] = 0;
		this->a2y[i] = 0;
		this->a2z[i] = 0;
		this->a3x[i] = 0;
		this->a3y[i] = 0;
		this->a3z[i] = 0;
	}
}

HermitePC::~HermitePC() {
	 delete[] rx_p;
	 delete[] ry_p;
	 delete[] rz_p;
	 delete[] vx_p;
	 delete[] vy_p;
	 delete[] vz_p;
	 delete[] ax_p;
	 delete[] ay_p;
	 delete[] az_p;
	 delete[] adotx_p;
	 delete[] adoty_p;
	 delete[] adotz_p;
	 delete[] a2x;
	 delete[] a2y;
	 delete[] a2z;
	 delete[] a3x;
	 delete[] a3y;
	 delete[] a3z;
}

void HermitePC::iterate_ind(int i) {
	// Calculate predicted velocities.
	this->vx_p[i] = this->nbs->vx[i] + this->nbs->ax[i]*this->dt + 0.5*this->nbs->adotx[i]*this->dt*this->dt;
	this->vy_p[i] = this->nbs->vy[i] + this->nbs->ay[i]*this->dt + 0.5*this->nbs->adoty[i]*this->dt*this->dt;
	this->vz_p[i] = this->nbs->vz[i] + this->nbs->az[i]*this->dt + 0.5*this->nbs->adotz[i]*this->dt*this->dt;
	// Calculate predicted positions.
	this->rx_p[i] = this->nbs->rx[i] + this->nbs->vx[i]*this->dt + 0.5*this->nbs->ax[i]*this->dt*this->dt + this->nbs->adotx[i]*this->dt*this->dt*this->dt/6;
	this->ry_p[i] = this->nbs->ry[i] + this->nbs->vy[i]*this->dt + 0.5*this->nbs->ay[i]*this->dt*this->dt + this->nbs->adoty[i]*this->dt*this->dt*this->dt/6;
	this->rz_p[i] = this->nbs->rz[i] + this->nbs->vz[i]*this->dt + 0.5*this->nbs->az[i]*this->dt*this->dt + this->nbs->adotz[i]*this->dt*this->dt*this->dt/6;
}

void HermitePC::iterate_correction_ind(int i) {
	// Calculate the correction values. Store corrected values in _next variables.

	// Calculate higher derivatives of a.
	this->calc_acc_derivatives(i);

	this->vx_next[i] = this->vx_p[i] + this->a2x[i]*this->dt*this->dt*this->dt/6 + this->a3x[i]*this->dt*this->dt*this->dt*this->dt/24;
	this->vy_next[i] = this->vy_p[i] + this->a2y[i]*this->dt*this->dt*this->dt/6 + this->a3y[i]*this->dt*this->dt*this->dt*this->dt/24;
	this->vz_next[i] = this->vz_p[i] + this->a2z[i]*this->dt*this->dt*this->dt/6 + this->a3z[i]*this->dt*this->dt*this->dt*this->dt/24;

	this->rx_next[i] = this->rx_p[i] + this->a2x[i]*this->dt*this->dt*this->dt*this->dt/24 + this->a3x[i]*this->dt*this->dt*this->dt*this->dt*this->dt/120;
	this->ry_next[i] = this->ry_p[i] + this->a2y[i]*this->dt*this->dt*this->dt*this->dt/24 + this->a3y[i]*this->dt*this->dt*this->dt*this->dt*this->dt/120;
	this->rz_next[i] = this->rz_p[i] + this->a2z[i]*this->dt*this->dt*this->dt*this->dt/24 + this->a3z[i]*this->dt*this->dt*this->dt*this->dt*this->dt/120;
}

void HermitePC::calc_acc_derivatives(int i) {
	// Calculate higher derivatives of a.
	this->a2x[i] = -6*(this->nbs->ax[i]-this->ax_p[i])/this->dt/this->dt - 2*(2*this->nbs->adotx[i]+this->adotx_p[i])/this->dt;
	this->a2y[i] = -6*(this->nbs->ay[i]-this->ay_p[i])/this->dt/this->dt - 2*(2*this->nbs->adoty[i]+this->adoty_p[i])/this->dt;
	this->a2z[i] = -6*(this->nbs->az[i]-this->az_p[i])/this->dt/this->dt - 2*(2*this->nbs->adotz[i]+this->adotz_p[i])/this->dt;

	this->a3x[i] = 12*(this->nbs->ax[i]-this->ax_p[i])/this->dt/this->dt/this->dt + 6*(this->nbs->adotx[i]+this->adotx_p[i])/this->dt/this->dt;
	this->a3y[i] = 12*(this->nbs->ay[i]-this->ay_p[i])/this->dt/this->dt/this->dt + 6*(this->nbs->adoty[i]+this->adoty_p[i])/this->dt/this->dt;
	this->a3z[i] = 12*(this->nbs->az[i]-this->az_p[i])/this->dt/this->dt/this->dt + 6*(this->nbs->adotz[i]+this->adotz_p[i])/this->dt/this->dt;
}


void HermitePC::start() {} // Nothing to do.

void HermitePC::stop() {} // Nothing to do.

void HermitePC::iterate() {

	// Calculate accelerations.
	this->nbs->calc_acc();
	this->nbs->calc_dot_acc();

	// Iterate all particles.
	for (int i=0; i<this->nbs->N; i++) {
		this->iterate_ind(i);
	}
	// Calculate new accelerations and time derivatives of them using predicted positions and velocities.
	for (int i=0; i<this->nbs->N; i++) {
		this->nbs->acc(i, this->rx_p, this->ry_p, this->rz_p, this->ax_p, this->ay_p, this->az_p);
		this->nbs->dot_acc(i, this->rx_p, this->ry_p, this->rz_p, this->vx_p, this->vy_p, this->vz_p, this->adotx_p, this->adoty_p, this->adotz_p);
	}
	// Calculate corrected values.
	for (int i=0; i<this->nbs->N; i++) {
		this->iterate_correction_ind(i);
	}

	// Store the new positions and velocities in the system.
	for (int i=0; i<this->nbs->N; i++) {
		this->nbs->rx[i] = this->rx_next[i];
		this->nbs->ry[i] = this->ry_next[i];
		this->nbs->rz[i] = this->rz_next[i];
		this->nbs->vx[i] = this->vx_next[i];
		this->nbs->vy[i] = this->vy_next[i];
		this->nbs->vz[i] = this->vz_next[i];
	}
}


HermiteIterated::HermiteIterated(std::string filename, double dt, int iteration_depth)
	: Integrator(filename) {
	// Create new Integrator.
	this->name = "HermiteIterated";
	this->iteration_depth = iteration_depth;
	this->dt = dt;
	this->rx_p = new double[this->nbs->N];
	this->ry_p = new double[this->nbs->N];
	this->rz_p = new double[this->nbs->N];
	this->vx_p = new double[this->nbs->N];
	this->vy_p = new double[this->nbs->N];
	this->vz_p = new double[this->nbs->N];
	this->ax_p = new double[this->nbs->N];
	this->ay_p = new double[this->nbs->N];
	this->az_p = new double[this->nbs->N];
	this->adotx_p = new double[this->nbs->N];
	this->adoty_p = new double[this->nbs->N];
	this->adotz_p = new double[this->nbs->N];
	this->a2x = new double[this->nbs->N];
	this->a2y = new double[this->nbs->N];
	this->a2z = new double[this->nbs->N];
	this->a3x = new double[this->nbs->N];
	this->a3y = new double[this->nbs->N];
	this->a3z = new double[this->nbs->N];

	for (int i=0; i<this->nbs->N; i++) {
		this->rx_p[i] = 0;
		this->ry_p[i] = 0;
		this->rz_p[i] = 0;
		this->vx_p[i] = 0;
		this->vy_p[i] = 0;
		this->vz_p[i] = 0;
		this->ax_p[i] = 0;
		this->ay_p[i] = 0;
		this->az_p[i] = 0;
		this->adotx_p[i] = 0;
		this->adoty_p[i] = 0;
		this->adotz_p[i] = 0;
		this->a2x[i] = 0;
		this->a2y[i] = 0;
		this->a2z[i] = 0;
		this->a3x[i] = 0;
		this->a3y[i] = 0;
		this->a3z[i] = 0;
	}
}

HermiteIterated::~HermiteIterated() {
	 delete[] rx_p;
	 delete[] ry_p;
	 delete[] rz_p;
	 delete[] vx_p;
	 delete[] vy_p;
	 delete[] vz_p;
	 delete[] ax_p;
	 delete[] ay_p;
	 delete[] az_p;
	 delete[] adotx_p;
	 delete[] adoty_p;
	 delete[] adotz_p;
	 delete[] a2x;
	 delete[] a2y;
	 delete[] a2z;
	 delete[] a3x;
	 delete[] a3y;
	 delete[] a3z;
}

void HermiteIterated::iterate_ind(int i) {
	// Calculate predicted velocities.
	this->vx_p[i] = this->nbs->vx[i] + this->nbs->ax[i]*this->dt + 0.5*this->nbs->adotx[i]*this->dt*this->dt;
	this->vy_p[i] = this->nbs->vy[i] + this->nbs->ay[i]*this->dt + 0.5*this->nbs->adoty[i]*this->dt*this->dt;
	this->vz_p[i] = this->nbs->vz[i] + this->nbs->az[i]*this->dt + 0.5*this->nbs->adotz[i]*this->dt*this->dt;
	// Calculate predicted positions.
	this->rx_p[i] = this->nbs->rx[i] + this->nbs->vx[i]*this->dt + 0.5*this->nbs->ax[i]*this->dt*this->dt + this->nbs->adotx[i]*this->dt*this->dt*this->dt/6;
	this->ry_p[i] = this->nbs->ry[i] + this->nbs->vy[i]*this->dt + 0.5*this->nbs->ay[i]*this->dt*this->dt + this->nbs->adoty[i]*this->dt*this->dt*this->dt/6;
	this->rz_p[i] = this->nbs->rz[i] + this->nbs->vz[i]*this->dt + 0.5*this->nbs->az[i]*this->dt*this->dt + this->nbs->adotz[i]*this->dt*this->dt*this->dt/6;
}

void HermiteIterated::iterate_correction_ind(int i) {
	// Calculate the correction values. Store corrected values in _next variables.

	// Calculate higher derivatives of a.
	this->calc_acc_derivatives(i);

	this->vx_next[i] = this->nbs->vx[i] + 0.5*(this->ax_p[i] + this->nbs->ax[i])*this->dt + (this->adotx_p[i] - this->nbs->adotx[i])/this->dt*this->dt/12;
	this->vy_next[i] = this->nbs->vy[i] + 0.5*(this->ay_p[i] + this->nbs->ay[i])*this->dt + (this->adoty_p[i] - this->nbs->adoty[i])/this->dt*this->dt/12;
	this->vz_next[i] = this->nbs->vz[i] + 0.5*(this->az_p[i] + this->nbs->az[i])*this->dt + (this->adotz_p[i] - this->nbs->adotz[i])/this->dt*this->dt/12;

	this->rx_next[i] = this->nbs->rx[i] + 0.5*(this->vx_next[i] + this->nbs->vx[i])*this->dt + (this->ax_p[i] - this->nbs->ax[i])*this->dt*this->dt/12;
	this->ry_next[i] = this->nbs->ry[i] + 0.5*(this->vy_next[i] + this->nbs->vy[i])*this->dt + (this->ay_p[i] - this->nbs->ay[i])*this->dt*this->dt/12;
	this->rz_next[i] = this->nbs->rz[i] + 0.5*(this->vz_next[i] + this->nbs->vz[i])*this->dt + (this->az_p[i] - this->nbs->az[i])*this->dt*this->dt/12;

}

void HermiteIterated::calc_acc_derivatives(int i) {
	// Calculate higher derivatives of a.
	this->a2x[i] = -6*(this->nbs->ax[i]-this->ax_p[i])/this->dt/this->dt - 2*(2*this->nbs->adotx[i]+this->adotx_p[i])/this->dt;
	this->a2y[i] = -6*(this->nbs->ay[i]-this->ay_p[i])/this->dt/this->dt - 2*(2*this->nbs->adoty[i]+this->adoty_p[i])/this->dt;
	this->a2z[i] = -6*(this->nbs->az[i]-this->az_p[i])/this->dt/this->dt - 2*(2*this->nbs->adotz[i]+this->adotz_p[i])/this->dt;

	this->a3x[i] = 12*(this->nbs->ax[i]-this->ax_p[i])/this->dt/this->dt/this->dt + 6*(this->nbs->adotx[i]+this->adotx_p[i])/this->dt/this->dt;
	this->a3y[i] = 12*(this->nbs->ay[i]-this->ay_p[i])/this->dt/this->dt/this->dt + 6*(this->nbs->adoty[i]+this->adoty_p[i])/this->dt/this->dt;
	this->a3z[i] = 12*(this->nbs->az[i]-this->az_p[i])/this->dt/this->dt/this->dt + 6*(this->nbs->adotz[i]+this->adotz_p[i])/this->dt/this->dt;
}


void HermiteIterated::start() {} // Nothing to do.

void HermiteIterated::stop() {} // Nothing to do.

void HermiteIterated::iterate() {

	// Calculate accelerations.
	this->nbs->calc_acc();
	this->nbs->calc_dot_acc();

	// Iterate all particles.
	for (int i=0; i<this->nbs->N; i++) {
		this->iterate_ind(i);
	}
	// Calculate new accelerations and time derivatives of them using predicted positions and velocities.
	for (int i=0; i<this->nbs->N; i++) {
		this->nbs->acc(i, this->rx_p, this->ry_p, this->rz_p, this->ax_p, this->ay_p, this->az_p);
		this->nbs->dot_acc(i, this->rx_p, this->ry_p, this->rz_p, this->vx_p, this->vy_p, this->vz_p, this->adotx_p, this->adoty_p, this->adotz_p);
	}
	// Calculate corrected values.
	for (int i=0; i<this->nbs->N; i++) {
		this->iterate_correction_ind(i);
	}

	// Iterate the Hermite step according to iteration_depth.
	for (int n=0; n<this->iteration_depth; n++) {
		// Calculate new accelerations and time derivatives of them using corrected positions and velocities.
		for (int i=0; i<this->nbs->N; i++) {
			this->nbs->acc(i, this->rx_next, this->ry_next, this->rz_next, this->ax_p, this->ay_p, this->az_p);
			this->nbs->dot_acc(i, this->rx_next, this->ry_next, this->rz_next, this->vx_next, this->vy_next, this->vz_next, this->adotx_p, this->adoty_p, this->adotz_p);
		}
		// Calculate new corrected values.
		for (int i=0; i<this->nbs->N; i++) {
			this->iterate_correction_ind(i);
		}

	}


	// Store the new positions and velocities in the system.
	for (int i=0; i<this->nbs->N; i++) {
		this->nbs->rx[i] = this->rx_next[i];
		this->nbs->ry[i] = this->ry_next[i];
		this->nbs->rz[i] = this->rz_next[i];
		this->nbs->vx[i] = this->vx_next[i];
		this->nbs->vy[i] = this->vy_next[i];
		this->nbs->vz[i] = this->vz_next[i];
	}
}



RungeKutta::RungeKutta(std::string filename, double dt)
	: Integrator(filename) {
	// Create new Integrator.
	this->name = "RungeKutta";
	this->dt = dt;
	this->rx_t = new double[this->nbs->N];
	this->ry_t = new double[this->nbs->N];
	this->rz_t = new double[this->nbs->N];
	this->vx_t = new double[this->nbs->N];
	this->vy_t = new double[this->nbs->N];
	this->vz_t = new double[this->nbs->N];
	this->ax_t = new double[this->nbs->N];
	this->ay_t = new double[this->nbs->N];
	this->az_t = new double[this->nbs->N];
	this->rx_aux = new double[this->nbs->N];
	this->ry_aux = new double[this->nbs->N];
	this->rz_aux = new double[this->nbs->N];
	for (int i=0; i<this->nbs->N; i++) {
		this->rx_t[i] = 0;
		this->ry_t[i] = 0;
		this->rz_t[i] = 0;
		this->vx_t[i] = 0;
		this->vy_t[i] = 0;
		this->vz_t[i] = 0;
		this->ax_t[i] = 0;
		this->ay_t[i] = 0;
		this->az_t[i] = 0;
		this->rx_aux[i] = 0;
		this->ry_aux[i] = 0;
		this->rz_aux[i] = 0;
	}
}

RungeKutta::~RungeKutta() {
	 delete[] rx_t;
	 delete[] ry_t;
	 delete[] rz_t;
	 delete[] vx_t;
	 delete[] vy_t;
	 delete[] vz_t;
	 delete[] ax_t;
	 delete[] ay_t;
	 delete[] az_t;
}

void RungeKutta::start() {}	//	Nothing to do.
void RungeKutta::stop() {}	//	Nothing to do.

void RungeKutta::iterate_ind(int i) {
	// Iterate velocity.
	this->vx_next[i] = this->nbs->vx[i] + this->nbs->ax[i]*this->dt;

	// Iterate position.
	this->rx_next[i] = this->nbs->rx[i] + this->nbs->vx[i]*this->dt;

}

void RungeKutta::iterate() {

	// 1st step
	this->nbs->calc_acc();
	for (int i=0; i<this->nbs->N; i++) {
		this->vx_t[i] = this->nbs->ax[i]*this->dt;
		this->vy_t[i] = this->nbs->ay[i]*this->dt;
		this->vz_t[i] = this->nbs->az[i]*this->dt;

		// Update v_next.
		this->vx_next[i] = this->nbs->vx[i] + this->vx_t[i]/6;
		this->vy_next[i] = this->nbs->vy[i] + this->vy_t[i]/6;
		this->vz_next[i] = this->nbs->vz[i] + this->vz_t[i]/6;

		this->rx_t[i] = this->nbs->vx[i]*this->dt;
		this->ry_t[i] = this->nbs->vy[i]*this->dt;
		this->rz_t[i] = this->nbs->vz[i]*this->dt;

		// Update r-next.
		this->rx_next[i] = this->nbs->rx[i] + this->rx_t[i]/6;
		this->ry_next[i] = this->nbs->ry[i] + this->ry_t[i]/6;
		this->rz_next[i] = this->nbs->rz[i] + this->rz_t[i]/6;
	}

	// 2nd step
	// Calculate acceleration.
	for (int i=0; i<this->nbs->N; i++) {
		this->rx_aux[i] = this->nbs->rx[i] + 0.5*this->rx_t[i];
		this->ry_aux[i] = this->nbs->ry[i] + 0.5*this->ry_t[i];
		this->rz_aux[i] = this->nbs->rz[i] + 0.5*this->rz_t[i];
	}
	for (int i=0; i<this->nbs->N; i++) {
		this->nbs->acc(i, this->rx_aux, this->ry_aux, this->rz_aux, this->ax_t, this->ay_t, this->az_t);
	}
	// Calculate new velcoties and positions.
	for (int i=0; i<this->nbs->N; i++) {
		this->rx_t[i] = (this->nbs->vx[i] + 0.5*this->vx_t[i])*this->dt;
		this->ry_t[i] = (this->nbs->vy[i] + 0.5*this->vy_t[i])*this->dt;
		this->rz_t[i] = (this->nbs->vz[i] + 0.5*this->vz_t[i])*this->dt;

		// Update r-next.
		this->rx_next[i] += this->rx_t[i]/3;
		this->ry_next[i] += this->ry_t[i]/3;
		this->rz_next[i] += this->rz_t[i]/3;

		this->vx_t[i] = this->ax_t[i]*this->dt;
		this->vy_t[i] = this->ay_t[i]*this->dt;
		this->vz_t[i] = this->az_t[i]*this->dt;

		// Update v_next.
		this->vx_next[i] += this->vx_t[i]/3;
		this->vy_next[i] += this->vy_t[i]/3;
		this->vz_next[i] += this->vz_t[i]/3;
	}

	// 3nd step (same as 2nd step)
		// Calculate acceleration.
		for (int i=0; i<this->nbs->N; i++) {
			this->rx_aux[i] = this->nbs->rx[i] + 0.5*this->rx_t[i];
			this->ry_aux[i] = this->nbs->ry[i] + 0.5*this->ry_t[i];
			this->rz_aux[i] = this->nbs->rz[i] + 0.5*this->rz_t[i];
		}
		for (int i=0; i<this->nbs->N; i++) {
			this->nbs->acc(i, this->rx_aux, this->ry_aux, this->rz_aux, this->ax_t, this->ay_t, this->az_t);
		}
		// Calculate new velcoties and positions.
		for (int i=0; i<this->nbs->N; i++) {
			this->rx_t[i] = (this->nbs->vx[i] + 0.5*this->vx_t[i])*this->dt;
			this->ry_t[i] = (this->nbs->vy[i] + 0.5*this->vy_t[i])*this->dt;
			this->rz_t[i] = (this->nbs->vz[i] + 0.5*this->vz_t[i])*this->dt;

			// Update r-next.
			this->rx_next[i] += this->rx_t[i]/3;
			this->ry_next[i] += this->ry_t[i]/3;
			this->rz_next[i] += this->rz_t[i]/3;

			this->vx_t[i] = this->ax_t[i]*this->dt;
			this->vy_t[i] = this->ay_t[i]*this->dt;
			this->vz_t[i] = this->az_t[i]*this->dt;

			// Update v_next.
			this->vx_next[i] += this->vx_t[i]/3;
			this->vy_next[i] += this->vy_t[i]/3;
			this->vz_next[i] += this->vz_t[i]/3;
		}

		// 4nd step
		// Calculate acceleration.
		for (int i=0; i<this->nbs->N; i++) {
			this->rx_aux[i] = this->nbs->rx[i] + this->rx_t[i];
			this->ry_aux[i] = this->nbs->ry[i] + this->ry_t[i];
			this->rz_aux[i] = this->nbs->rz[i] + this->rz_t[i];
		}
		for (int i=0; i<this->nbs->N; i++) {
			this->nbs->acc(i, this->rx_aux, this->ry_aux, this->rz_aux, this->ax_t, this->ay_t, this->az_t);
		}
		// Calculate new velcoties and positions.
		for (int i=0; i<this->nbs->N; i++) {
			this->rx_t[i] = (this->nbs->vx[i] + this->vx_t[i])*this->dt;
			this->ry_t[i] = (this->nbs->vy[i] + this->vy_t[i])*this->dt;
			this->rz_t[i] = (this->nbs->vz[i] + this->vz_t[i])*this->dt;

			// Update r-next.
			this->rx_next[i] += this->rx_t[i]/6;
			this->ry_next[i] += this->ry_t[i]/6;
			this->rz_next[i] += this->rz_t[i]/6;

			this->vx_t[i] = this->ax_t[i]*this->dt;
			this->vy_t[i] = this->ay_t[i]*this->dt;
			this->vz_t[i] = this->az_t[i]*this->dt;

			// Update v_next.
			this->vx_next[i] += this->vx_t[i]/6;
			this->vy_next[i] += this->vy_t[i]/6;
			this->vz_next[i] += this->vz_t[i]/6;
		}



	// Store the new positions and velocities in the system.
	for (int i=0; i<this->nbs->N; i++) {
		this->nbs->rx[i] = this->rx_next[i];
		this->nbs->ry[i] = this->ry_next[i];
		this->nbs->rz[i] = this->rz_next[i];
		this->nbs->vx[i] = this->vx_next[i];
		this->nbs->vy[i] = this->vy_next[i];
		this->nbs->vz[i] = this->vz_next[i];
	}
}
