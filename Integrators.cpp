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


Euler::Euler(std::string filename, double dt)
	: Integrator(filename) {
	// Create new Integrator.
	this->dt = dt;
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

Euler_Comer::Euler_Comer(std::string filename, double dt)
	: Integrator(filename) {
	// Create new Integrator.
	this->dt = dt;
}

Euler_Comer::~Euler_Comer() {

}

void Euler_Comer::start() {}	//	Nothing to do.
void Euler_Comer::stop() {}	//	Nothing to do.


void Euler_Comer::iterate_ind(int i) {
	// Iterate velocity.
	this->vx_next[i] = this->nbs->vx[i] + this->nbs->ax[i]*this->dt;
	this->vy_next[i] = this->nbs->vy[i] + this->nbs->ay[i]*this->dt;
	this->vz_next[i] = this->nbs->vz[i] + this->nbs->az[i]*this->dt;
	// Iterate position.
	this->rx_next[i] = this->nbs->rx[i] + this->vx_next[i]*this->dt;
	this->ry_next[i] = this->nbs->ry[i] + this->vy_next[i]*this->dt;
	this->rz_next[i] = this->nbs->rz[i] + this->vz_next[i]*this->dt;
}

void Euler_Comer::iterate() {
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
	double a2x, a2y, a2z, a3x, a3y, a3z;

	// Calculate higher derivatives of a.
	a2x = -6*(this->nbs->ax[i]-this->ax_p[i])/this->dt/this->dt - 2*(2*this->nbs->adotx[i]+this->adotx_p[i])/this->dt;
	a2y = -6*(this->nbs->ay[i]-this->ay_p[i])/this->dt/this->dt - 2*(2*this->nbs->adoty[i]+this->adoty_p[i])/this->dt;
	a2z = -6*(this->nbs->az[i]-this->az_p[i])/this->dt/this->dt - 2*(2*this->nbs->adotz[i]+this->adotz_p[i])/this->dt;

	a3x = 12*(this->nbs->ax[i]-this->ax_p[i])/this->dt/this->dt/this->dt + 6*(this->nbs->adotx[i]+this->adotx_p[i])/this->dt/this->dt;
	a3y = 12*(this->nbs->ay[i]-this->ay_p[i])/this->dt/this->dt/this->dt + 6*(this->nbs->adoty[i]+this->adoty_p[i])/this->dt/this->dt;
	a3z = 12*(this->nbs->az[i]-this->az_p[i])/this->dt/this->dt/this->dt + 6*(this->nbs->adotz[i]+this->adotz_p[i])/this->dt/this->dt;

	this->vx_next[i] = this->vx_p[i] + a2x*this->dt*this->dt*this->dt/6 + a3x*this->dt*this->dt*this->dt*this->dt/24;
	this->vy_next[i] = this->vy_p[i] + a2y*this->dt*this->dt*this->dt/6 + a3y*this->dt*this->dt*this->dt*this->dt/24;
	this->vz_next[i] = this->vz_p[i] + a2z*this->dt*this->dt*this->dt/6 + a3z*this->dt*this->dt*this->dt*this->dt/24;

	this->rx_next[i] = this->rx_p[i] + a2x*this->dt*this->dt*this->dt*this->dt/24 + a3x*this->dt*this->dt*this->dt*this->dt*this->dt/120;
	this->ry_next[i] = this->ry_p[i] + a2y*this->dt*this->dt*this->dt*this->dt/24 + a3y*this->dt*this->dt*this->dt*this->dt*this->dt/120;
	this->rz_next[i] = this->rz_p[i] + a2z*this->dt*this->dt*this->dt*this->dt/24 + a3z*this->dt*this->dt*this->dt*this->dt*this->dt/120;
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


