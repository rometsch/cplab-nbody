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
