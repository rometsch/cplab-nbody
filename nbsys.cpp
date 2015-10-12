/*
 * nbodysys.cpp
 *
 *  Caeated on: Oct 12, 2015
 *      Author: thomas
 */

#include "nbsys.h"

nbsys::nbsys(std::string filename) {
	// Constructor for N-body system with initial data given in file.

	// Initialize system quantities.
	Rx = 0, Ry = 0, Rz = 0;			// Center of mass (COM) coordinate.
	P = 0, Px = 0, Py = 0, Pz = 0;	// Total momentum.
	Pcom = 0, Pcomx = 0, Pcomy = 0, Pcomz = 0;	// Center of mass momentum
	J = 0, Jx = 0, Jy = 0, Jz = 0;	// Total angular momentum.
	E = 0;		// Total energy.

	// Quantities for 2-body systems.
	j = 0, jx = 0, jy = 0, jz = 0;	// Angular momentum.
	e = 0, ex = 0, ey = 0, ez = 0;	// Runge-Lenz vector.
	ae = 0;							// Major axis of ellipse.


	// Load data from a sample and store it in the systems' variables.
	std::ifstream dfile(filename);
	std::string line;

	// Read system parameters from first line. Must be: " N   t_max   eta".
	dfile >> this->N;	// Number of particles.
	dfile >> this->t_max;	// Maximal integration time.
	dfile >> this->eta;		// Timestep parameter.


	// Allocate new arrays.
	// Mass.
	this->m = new double[this->N];
	// Position.
	this->rx = new double[this->N];
	this->ry = new double[this->N];
	this->rz = new double[this->N];
	// Velocity.
	this->vx = new double[this->N];
	this->vy = new double[this->N];
	this->vz = new double[this->N];
	// Acceleration.
	this->ax = new double[this->N];
	this->ay = new double[this->N];
	this->az = new double[this->N];
	for (int i=0; i<this->N; i++) {
		this->ax[i] = 0;
		this->ay[i] = 0;
		this->az[i] = 0;
	}

	double M=0;
	// Read masses from file. Stored as a single value per line.
	for (int i=0; i<this->N; i++) {
		dfile >> this->m[i];
		M += this->m[i];
	}
	// Normalize total mass to 1.
	for (int i=0; i<this->N; i++) this->m[i] /= M;


	// Read positions from file. Stored as three values per line, "x y z" components.
	for (int i=0; i<this->N; i++) {
		dfile >> this->rx[i] >> this->ry[i] >> this->rz[i];
	}

	// Read positions from file. Stored as three values per line, "x y z" components.
	for (int i=0; i<this->N; i++) {
		dfile >> this->vx[i] >> this->vy[i] >> this->vz[i];
	}

	// Close file.
	dfile.close();


	// Calculate com coordinates and com momentum (= velocity because total mass normalized to 1).
	this->calc_com();
	this->calc_momentum();
	// Store the com momentum in the appropriate variable.
	this->Pcom = this->P;
	this->Pcomx = this->Px;
	this->Pcomy = this->Py;
	this->Pcomz = this->Pz;

	// Transform the loaded coordinates to the com system.
	this->transform_to_com();

}

nbsys::~nbsys() {
	// Destructor.

	// Delete all the allocated arrays.
	// Mass.
	delete[] this->m;
	// Position.
	delete[] this->rx;
	delete[] this->ry;
	delete[] this->rz;
	// Velocity.
	delete[] this->vx;
	delete[] this->vy;
	delete[] this->vz;
	// Acceleration.
	delete[] this->ax;
	delete[] this->ay;
	delete[] this->az;
}


void nbsys::print_system() {
	// Print r,v,a for all particles.

	std::cout << "N-body system with " << this->N << " particles:" << std::endl << std::endl;

	// Print information for all particles.
	for (int i=0; i<this->N; i++) {
		std::cout << "Particle " << i << std::endl;
		std::cout << "m = " << this->m[i] << std::endl;
		std::cout << "r = (" << this->rx[i] << "," << this->ry[i] << "," << this->rz[i] << ")" << std::endl;
		std::cout << "v = (" << this->vx[i] << "," << this->vy[i] << "," << this->vz[i] << ")" << std::endl;
		std::cout << "a = (" << this->ax[i] << "," << this->ay[i] << "," << this->az[i] << ")" << std::endl;
	}

}


void nbsys::calc_energy() {
	// Calculate the total energy.
	double E = 0;
	for (int i=0; i<this->N; i++) {
		// Kinetic energy.
		E += 0.5*this->m[i]*(this->vx[i]*this->vx[i] + this->vy[i]*this->vy[i] + this->vz[i]*this->vz[i]);
		// Potential energy.
		for (int j=i+1; j<this->N; j++) {
			double dr_sq = (this->rx[i] - this->rx[j])*(this->rx[i] - this->rx[j]) + (this->ry[i] - this->ry[j])*(this->ry[i] - this->ry[j]) + (this->rz[i] - this->rz[j])*(this->rz[i] - this->rz[j]);
			E -= this->m[i]*this->m[j]/sqrt(dr_sq);
		}
	}
	// Update the system variable.
	this->E = E;
}

void nbsys::calc_com() {
	// Calculate the center of mass, its velocity.
	double X=0, Y=0, Z=0;
	for (int i=0; i<this->N; i++) {
		X += this->rx[i]*this->m[i];
		Y += this->ry[i]*this->m[i];
		Z += this->rz[i]*this->m[i];
	}
	// Update the system variable.
	this->Rx = X;
	this->Ry = Y;
	this->Rz = Z;
}

void nbsys::calc_momentum() {
	// Calculate the total momentum and the velocity of the center of mass.
	double Px=0, Py=0, Pz=0;
	for (int i=0; i<this->N; i++) {
		Px += this->vx[i]*this->m[i];
		Py += this->vy[i]*this->m[i];
		Pz += this->vz[i]*this->m[i];
	}
	// Update the system variable.
	this->Px = Px;
	this->Py = Py;
	this->Pz = Pz;

	double P_sq = this->Px*this->Px + this->Py*this->Py + this->Pz*this->Pz;
	this->P = sqrt(P_sq);
}
void nbsys::calc_angular_momentum() {
	// Calculate the total angular momentum.
	double Jx=0, Jy=0, Jz=0;
	for (int i=0; i<this->N; i++) {
		Jx += this->ry[i]*this->vz[i] - this->rz[i]*this->vy[i];
		Jy += this->rz[i]*this->vx[i] - this->rx[i]*this->vz[i];
		Jz += this->rx[i]*this->vy[i] - this->ry[i]*this->vx[i];
	}
	this->Jx = Jx;
	this->Jy = Jy;
	this->Jz = Jz;

	double J_sq = this->Jx*this->Jx + this->Jy*this->Jy + this->Jz*this->Jz;
	this->J = sqrt(J_sq);
}

void nbsys::transform_to_com() {
	// Transforms all coordinates and velocities to the center of mass frame using the com position in Rx,Ry,Rz and the velocity in Pcomx,Pcomy,Pcomz (total mass normalized to 1).
	for (int i=0; i<this->N; i++) {
		this->rx[i] -= this->Rx;
		this->ry[i] -= this->Ry;
		this->rz[i] -= this->Rz;
		this->vx[i] -= this->Pcomx;
		this->vy[i] -= this->Pcomy;
		this->vz[i] -= this->Pcomz;
	}

}

