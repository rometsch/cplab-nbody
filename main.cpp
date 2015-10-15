/*
 * main.cpp
 *
 *  Created on: Oct 12, 2015
 *      Author: thomas
 */

#include "Integrators.h"
#include "Logger.h"
#include <iostream>


int main() {
	Euler *integrator = new Euler("/Users/thomas/GoogleDrive/uni/ComPhy/n-body/data/in2", 0.001);

	Logger log(integrator);
	log.iterate(30000,20);
	log.plot_trajectories();

	delete integrator;

}
