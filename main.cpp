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
	Integrator *integrator = new HermitePC("/Users/thomas/GoogleDrive/uni/ComPhy/n-body/data/pla3", 0.01);

	Logger log(integrator);
	log.iterate(5000,50);
	log.plot_trajectories();

	delete integrator;

}
