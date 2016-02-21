/*
 * main.cpp
 *
 *  Created on: Oct 12, 2015
 *      Author: thomas
 */

#include "Integrators.h"
#include "Logger.h"
#include <iostream>
#include <fstream>
#include <chrono>

void aufgabe2a();
void aufgabe2b();
void aufgabe3();
void save_two_body_sim_to_file(Integrator *integrator, double dt, int num_of_tmax, int num_data_pts, std::string outfile_path, std::string outfile_extension);
void save_N_body_sim_to_file(Integrator *integrator, double dt, double num_of_tmax, int num_data_pts, std::string outfile_path, std::string outfile_extension);



int main() {

	aufgabe2a();
	aufgabe2b();
	aufgabe3();

}

void aufgabe3() {

	double eta[6] = {0.1, 0.05, 0.01, 0.005, 0.001, 0.0005};
	double num_of_tmax = 1;

	int num_data_pts = 1000;

	std::string datafile_100 = "./data/pl.100";
	std::string datafile_1k = "./data/pl.1k";

	for (int i=0; i<6	; i++) {
//		save_N_body_sim_to_file(new Euler(datafile_100, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_100_"+std::to_string(eta[i])+".txt");
//		save_N_body_sim_to_file(new Euler(datafile_1k, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_1000_"+std::to_string(eta[i])+".txt");
//		save_N_body_sim_to_file(new EulerComer(datafile_100, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_100_"+std::to_string(eta[i])+".txt");
//		save_N_body_sim_to_file(new EulerComer(datafile_1k, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_1000_"+std::to_string(eta[i])+".txt");
//		save_N_body_sim_to_file(new HermiteIterated(datafile_100, 1.0, 2), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_100_"+std::to_string(eta[i])+".txt");
//		save_N_body_sim_to_file(new HermiteIterated(datafile_1k, 1.0, 2), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_1000_"+std::to_string(eta[i])+".txt");
		save_N_body_sim_to_file(new HermitePC(datafile_100, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_100_"+std::to_string(eta[i])+".txt");
		save_N_body_sim_to_file(new HermitePC(datafile_1k, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_1000_"+std::to_string(eta[i])+".txt");
		save_N_body_sim_to_file(new LeapFrog(datafile_100, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_100_"+std::to_string(eta[i])+".txt");
		save_N_body_sim_to_file(new LeapFrog(datafile_1k, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_1000_"+std::to_string(eta[i])+".txt");
//		save_N_body_sim_to_file(new Mittelung(datafile_100, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_100_"+std::to_string(eta[i])+".txt");
//		save_N_body_sim_to_file(new Mittelung(datafile_1k, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_1000_"+std::to_string(eta[i])+".txt");
//		save_N_body_sim_to_file(new RungeKutta(datafile_100, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_100_"+std::to_string(eta[i])+".txt");
//		save_N_body_sim_to_file(new RungeKutta(datafile_1k, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_1000_"+std::to_string(eta[i])+".txt");
//		save_N_body_sim_to_file(new Verlet(datafile_100, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_100_"+std::to_string(eta[i])+".txt");
//		save_N_body_sim_to_file(new Verlet(datafile_1k, 1.0), eta[i], num_of_tmax, num_data_pts, "aufgabe3/", "_1000_"+std::to_string(eta[i])+".txt");
	}
}

void aufgabe2b() {
	int num_of_tmax = 2;
	int num_data_pts = 20;
	double eta = 0.01;
	std::vector<std::string> filenames;
	filenames.push_back("./data/in2");
	filenames.push_back("./data/in2i");
	filenames.push_back("./data/in2ii");

	// Simulation für alle drei Datenfiles und alle Integratoren durchführen.
	for (int i=0; i<3; i++) {
		save_two_body_sim_to_file(new Euler(filenames[i], 1.0), eta, num_of_tmax, num_data_pts, "aufgabe2b/", "_"+std::to_string(i)+".txt");
		save_two_body_sim_to_file(new EulerComer(filenames[i], 1.0), eta, num_of_tmax, num_data_pts, "aufgabe2b/", "_"+std::to_string(i)+".txt");
		save_two_body_sim_to_file(new HermiteIterated(filenames[i], 1.0, 2), eta, num_of_tmax, num_data_pts, "aufgabe2b/", "_"+std::to_string(i)+".txt");
		save_two_body_sim_to_file(new HermitePC(filenames[i], 1.0), eta, num_of_tmax, num_data_pts, "aufgabe2b/", "_"+std::to_string(i)+".txt");
		save_two_body_sim_to_file(new LeapFrog(filenames[i], 1.0), eta, num_of_tmax, num_data_pts, "aufgabe2b/", "_"+std::to_string(i)+".txt");
		save_two_body_sim_to_file(new Mittelung(filenames[i], 1.0), eta, num_of_tmax, num_data_pts, "aufgabe2b/", "_"+std::to_string(i)+".txt");
		save_two_body_sim_to_file(new RungeKutta(filenames[i], 1.0), eta, num_of_tmax, num_data_pts, "aufgabe2b/", "_"+std::to_string(i)+".txt");
		save_two_body_sim_to_file(new Verlet(filenames[i], 1.0), eta, num_of_tmax, num_data_pts, "aufgabe2b/", "_"+std::to_string(i)+".txt");
	}

}


void aufgabe2a() {
	double eta[5] = {0.5, 0.1, 0.05, 0.01, 0.005};
	int num_eta = 5;
	int num_data_pts = 20;
	std::string datafile_path = "./data/in2";

	// Für alle Integratoren und für alle eta Log in Datei speichern.
	for (int i=0; i<num_eta ; i++ ) {
		save_two_body_sim_to_file(new Euler(datafile_path, 1.0), eta[i], 10, num_data_pts, "aufgabe2a/", "_"+std::to_string(eta[i])+".txt");
		save_two_body_sim_to_file(new EulerComer(datafile_path, 1.0), eta[i], 10, num_data_pts, "aufgabe2a/", "_"+std::to_string(eta[i])+".txt");
		save_two_body_sim_to_file(new HermiteIterated(datafile_path, 1.0, 2), eta[i], 10, num_data_pts, "aufgabe2a/", "_"+std::to_string(eta[i])+".txt");
		save_two_body_sim_to_file(new HermitePC(datafile_path, 1.0), eta[i], 10, num_data_pts, "aufgabe2a/", "_"+std::to_string(eta[i])+".txt");
		save_two_body_sim_to_file(new LeapFrog(datafile_path, 1.0), eta[i], 10, num_data_pts, "aufgabe2a/", "_"+std::to_string(eta[i])+".txt");
		save_two_body_sim_to_file(new Mittelung(datafile_path, 1.0), eta[i], 10, num_data_pts, "aufgabe2a/", "_"+std::to_string(eta[i])+".txt");
		save_two_body_sim_to_file(new RungeKutta(datafile_path, 1.0), eta[i], 10, num_data_pts, "aufgabe2a/", "_"+std::to_string(eta[i])+".txt");
		save_two_body_sim_to_file(new Verlet(datafile_path, 1.0), eta[i], 10, num_data_pts, "aufgabe2a/", "_"+std::to_string(eta[i])+".txt");
	}


}

void save_two_body_sim_to_file(Integrator *integrator, double dt, int num_of_tmax, int num_data_pts, std::string outfile_path, std::string outfile_extension) {
	// Adjust integrator to new dt.
	integrator->set_dt(dt);

	double tmax = integrator->nbs->t_max;
	int step = int(num_of_tmax*tmax/dt/num_data_pts);
	if (step==0) step = 1;

	std::string integrator_name = integrator->name;
	Logger log(integrator);
	auto start = std::chrono::system_clock::now();
	log.iterate(num_of_tmax,step);
	auto end = std::chrono::system_clock::now();
	// New data file to write to.
	std::string outfile_name = outfile_path+integrator_name+outfile_extension;
	std::ofstream data_file(outfile_name, std::ios::trunc	);
	// Save ellapsed time to file
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	data_file << "#t=" << elapsed.count() << std::endl;
	// Save logged data to file.
	log.output_drift_two_body(data_file);
	data_file.close();
}

void save_N_body_sim_to_file(Integrator *integrator, double dt, double num_of_tmax, int num_data_pts, std::string outfile_path, std::string outfile_extension) {
	// Adjust integrator to new dt.
	integrator->set_dt(dt);

	double tmax = integrator->nbs->t_max;
	int number_of_datapts = 20;
	int step = int(num_of_tmax*tmax/dt/num_data_pts);
	if (step==0) step = 1;

	std::string integrator_name = integrator->name;
	Logger log(integrator);
	auto start = std::chrono::system_clock::now();
	log.iterate(num_of_tmax,step);
	auto end = std::chrono::system_clock::now();
	// New data file to write to.
	std::string outfile_name = outfile_path+integrator_name+outfile_extension;
	std::ofstream data_file(outfile_name, std::ios::trunc	);
	// Save ellapsed time to file
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	data_file << "#t=" << elapsed.count() << std::endl;
	// Save logged data to file.
	log.output_energy_drift(data_file);
	data_file.close();
}


