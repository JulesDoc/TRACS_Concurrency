/*
 * @ Copyright 2014-2017 CERN and Instituto de Fisica de Cantabria - Universidad de Cantabria. All rigths not expressly granted are reserved [tracs.ssd@cern.ch]
 * This file is part of TRACS.
 *
 * TRACS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the Licence.
 *
 * TRACS is distributed in the hope that it will be useful , but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with TRACS. If not, see <http://www.gnu.org/licenses/>
 */

/************************************DoTRACSOnly***********************************
 *
 * Basic main TRACS execution which gives back the Neff of a defined detector and some extra information related to diffusion, if swithed ON.
 *
 */

#include <thread>
#include <boost/asio.hpp>
#include <TRACSFit.h>
#include <TRACSInterface.h>
#include <stdio.h>
#include "../include/Threading.h"
#include <Utilities.h>

std::vector<TRACSInterface*> TRACSsim;
std::vector<std::thread> t;
std::vector<double> vector_zValues, vector_voltValues, vector_yValues;
std::string scanType;
double dTime;
double max_time;

void spread_into_threads();

int main( int argc, char *argv[]) {

	std::string carrierFile;
	std::vector<std::string> carrierThread_fileNames;
	std::valarray<std::valarray<double>> i_total;
	TH1D *i_rc;


	int counted_numLines = 0;
	int number_of_lines = 0;
	int nns, count2;
	std::string line;
	std::string temp;


	if(argc==1){
		num_threads = std::thread::hardware_concurrency(); // No. of threads = No. of cores

	}

	if(argc==2){
		num_threads = atoi(argv[1]);
		if (num_threads == 0){
			num_threads = 1;
		}
	}

	if(argc==3){
		num_threads = atoi(argv[1]);
		fnm = argv[2];

	}

	//file creation
	for (int i = 0; i < num_threads; ++i) {
		temp = std::to_string(i);
		carrierThread_fileNames.push_back("file" + temp);
	}

	utilities::parse_config_file(fnm, carrierFile);
	std::ifstream in(carrierFile);
	while (in){

		for (int i = 0 ; i < num_threads ; i++){
			std::ofstream outfile;
			const char * c = carrierThread_fileNames[i].c_str();
			outfile.open(c, std::ofstream::out | std::ofstream::app);
			if (in) std::getline(in, line);
			outfile << line << std::endl;
			outfile.close();
		}

	}
	in.close();




	spread_into_threads();
	double timeSteps = (int) std::floor(max_time / dTime);
	int total_sizeZ = vector_voltValues.size() * vector_zValues.size();
	int total_sizeY = vector_voltValues.size() * vector_yValues.size();


	if (scanType == "edge"){
		if (vector_yValues.size() > 1){
					std::cout << "This execution is wrongly configured with edge-TCT; Check parameters." << std::endl;
					std::quick_exit(1);
				}
		//i_rc_array.resize(vector_voltValues.size());
		i_rc_array.resize(total_sizeZ);
		vItotals.resize(total_sizeZ);
		for (int i = 0; i < total_sizeZ ; i++)
			vItotals[i].resize(timeSteps);
			//i_rc_array[i].resize(vector_zValues.size());
	}


	if (scanType == "top" || scanType == "bottom"){
		if (vector_zValues.size() > 1){
					std::cout << "This execution is wrongly configured with top/bottom-TCT; Check parameters." << std::endl;
					std::quick_exit(1);
				}
		//i_rc_array.resize(vector_voltValues.size());
		i_rc_array.resize(total_sizeY);
		vItotals.resize(total_sizeY);
		for (int i = 0; i < total_sizeY ; i++)
			//i_rc_array[i].resize(vector_yValues.size());
			vItotals[i].resize(timeSteps);

	}

	for (int i = 0 ; i < vItotals.size() ; i++){
		for (int j = 0 ; j < vItotals[i].size() ; j++)
			vItotals[i][j] = 0;
	}

	/*for (int i = 0 ; i < vItotals.size() ; i++){
						for (int j = 0 ; j < vItotals[i].size() ; j++)
							std::cout << "i " << i << "; j " << j << "    " <<  vItotals[i][j] << std::endl;
					}*/


	TRACSsim.resize(num_threads);
	t.resize(num_threads);
	for (int i = 0; i < num_threads; ++i) {
		t[i] = std::thread(call_from_thread, i, std::ref(carrierThread_fileNames[i]), vector_zValues, vector_yValues, vector_voltValues);
	}
	for (int i = 0; i < num_threads; ++i) {
		t[i].join();
	}


	//Current to rc array of TH1D -> root file
	if (scanType == "edge"){
		for (int i = 0 ; i < i_rc_array.size(); i++ ){

			//for (int j = 0 ; j <= vector_zValues.size(); j++ ){
				TString htit, hname;
				htit.Form("ramo_rc%d%d", 0, count2);
				hname.Form("Ramo_current_%d_%d", 0, count2);
				i_rc = new TH1D(htit,hname, timeSteps, 0.0, max_time);
				for (int k = 1 ; k < timeSteps; k++ ){

					i_rc->SetBinContent(k+1, vItotals[i][k]);

				}
				i_rc_array[i] = i_rc;
				i_rc = nullptr;
				count2++;
			//}

		}

	}


	//Current to rc array of TH1D -> root file
	if (scanType == "top" || scanType == "bottom"){
		for (int i = 0 ; i < i_rc_array.size(); i++ ){

			//for (int j = 0 ; j <= vector_yValues.size(); j++ ){
				TString htit, hname;
				htit.Form("ramo_rc%d%d", 0, count2);
				hname.Form("Ramo_current_%d_%d", 0, count2);
				i_rc = new TH1D(htit,hname, timeSteps, 0.0, max_time);
				for (int k = 1 ; k < timeSteps; k++ ){

					i_rc->SetBinContent(k+1, vItotals[i][k]);
				}
				i_rc_array[i] = i_rc;
				i_rc = nullptr;
				count2++;
			//}

		}

	}

	//write output to single file!
	TRACSsim[0]->write_to_file(0);
	for (int i = 0 ; i < num_threads ; i++){
			const char * c = carrierThread_fileNames[i].c_str();
			remove(c);
		}

	for (uint i = 0; i < TRACSsim.size(); i++)	{
		delete TRACSsim[i];
	}

	std::quick_exit(1);


}


//-----------


Double_t TRACSFit::operator() ( const std::vector<Double_t>& par  ) const {

	return 0;
}


void spread_into_threads(){

	double vInit, deltaV, vMax, v_depletion, zInit, zMax, deltaZ, yInit, yMax, deltaY;
	int zSteps, ySteps, vSteps, seeker;
	bool lastElement = false;


	uint index = 1;

	//Reading file to get values that will determine scan vectors
	utilities::parse_config_file(fnm, scanType, vInit, deltaV, vMax, v_depletion, zInit, zMax, deltaZ, yInit, yMax, deltaY, dTime, max_time);

	//Reserve to 300 since is a common measurement for a detector Z dimension.
	//Reserve does not implies anything, better for performance.
	//For sure there is at least one value of Z coordinates in an edge-TCT scan
	vector_zValues.reserve(300);
	vector_zValues.push_back(zInit);
	seeker = zInit + deltaZ;

	//Filling Z scan vector
	while (seeker < zMax){
		vector_zValues.push_back(seeker);
		++index;
		seeker = seeker + deltaZ;
		if (seeker > zMax){
			vector_zValues.push_back(zMax);
			lastElement = true;
		}

	}
	if (seeker >= zMax && !lastElement && zMax !=zInit) vector_zValues.push_back(zMax);
	lastElement = false;

	//Reserve to 100 for differente voltages.
	//There will be at least one bias voltage
	vector_voltValues.reserve(300);
	vector_voltValues.push_back(vInit);
	seeker = vInit + deltaV;

	//Filling voltage scan vector
	while (seeker < vMax){
		vector_voltValues.push_back(seeker);
		++index;
		seeker = seeker + deltaV;
		if (seeker > vMax) {
			vector_voltValues.push_back(vMax);
			lastElement = true;
		}
	}
	if (seeker >= vMax && !lastElement && vMax != vInit) vector_voltValues.push_back(vMax);
	lastElement = false;

	//Reserve to 400 since is a common measurement for a detector 4 dimension.
	//For sure there is at least one value of Z coordinates in an edge-TCT scan
	vector_yValues.reserve(400);
	vector_yValues.push_back(yInit);
	seeker = yInit + deltaY;

	//Filling Y scan vector
	while (seeker < yMax){
		vector_yValues.push_back(seeker);
		++index;
		seeker = seeker + deltaY;
		if (seeker > yMax) {
			vector_yValues.push_back(yMax);
			lastElement = true;
		}
	}
	if (seeker >= yMax && !lastElement && yMax != yInit) vector_yValues.push_back(yMax);
	lastElement = false;

	//Checking for a right user configuration
	if ((vector_yValues.size() > 1) && (vector_zValues.size() > 1)){
		std::cout << "TRACS does not allow this configuration, please, choose zScan OR yScan" << std::endl;
		std::quick_exit(1);
	}



}

