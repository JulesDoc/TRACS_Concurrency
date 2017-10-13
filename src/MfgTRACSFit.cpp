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

/*
  Example

    DoTRACSFit MeasurementFile TRACS.conf "Vbias==200 && Tset == 20"

 */

//#include <TApplication.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnMachinePrecision.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnPlot.h>
#include <Minuit2/MinosError.h>
#include <Minuit2/FCNBase.h>
#include <Math/MinimizerOptions.h>

#include <boost/asio.hpp>

#include <TRACSFit.h>
#include <TRACSInterface.h>
#include <TString.h>
#include <stdio.h>

#include "../include/Global.h"

std::vector<TRACSInterface*> TRACSsim;
std::vector<std::thread> t;
boost::posix_time::time_duration total_timeTaken ;
TRACSFit *fit ;
std::string neffType;
bool irradiated;

std::vector<double> vector_zValues, vector_voltValues, vector_yValues;
std::vector<std::string> carrierThread_fileNames;
std::string scanType;
double dTime;
double max_time;

void spread_into_threads();

using namespace ROOT::Minuit2;

int main( int argc, char *argv[]) {


	std::string carrierFile;
	std::valarray<std::valarray<double>> i_total;
	TH1D *i_rc;
	int counted_numLines = 0;
	int number_of_lines = 0;
	int nns, count2;
	std::string line;
	std::string temp;

	double fitParamVdep;
	double fitParamNorm;
	double fitParamDepth;
	double fitParamCapac;
	double vBias;
	double vDep;
	double fluence;
	irradiated = false;
	vector<Double_t> parIni;
	Int_t parIniSize;
	vector<Double_t> parErr;

	//TApplication theApp("DoTRACSFit", 0, 0);

	//Number of threads
	num_threads = atoi(argv[1]);

	//Measurement file
	TString FileMeas = TString( argv[2] ) ;

	//Configuration file
	TString FileConf = TString( argv[3] ) ;
	std::string lfnm(argv[3]) ;
	fnm = lfnm;

	//Restrictions for fits
	TString how="";
	if (argc>2) how = TString( argv[4] ) ;

	//file creation
	for (int i = 0; i < num_threads; ++i) {
		temp = std::to_string(i);
		carrierThread_fileNames.push_back("file" + temp);
	}

	utilities::parse_config_file(fnm, carrierFile);

	std::ifstream infile(carrierFile);
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
	/*while (std::getline(infile, line))
	{
		++number_of_lines;
	}
	int resto_carriers = number_of_lines % num_threads;
	int carriers_per_thr = number_of_lines / num_threads;
	infile.close();

	counted_numLines = counted_numLines - resto_carriers;
	std::ifstream in(carrierFile);
	for (int i = 0; i < num_threads; ++i) {
		std::ofstream out(carrierThread_fileNames[i]);
		while (in){
			counted_numLines++;
			std::getline(in, line);
			out << line << std::endl;
			if (counted_numLines == carriers_per_thr){
				out.close();
				counted_numLines = 0;
				break;
			}

		}
	}
	in.close();*/



	spread_into_threads();
	double timeSteps = (int) std::floor(max_time / dTime);
	int total_sizeZ = vector_voltValues.size() * vector_zValues.size();
	int total_sizeY = vector_voltValues.size() * vector_yValues.size();

	if (scanType == "edge"){
		if (vector_yValues.size() > 1){
			std::cout << "This execution is wrongly configured with edge-TCT; Check parameters." << std::endl;
			std::quick_exit(1);
		}
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


	TRACSsim.resize(num_threads);
	t.resize(num_threads);
	for (int i = 0; i < num_threads; ++i) {
		t[i] = std::thread(call_from_thread, i, std::ref(carrierThread_fileNames[i]), vector_zValues, vector_yValues, vector_voltValues);
	}
	for (int i = 0; i < num_threads; ++i) {
		t[i].join();
	}


	fit = new TRACSFit( FileMeas, FileConf , how ) ;

	neffType = TRACSsim[0]->get_neff_type();
	vBias = TRACSsim[0]->get_vBias();
	vDep = TRACSsim[0]->get_vDep();
	fluence = TRACSsim[0]->get_fluence();

	/*********Begin Fit. For irradiated or undepleted dectectors****************/
	/***********************************************************************/
	//Fitting Neff and normalizator

	if ((fluence > 0) || (vBias < vDep)){

		irradiated = true;

		//Define parameters and their errors to Minuit

		parIni = TRACSsim[0]->get_NeffParam();
		parIni.push_back(TRACSsim[0]->get_fitNorm());
		parIni.push_back(TRACSsim[0]->get_depth());
		parIniSize = parIni.size() ;
		//parErr = parErr(nNeff, 60.) ;
		//To the vector initialization correctly using a preallocate variable
		//that can be used everywhere in main.
		parErr.resize(parIniSize);
		for (size_t i=0; i<parIniSize; i++)
			parErr[i] = 60.;

		//Pass parameters to Minuit

		MnUserParameters upar(parIni,parErr) ;
		for ( int i=0 ; i < parIniSize ; i++ ) {
			char pname[parIniSize]; sprintf( pname , "p%d" , i);
			upar.SetName( i , pname );
		}

		//Fix parameters
		upar.Fix(0) ;
		upar.Fix(1) ; upar.Fix(2) ;
		//upar.Fix(3) ;
		upar.Fix(4) ; upar.Fix(5); upar.Fix(6) ; upar.Fix(7);
		//upar.Fix(8); //Normalizator
		//upar.Fix(9); //Depth
		std::cout << "=============================================" << std::endl;
		std::cout << "Initial parameters: "<<upar<<std::endl;
		std::cout << "=============================================" << std::endl;
		std::cout << "=============================================" << std::endl;
		std::cout << "tolerance= " << TRACSsim[0]->GetTolerance()    << std::endl;
		std::cout << "chiFinal= " << TRACSsim[0]->GetchiFinal()      << std::endl;
		std::cout << "=============================================" << std::endl;
		std::cout << "=============================================" << std::endl;


		//Do the minimization

		ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);
		ROOT::Math::MinimizerOptions::SetDefaultTolerance(TRACSsim[0]->GetTolerance());
		MnMigrad mn( *fit , upar , MnStrategy(1)) ;
		FunctionMinimum min = mn() ;

		//Status report
		std::cout << "Total time: " << total_timeTaken.total_seconds() << std::endl ;
		std::cout << "MINIMIZATION OUTCOME: " <<  min  << std::endl ;

		//Release parameter 3, fix 0, minimize again
		upar.SetValue(3, min.UserState().Value(3) ) ; upar.SetError(3, min.UserState().Error(3)) ; upar.Fix(3) ;
		upar.Release(0) ; upar.SetError(0,60.);

		std::cout << "=============================================" << std::endl;
		std::cout<<"Second Minimization: "<<upar<<std::endl;
		std::cout << "=============================================" << std::endl;
		MnMigrad mnr( *fit , upar , MnStrategy(0)) ;
		min = mnr() ;

		//Status report
		if (min.IsValid()) std::cout << "Fit success"         << std::endl ;
		else               std::cout << "Fit failed"   << std::endl ;
		std::cout << "Total time: " << total_timeTaken.total_seconds() << std::endl ;
		std::cout << "MINIMIZATION OUTCOME: " <<  min  << std::endl ;

		//Release parameter 0,1,2,3 minimize again
		for ( int i=0 ; i < 4 ; i++ ) { upar.Release(i) ; upar.SetError(i,60.); }

		std::cout << "=============================================" << std::endl;
		std::cout<<"0-3 par free: "<<upar<<std::endl;
		std::cout << "=============================================" << std::endl;
		MnMigrad mn3( *fit , upar , MnStrategy(0)) ;
		min = mn3() ;

		//Status report
		if (min.IsValid()) std::cout << "Fit success"         << std::endl ;
		else               std::cout << "Fit failed"   << std::endl ;
		std::cout << "Total time: " << total_timeTaken.total_seconds() << std::endl ;
		std::cout << "MINIMIZATION OUTCOME: " <<  min  << std::endl ;



		//Get the fitting parameters
		for (uint i=0; i < parIniSize;i++) {
			parIni[i]=min.UserState().Value(i);
			parErr[i]=min.UserState().Error(i);
		}
	}
	/*********Finish Irradiated/undepleted fit******************************************/
	/***********************************************************************/


	/*********Begin depleted fit. Non-irradiated dectectors**************/
	/***********************************************************************/
	//Fitting normalizator, Vdep, depth, Capacitance if detector is depleted and non-irradiated


	if (fluence == 0) /*&& (vBias >= vDep)*/{

		irradiated = false;

		fitParamNorm = TRACSsim[0]->get_fitNorm();
		fitParamVdep = TRACSsim[0]->get_vDep();
		fitParamDepth = TRACSsim[0]->get_depth();
		fitParamCapac = TRACSsim[0]->get_capacitance();

		parIni = {fitParamNorm, fitParamVdep, fitParamDepth, fitParamCapac};
		parIniSize = parIni.size() ;

		parErr.resize(parIniSize);
		for (size_t i=0; i < parIniSize; i++)
			parErr[i] = 60.;

		//Pass parameters to Minuit
		MnUserParameters upar(parIni,parErr) ;
		for ( int i=0 ; i < parIniSize ; i++ ) {
			char pname[parIniSize]; sprintf( pname , "p%d" , i);
			upar.SetName( i , pname );
		}

		//upar.Fix(0) ; //Normalizator
		upar.Fix(1); //Depletion voltage
		//upar.Fix(2); //Detector depth
		//upar.Fix(3); //Capacitance

		std::cout << "=============================================" << std::endl;
		std::cout << "Initial parameters: "<<upar<<std::endl;
		std::cout << "=============================================" << std::endl;
		std::cout << "=============================================" << std::endl;
		std::cout << "tolerance= " << TRACSsim[0]->GetTolerance()    << std::endl;
		std::cout << "chiFinal= " << TRACSsim[0]->GetchiFinal()      << std::endl;
		std::cout << "=============================================" << std::endl;
		std::cout << "=============================================" << std::endl;


		//Do the minimization

		ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);
		ROOT::Math::MinimizerOptions::SetDefaultTolerance(TRACSsim[0]->GetTolerance());
		MnMigrad mn( *fit , upar , MnStrategy(0)) ;
		FunctionMinimum min = mn() ;

		//Status report
		if (min.IsValid()) std::cout << "Fit success"         << std::endl ;
		else               std::cout << "Fit failed"   << std::endl ;
		std::cout << "Total time: " << total_timeTaken.total_seconds() << std::endl ;
		std::cout << "MINIMIZATION OUTCOME: " <<  min  << std::endl ;

		//Get the fitting parameters
		for (uint i=0; i < parIniSize;i++) {
			parIni[i]=min.UserState().Value(i);
			parErr[i]=min.UserState().Error(i);
		}
	}


	/*********Finish depleted non-irradiated**********************************************/
	/***********************************************************************/

	//Calculate TCT pulses with the fit output parameters
	for (int i = 0; i < num_threads; ++i) {
		if (irradiated) t[i] = std::thread(call_from_thread_FitPar, i, std::ref(carrierThread_fileNames[i]), vector_zValues, vector_yValues, vector_voltValues ,parIni);
		else t[i] = std::thread(call_from_thread_FitNorm, i, std::ref(carrierThread_fileNames[i]), vector_zValues, vector_yValues, vector_voltValues ,parIni);

	}

	for (int i = 0; i < num_threads; ++i) {
		t[i].join();
	}

	//Dump tree to disk
	TFile fout("output.root","RECREATE") ;
	TTree *tout = new TTree("edge","Fitting results");

	TMeas *emo = new TMeas( );
	emo->Nt   = TRACSsim[0]->GetnSteps() ;
	emo->volt = new Double_t [emo->Nt] ;
	emo->time = new Double_t [emo->Nt] ;
	emo->Qt   = new Double_t [emo->Nt] ;

	// Create branches
	tout->Branch("raw", &emo,32000,0);

	//Read RAW file
	TRACSsim[0]->DumpToTree( emo , tout ) ;

	//TRACSsim[0]->GetTree( tsim );
	fout.Write();
	delete tout ;
	fout.Close();
	delete emo ;

	//Clean
	for (int i = 0 ; i < num_threads ; i++){
				const char * c = carrierThread_fileNames[i].c_str();
				remove(c);
			}
	for (uint i = 0; i < TRACSsim.size(); i++)	{
		delete TRACSsim[i];
	}

	delete fit;
	std::quick_exit(1);
}

//_____________________________________________________________________

Double_t TRACSFit::operator() ( const std::vector<Double_t>& par  ) const {

	static int icalls ;
	boost::posix_time::ptime start = boost::posix_time::second_clock::local_time();

	for (int i = 0 ; i < vItotals.size() ; i++){
		for (int j = 0 ; j < vItotals[i].size() ; j++)
			vItotals[i][j] = 0;
	}

	for (int i = 0; i < num_threads; ++i) {
		if (irradiated) t[i] = std::thread(call_from_thread_FitPar, i, std::ref(carrierThread_fileNames[i]), vector_zValues, vector_yValues, vector_voltValues ,par);
		else t[i] = std::thread(call_from_thread_FitNorm, i, std::ref(carrierThread_fileNames[i]), vector_zValues, vector_yValues, vector_voltValues ,par);
	}

	for (int i = 0; i < num_threads; ++i) {
		t[i].join();
	}

	Double_t chi2 = fit->LeastSquares( ) ;
	boost::posix_time::ptime end = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration timeTaken = end - start;
	total_timeTaken += timeTaken;

	std::cout << "-------------------------------------------------------------------------------------> " << std::endl;
	std::cout << "----------------------------> Time taken for chi2 calculation (milliseconds): " << timeTaken.total_milliseconds() << std::endl;
	std::cout << "-------------------------------------------------------------------------------------> " << std::endl;
	std::cout << "----------------------------> icalls="<<icalls<<" chi2=" << chi2 << "\t" ;
	for (uint ipar=0 ; ipar<par.size() ; ipar++) std::cout << "p["<<ipar<<"]="<<par[ipar]<<"\t" ; std::cout << std::endl;
	std::cout << "-------------------------------------------------------------------------------------> " << std::endl;
	icalls++;

	for (uint i = 0; i < TRACSsim.size(); i++)	{
		delete TRACSsim[i];
	}

	return chi2 ;

	for (int i = 0 ; i < vItotals.size() ; i++){
		for (int j = 0 ; j < vItotals[i].size() ; j++)
			vItotals[i][j] = 0;
	}


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


