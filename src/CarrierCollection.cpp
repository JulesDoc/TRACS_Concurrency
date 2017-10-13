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

/************************************CARRIER COLLECTION***********************************
 *
 * Carrier collection get carriers information from a file and store it in a vector of carriers.
 * This vector of carriers will be used by each of the threads in each of the steps of the TCT scan.
 *
 *
 */

#include "../include/CarrierCollection.h"

double extra_y;

CarrierCollection::CarrierCollection(SMSDetector * detector) :
						_detector(detector)
{

}
/*
 * Adds carriers from the file. The name of the file must be set in the steering file, in the field CarrierFile, more instructions can be found there.
 * The calculation of average positions is neccesary for fitting purposes, DumpToTree function.
 */
/**
 *
 * @param filename
 */
void CarrierCollection::add_carriers_from_file(const std::string& filename, const std::string& scanType, double depth)
{
	// get char representation and make ifstream
	//char * char_fn = filename.toLocal8Bit().data();
	std::ifstream infile(filename);
	bool once = true;
	char carrier_type;
	double q, x_init, y_init, gen_time;

	// process line by line
	std::string line;


	//Preprocessing for fitting bottom scan_type: Fixing mismatch between detector depth and y_init.
	//Outside the loop for performance purpose
	if (scanType == "bottom"){
		std::getline(infile, line);
		std::istringstream iss(line);
		if (!(iss >> carrier_type >> q >> x_init >> y_init >> gen_time)) {
			//std::cout << "Error while reading file" << std::endl;
		}
		//Extra in micrometers to shift y_init position
		extra_y = depth - y_init;

		//Calculate average beam position
		beamy += x_init;
		beamz += y_init + extra_y;
		Carrier carrier(carrier_type, q, x_init, y_init + extra_y, _detector, gen_time);
		_carrier_list_sngl.push_back(carrier);
		if ( _carrier_list_sngl.size()!=0 ) {
			beamy = beamy / _carrier_list_sngl.size();
			beamz = beamz / _carrier_list_sngl.size();
		}
	}
	while (std::getline(infile, line))
	{

		std::istringstream iss(line);
		if (!(iss >> carrier_type >> q >> x_init >> y_init >> gen_time)) { 
			//std::cout << "Error while reading file" << std::endl;
			break;
		} 


		//Calculate average beam position
		beamy += x_init;
		beamz += y_init + extra_y;
		Carrier carrier(carrier_type, q, x_init, y_init + extra_y, _detector, gen_time);
		_carrier_list_sngl.push_back(carrier);
	}
	if ( _carrier_list_sngl.size()!=0 ) {
		beamy = beamy / _carrier_list_sngl.size();
		beamz = beamz / _carrier_list_sngl.size();
	}



}

/*
 * From the carrier list taken from the file of carriers, this method runs through the vector of carriers and calls the simulate.drift for each one of them, depending if the carrier
 * is a hole or an electron. Trapping effects are included at the end directly in the valarry of currents.
 * Info related to diffusion is displayed using this method as well. Whenever a carrier crosses to depleted region, the acumulative variable totalCross grows.
 */
/**
 *
 * @param dt
 * @param max_time
 * @param shift_x
 * @param shift_y
 * @param curr_elec
 * @param curr_hole
 * @param totalCrosses
 */
void CarrierCollection::simulate_drift( double dt, double max_time, double shift_x /*yPos*/, double shift_y /*zPos*/,
		std::valarray<double>&curr_elec, std::valarray<double> &curr_hole, int &totalCrosses)
{

	double x_init, y_init;
	int totalCross = 0;
	bool control = true;

	//fileDiffDrift.open ("fileDiffDrift");
	// range for through the carriers
	for (auto carrier : _carrier_list_sngl)
	{
		char carrier_type = carrier.get_carrier_type();
		// simulate drift and add to proper valarray
		if (carrier_type == 'e')
		{

			// get and shift carrier position
			std::array< double,2> x = carrier.get_x();
			x_init = x[0] + shift_x;
			y_init = x[1] + shift_y;

			//x[0] represents the X position read from the carriers file
			//shift_x represents the shift applied to X read from the steering file, namely, where the laser points to.
			//x[1] represents the Y position read from the carriers file. Y is seen as Z.
			//shift_y represents the Z (y) position read from the steering file. Since the program (usually) does edge-TCT, Z can be defined in one or more steps.

			curr_elec += carrier.simulate_drift( dt , max_time, x_init, y_init);
		}

		else if (carrier_type =='h')
		{
			// get and shift carrier position
			std::array< double,2> x = carrier.get_x();
			double x_init = x[0] + shift_x;
			double y_init = x[1] + shift_y ;

			curr_hole += carrier.simulate_drift( dt , max_time, x_init, y_init);
		}
		//Let's see how many carriers from the carrier list for this step in Z have crossed to the depleted region.
		//A flag is switched to true on carrier.simulate_drift when a carrier filfull the requirements. See carrier.simulate_drift
		if (carrier.crossed()){
			totalCross += 1;
		}
	}
	//fileDiffDrift.close();
	//std::cout << "Number of carriers crossed to DR in last Z step with Height " << shift_y << ": " << totalCross << std::endl;
	totalCrosses += totalCross;

	double trapping_time = _detector->get_trapping_time();

	for (double i = 0.; i < curr_hole.size(); i ++)
	{
		double elapsedT = i*dt;
		curr_elec[i] *= exp(-elapsedT/trapping_time);
		curr_hole[i] *= exp(-elapsedT/trapping_time);
	}

}
/**
 *
 * @param n_bins_x
 * @param n_bins_y
 * @param hist_name
 * @param hist_title
 * @return
 */
TH2D CarrierCollection::get_e_dist_histogram(int n_bins_x, int n_bins_y,  TString hist_name, TString hist_title)
{
	// get detector limits
	double x_min = _detector->get_x_min();
	double x_max = _detector->get_x_max();
	double y_min = _detector->get_y_min();
	double y_max = _detector->get_y_max();

	// create histogram object
	TH2D e_dist = TH2D(hist_name, hist_title, n_bins_x , x_min, x_max, n_bins_y, y_min, y_max);

	// range for through the carriers and fill the histogram
	for (auto carrier : _carrier_list_sngl)
	{
		char carrier_type = carrier.get_carrier_type();
		if (carrier_type == 'e')
		{
			std::array< double,2> x = carrier.get_x();
			double q = carrier.get_q();
			e_dist.Fill(x[0], x[1], q);
		}
	}
	return e_dist;
}
/**
 *
 * @param n_bins_x
 * @param n_bins_y
 * @param shift_x
 * @param shift_y
 * @param hist_name
 * @param hist_title
 * @return
 */
TH2D CarrierCollection::get_e_dist_histogram(int n_bins_x, int n_bins_y, double shift_x, double shift_y, TString hist_name, TString hist_title)
{
	// get detector limits
	double x_min = _detector->get_x_min();
	double x_max = _detector->get_x_max();
	double y_min = _detector->get_y_min();
	double y_max = _detector->get_y_max();

	// create histogram object
	TH2D e_dist = TH2D(hist_name, hist_title, n_bins_x , x_min, x_max, n_bins_y, y_min, y_max);

	// range for through the carriers and fill the histogram
	for (auto carrier : _carrier_list_sngl)
	{
		char carrier_type = carrier.get_carrier_type();
		if (carrier_type == 'e')
		{
			std::array< double,2> x = carrier.get_x();
			double q = carrier.get_q();
			e_dist.Fill(x[0]+shift_x, x[1]+shift_y, q);
		}
	}
	return e_dist;
}

/*
 ********************** DESTRUCTOR OF THE CLASS CARRIER	COLLECTION **************************
 */
CarrierCollection::~CarrierCollection()
{

}
