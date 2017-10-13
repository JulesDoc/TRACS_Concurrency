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

/************************************TRACSInterface***********************************
 *
 * High-level class used as interface between files and functionalities.
 *
 */


#include <TRACSInterface.h>
#include <mutex>          // std::mutex


std::mutex mtx2;

/*
 * The constructor mainly initializes all the values that will be used during the execution. Firstly it read most of them from the steering file by means of a parsing method inside utilities class.
 * Another important task carrying out here is the definition of the vectors and coordinates. Vectors to store currents and coordinates to define the scanning, positions, steps...
 * When the information is set up, a new instance of a detector can be launched, one per TRACSInterface object, avoiding data races and concurrent access to same positions memory when solving field equations in the
 * detector. Making copy constructor and passing copies of detector instances is expensive for the program and do not work, neither improve the performance. Same behavior is used with the current collection.
 * The only global variable shared by threads is the one that stores the total induce current vItotals, but it is important to note that threads access uniquely to their tid position in the vector.
 */
/**
 *
 * @param filename
 */
TRACSInterface::TRACSInterface(std::string filename, const std::string& carrFile, const std::vector<double>& zVector, const std::vector<double>& yVector,
		const std::vector<double>& voltVector):zVector(zVector), yVector(yVector), voltVector(voltVector), carrierFile(carrFile)
{
	neff_param = std::vector<double>(8,0);
	total_crosses = 0;

	utilities::parse_config_file(filename, depth, width,  pitch, nns, temp, trapping, fluence, nThreads, n_cells_x, n_cells_y, bulk_type,
			implant_type, waveLength, scanType, C, dt, max_time, vInit, deltaV, vMax, vDepletion, zInit, zMax, deltaZ, yInit, yMax, deltaY, neff_param, neffType,
			tolerance, chiFinal, diffusion, fitNorm/*, gen_time*/);

	if (fluence == 0) // if no fluence -> no trapping
	{
		//Trapping configuration
		trapping = std::numeric_limits<double>::max();
		trap = "NOtrapping";
		start = "NOirrad";
	}
	else{
		//trapping configuration
		trap = std::to_string((int) std::floor(1.e9*trapping));
		start = "irrad";

	}

	//Following variables are used for writing data in the output file.
	//Convert relevant simulation numbers to string for fileNaming
	dtime = std::to_string((int) std::floor(dt*1.e12));
	neighbrg_strips = std::to_string(nns);
	stepV = std::to_string((int) std::floor(deltaV));
	stepZ = std::to_string((int) std::floor(deltaZ));
	stepY = std::to_string((int) std::floor(deltaY));
	capacitance = std::to_string((int) std::floor(C*1.e12));
	voltage = std::to_string((int) std::floor(vInit));




	//Dolfin instruction por mesh boundary extrapolation
	parameters["allow_extrapolation"] = true;

//	mtx2.lock(); //Thread-safe for dolfin data type
	std::unique_lock<std::mutex> guard(mtx2);
	detector = new SMSDetector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y, temp, trapping, fluence, neff_param, neffType, diffusion, dt);
	guard.unlock();

	carrierCollection = new CarrierCollection(detector);
	//QString carrierFileName = QString::fromUtf8(carrFile.c_str());
	carrierCollection->add_carriers_from_file(carrierFile, scanType, depth);

	n_tSteps = (int) std::floor(max_time / dt);
	//currents vectors used to store temporal and final values.
	i_elec.resize((size_t) n_tSteps);
	i_hole.resize ((size_t) n_tSteps);
	i_total.resize((size_t) n_tSteps);
	i_shaped.resize((size_t) n_tSteps);
	i_temp.resize((size_t) n_tSteps);

	for (int k = 1 ; k < n_tSteps; k++ ){
		i_temp[k] = 0 ;
	}




	vBias = vInit;
	set_tcount(0);

	i_ramo  = NULL;
	i_rc    = NULL;
	i_conv  = NULL;

}


// Destructor
TRACSInterface::~TRACSInterface()
{
	delete i_ramo;
	delete i_rc;
	delete i_conv;
	delete carrierCollection;
	delete detector;
}

/*
 * Convert i_total to TH1D. ROOT based method.
 */
TH1D * TRACSInterface::GetItRamo()
{
	if (i_ramo != NULL) 
	{
	}
	else
	{
		//TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,10);
		//float r = f1->GetRandom();
		TString htit, hname;
		htit.Form("ramo_%d_%d", tcount, count1);
		hname.Form("Ramo_current_%d_%d", tcount, count1);
		i_ramo  = new TH1D(htit,hname,n_tSteps, 0.0, max_time);
		//std::cout << htit << std::endl;

		// Compute time + format vectors for writting to file
		for (int j=0; j < n_tSteps; j++)
		{
			i_ramo->SetBinContent(j+1, i_total[j] );
		}
		count1++;

	}
	return i_ramo;
}

/*
 * Convert i_total to TH1D after simulating simple RC circuit. ROOT based method.
 */
TH1D * TRACSInterface::GetItRc()
{

	double RC = 50.*C; // Ohms*Farad
	double alfa = dt/(RC+dt);


	for (int j = 1; j <n_tSteps; j++)
	{
		i_shaped[j]=i_shaped[j-1]+alfa*(i_total[j]-i_shaped[j-1]);

	}
	count2++;


}

/*
 * Convert i_total to TH1D after convolution with the amplifier TransferFunction. ROOT based method.
 */
TH1D * TRACSInterface::GetItConv()
{
	if (i_conv != NULL) 
	{
	}
	else
	{
		//TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,10);
		//float r = f1->GetRandom();
		TString htit, hname;
		htit.Form("ramo_conv_%d_%d", tcount, count3);
		hname.Form("Ramo current_%d_%d", tcount, count3);
		i_conv  = new TH1D(htit,hname,n_tSteps, 0.0, max_time);
		i_ramo = GetItRamo();
		//mtx2.lock();
		i_conv = H1DConvolution( i_ramo , C*1.e12, tcount );
		//mtx2.unlock();
		count3++;
	}
	return i_conv;
}

/*
 * Performs the simulation for all given carriers and stores the current in a valarray.
 * No variable is returned. To get the current one must choose the apropiate Getter for 
 * one's needs
 */
void TRACSInterface::simulate_ramo_current()
{
	i_rc = nullptr;
	i_ramo = nullptr;
	i_conv = nullptr;
	i_hole = 0;
	i_elec = 0;
	i_total = 0;

	/*Important for Diffusion:
	yPos is later transformed into X.
	zPos is later transformed into Y.*/
	carrierCollection->simulate_drift( dt, max_time, yPos, zPos, i_elec, i_hole, total_crosses);
	i_total = i_elec + i_hole;
}

double TRACSInterface::get_Itemp(int i , int j){

	//return i_temp[i][j];
	return 0;
}

double TRACSInterface::get_vDep(){
	return vDepletion;
}
double TRACSInterface::get_fitNorm(){
	return fitNorm;
}

double TRACSInterface::get_depth(){
	return depth;
}

double TRACSInterface::get_capacitance(){
	return C;
}

double TRACSInterface::get_fluence(){
	return fluence;
}

double TRACSInterface::get_vBias(){
	return vInit;
}

std::string TRACSInterface::get_neff_type(){
	return neffType;
}

double TRACSInterface::get_dt(){
	return dt;
}

/*
 *
 * Returns Neff parametrization
 *
 */
std::vector<double>TRACSInterface::get_NeffParam()
{
	return neff_param;
} 

int TRACSInterface::GetnSteps(){
	return n_tSteps;
}

double TRACSInterface::GetTolerance(){
	return tolerance;
}

double TRACSInterface::GetchiFinal(){
	return chiFinal;
}

int TRACSInterface::GettotalCrosses(){
	return total_crosses;
}

UShort_t TRACSInterface::GetYear(){
	time_t currentTime;
	struct tm *localTime;

	time( &currentTime );                   // Get the current time
	localTime = localtime( &currentTime );  // Convert the current time to the local time
	year = localTime->tm_year + 1900;
	return year;
}

UShort_t TRACSInterface::GetMonth(){
	time_t currentTime;
	struct tm *localTime;

	time( &currentTime );                   // Get the current time
	localTime = localtime( &currentTime );  // Convert the current time to the local time
	month = localTime->tm_mon + 1;
	return month;
}
UShort_t TRACSInterface::GetDay(){
	time_t currentTime;
	struct tm *localTime;

	time( &currentTime );                   // Get the current time
	localTime = localtime( &currentTime );  // Convert the current time to the local time
	day = localTime->tm_mday;
	return day;
}
UShort_t TRACSInterface::GetHour(){
	time_t currentTime;
	struct tm *localTime;

	time( &currentTime );                   // Get the current time
	localTime = localtime( &currentTime );  // Convert the current time to the local time
	hour = localTime->tm_hour;
	return hour;
}
UShort_t TRACSInterface::GetMinute(){
	time_t currentTime;
	struct tm *localTime;

	time( &currentTime );                   // Get the current time
	localTime = localtime( &currentTime );  // Convert the current time to the local time
	min = localTime->tm_min;
	return min;
}
UShort_t TRACSInterface::GetSecond(){
	time_t currentTime;
	struct tm *localTime;

	time( &currentTime );                   // Get the current time
	localTime = localtime( &currentTime );  // Convert the current time to the local time
	sec = localTime->tm_sec;
	return sec;
}

/*
 * Sets the trapping time in the detector to the input value. 
 * Remember that the trapping time must be a positive number.
 * The smaller the trapping time the bigger the signal loss.
 */
/**
 *
 * @param newTrapTime
 */
void TRACSInterface::set_trappingTime(double newTrapTime)
{
	trapping = newTrapTime;
	detector->set_trapping_time(trapping);
}

/*
 * Sets how much the carriers will be displaced in the Z axis from its original 
 * position in the file read by TRACS. Note that if the carriers are not inside 
 * the detector they will not produce current. This is relevant mainly for 
 * edge-TCT simulations.
 */
/**
 *
 * @param newZPos
 */
void TRACSInterface::set_zPos(double newZPos)
{
	zPos = newZPos;
}

/*
 * Sets how much the carriers will be displaced in the Y axis from its original 
 * position in the file read by TRACS. Note that if the carriers are not inside 
 * the detector they will not produce current. This is used to center the red 
 * pulse in redTCT and to center the focus in TPA and edgeTCT
 */
/**
 *
 * @param newYPos
 */
void TRACSInterface::set_yPos(double newYPos)
{
	if (std::abs(newYPos) > (2*nns+1)*pitch)
	{
		std::cout << "Watch out! You probably set the laser out of the detector" << std::endl;
	}
	yPos = newYPos;
}

/*
 * Sets bias voltages in the detector, fields should be recalculated again 
 * before simulating any transients
 */
/**
 *
 * @param newVBias
 */
void TRACSInterface::set_vBias(double newVBias)
{
	vBias = newVBias;
	detector->set_voltages(vBias, vDepletion);
}

/*
 * Sets a number (current thread).
 * Used to index the different output files 
 *
 */
/**
 *
 * @param tid
 */
void TRACSInterface::set_tcount(int tid)
{
	tcount = tid;
}

/*
 * Sets the desired Neff parameters in the detector. Fields should be calculated
 * again before simulating any current. Note that different neff parametrizations
 * use different parameters so not all may be used at once.
 */
/**
 *
 * @param newParam
 */
void TRACSInterface::set_FitParam(std::vector<double> newFitParam)
{

	for (uint i = 0 ; neff_param.size(); i++)
	{
		neff_param[i] = newFitParam[i];
	}
	fitNorm = newFitParam[8];
	depth = newFitParam[9];
	delete detector;
	std::unique_lock<std::mutex> guard(mtx2);
	detector = new SMSDetector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y, temp, trapping, fluence, neff_param, neffType, diffusion, dt);
	guard.unlock();
	//detector->setFitParameters(newFitParam);

}

void TRACSInterface::set_Fit_Norm(std::vector<double> vector_fitTri)
{

	fitNorm = vector_fitTri[0];
	vDepletion = vector_fitTri[1];
	depth = vector_fitTri[2];
	C = vector_fitTri[3];
	delete detector;
	std::unique_lock<std::mutex> guard(mtx2);
	detector = new SMSDetector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y, temp, trapping, fluence, neff_param, neffType, diffusion, dt);
	guard.unlock();
}



/*
 * Calculates the electric field and potential inside the detector. It is 
 * required after any modification of the Neff or the bias voltage applied. 
 * Weighting field and potential need not be calculated again since they 
 * are independent on those parameters.
 */



void TRACSInterface::calculate_fields()
{
	// Get detector ready

	detector->solve_w_u();
	detector->solve_d_u();
	detector->solve_w_f_grad();
	detector->solve_d_f_grad();
	detector->get_mesh()->bounding_box_tree();

}

/*
 * Change parametrization of the Neff. Possibilities are: Trilinear (default)
 * Linear, Triconstant. More information on this three parametrizarions can
 * be found in the documentation (Config.TRACS and README.md)
 */
/**
 *
 * @param newParametrization
 */
void TRACSInterface::set_neffType(std::string newParametrization)
{
	neffType = newParametrization;
	detector->set_neff_type(neffType);

}

/*
 * Allows the user to select a new carrier distribution from a file that
 * complies with the TRACS format. 
 */
/**
 *
 * @param newCarrFile
 */
void TRACSInterface::set_carrierFile(std::string newCarrFile)
{
	//QString carrierFileName = QString::fromUtf8(newCarrFile.c_str());
	carrierCollection->add_carriers_from_file(newCarrFile, scanType, depth);
}

/**
 *
 * @param tid
 */


void TRACSInterface::loop_on(int tid)
{

	int index_total = 0;
	int i,j;
	i = 0; j = 0;
	std::unique_lock<std::mutex> guard(mtx2);
	detector->solve_w_u();
	guard.unlock();

	if (scanType == "edge"){
		//Voltage scan
		for (int index_volt = 0; index_volt < voltVector.size() ; index_volt++){

			detector->set_voltages(voltVector[index_volt], vDepletion);
			//mtx2.lock();
			//calculate_fields();
			//mtx2.unlock();
			
			std::unique_lock<std::mutex> guard(mtx2);
			detector->solve_d_u();
			guard.unlock();

			detector->solve_w_f_grad();
			detector->solve_d_f_grad();
			detector->get_mesh()->bounding_box_tree();


			for (int index_zscan = 0; index_zscan < zVector.size(); index_zscan++){

				//simulate_ramo_current();
				carrierCollection->simulate_drift( dt, max_time, yInit, zVector[index_zscan], i_elec, i_hole, total_crosses);
				i_total = i_elec + i_hole;

				GetItRc();
				std::unique_lock<std::mutex> guard(mtx2);
				i_temp = vItotals[index_total];
				vItotals[index_total] = i_shaped + i_temp;
				guard.unlock();

				index_total++;
				i_total = 0 ; i_elec = 0; i_hole = 0; i_shaped = 0; i_temp = 0;


			}

			if (tid == 0) fields_hist_to_file(tid, index_volt);
		}
		//for (int i = 0 ; i < vItotals.size() ; i++){
		//			for (int j = 0 ; j < vItotals[i].size() ; j++)
		//				std::cout << "i " << i << "; j " << j << "    " <<  vItotals[i][j] << std::endl;
		//		}


	}


	if (scanType == "top" || scanType == "bottom"){
		//Voltage scan
		for (int index_volt = 0; index_volt < voltVector.size() ; index_volt++){

			detector->set_voltages(voltVector[index_volt], vDepletion);
			
			std::unique_lock<std::mutex> guard(mtx2);
			detector->solve_d_u();
			guard.unlock();

			detector->solve_w_f_grad();
			detector->solve_d_f_grad();
			detector->get_mesh()->bounding_box_tree();


			for (int index_yscan = 0; index_yscan < yVector.size(); index_yscan++){

				carrierCollection->simulate_drift( dt, max_time, yVector[index_yscan], zInit, i_elec, i_hole, total_crosses);
				i_total = i_elec + i_hole;

				GetItRc();
				
				std::unique_lock<std::mutex> guard(mtx2);
				i_temp = vItotals[index_total];
				vItotals[index_total] = i_shaped + i_temp;
				guard.unlock();

				index_total++;
				i_total = 0 ; i_elec = 0; i_hole = 0; i_shaped = 0; i_temp = 0;

			}

			if (tid == 0) fields_hist_to_file(tid, index_volt);
		}




	}

}

/*
 *
 * Write to file header. The input int is used to label files (multithreading)!
 *
 *
 *
 */
/**
 *
 * @param tid
 */
void TRACSInterface::write_header(int tid)
{
	// Convert Z to milimeters
	std::vector<double> aux_zsh(zVector.size());
	aux_zsh = zVector;
	std::transform(aux_zsh.begin(), aux_zsh.end(), aux_zsh.begin(), std::bind1st(std::multiplies<double>(),(1./1000.)));

	// Convert Z to milimeters
	std::vector<double> aux_ysh(yVector.size());
	aux_ysh = yVector;
	std::transform(aux_ysh.begin(), aux_ysh.end(), aux_ysh.begin(), std::bind1st(std::multiplies<double>(),(1./1000.)));
	hetct_rc_filename = start+"_dt"+dtime+"ps_"+capacitance+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neighbrg_strips+"nns_"+scanType+"_"+std::to_string(tcount)+"_rc.hetct";
	// write header for data analysis
	utilities::write_to_hetct_header(hetct_rc_filename, detector, C, dt, aux_ysh, aux_zsh, waveLength, scanType, carrierFile, voltVector);

}
/*
 *Resizing for storing data in memory 
 *and using it later/writing to a single output file.
 *
 */
void TRACSInterface::resize_array()
{

}
/*
 * Writing to a single file
 *
 *
 */
/**
 *
 * @param tid
 */
void TRACSInterface::write_to_file(int tid)
{

	int index_rc = 0;
	if (scanType == "edge"){

		for (int i = 0 ; i < voltVector.size() ;  i++){
			for (int j = 0 ; j < zVector.size() ; j++){
				utilities::write_to_file_row(hetct_rc_filename, i_rc_array[index_rc], detector->get_temperature(), yInit, zVector[j], voltVector[i]);
				index_rc++;
			}

		}

	}



	if (scanType == "top" || scanType == "bottom"){

		for (int i = 0 ; i < voltVector.size() ;  i++){
			for (int j = 0 ; j < yVector.size() ; j++){
				utilities::write_to_file_row(hetct_rc_filename, i_rc_array[index_rc], detector->get_temperature(), yVector[j], zInit, voltVector[i]);
				index_rc++;
			//std::cout << i_rc_array[i][j] << std::endl;
			}

		}

	}


}

void TRACSInterface::fields_hist_to_file(int tid, int vPos)
{

	/*Exporting 2D Histograms to file*/
	// get plot and set new data
	TString file_name;
	file_name.Form("wf%dV", TMath::Nint(voltVector[vPos]));//"2Dhistos"+voltages[vPos]+"V"+".root";
	TFile *fout = new TFile(file_name,"RECREATE");

	TString wpm;
	wpm.Form("WP_%d_V", TMath::Nint(voltVector[vPos]) ) ;
	TString wfm;
	wfm.Form("WF_%d_V", TMath::Nint(voltVector[vPos]) ) ;
	TString efm;
	efm.Form("EF_%d_V", TMath::Nint(voltVector[vPos]) ) ;

	TH2D h_w_u      = utilities::export_to_histogram    ( *detector->get_w_u(), "h_w_u", wpm /*"Weighting potential"*/, detector->get_n_cells_x(), detector->get_x_min(), detector->get_x_max(), detector->get_n_cells_y(), detector->get_y_min(), detector->get_y_max());
	TH2D h_w_f_grad = utilities::export_mod_to_histogram( *detector->get_w_f_grad(), "h_w_f_grad", wfm /*"Weighting field"*/, detector->get_n_cells_x(), detector->get_x_min(), detector->get_x_max(), detector->get_n_cells_y(), detector->get_y_min(), detector->get_y_max());
	TH2D h_d_f_grad = utilities::export_mod_to_histogram( *detector->get_d_f_grad(), "h_d_f_grad", efm /*"Electric field"*/, detector->get_n_cells_x(), detector->get_x_min(), detector->get_x_max(), detector->get_n_cells_y(), detector->get_y_min(), detector->get_y_max());

	h_w_u.Write();
	h_w_f_grad.Write();
	h_d_f_grad.Write();
	fout->Close();
}






/* ---------------------------------------------------------------- */
/**
 *
 * @param tree
 */
void TRACSInterface::GetTree( TTree *tree ) {

	//TFile *f=new TFile("test.root","RECREATE") ;

	// Create a ROOT Tree
	//tree->SetDirectory(0);

	// Create a pointer to an raw data object
	TMeas *em = new TMeas( );
	em->Nt = n_tSteps ;
	em->volt = new Double_t [n_tSteps] ;
	em->time = new Double_t [n_tSteps] ;
	em->Qt = new Double_t [n_tSteps] ;

	// Create branches
	tree->Branch("raw", &em,32000,0);

	//Read RAW file
	DumpToTree( em , tree ) ;
	//tree->Draw("volt-BlineMean:time","event==0","l"); gPad->Modified();gPad->Update();
	//f->Write();
	//f->Close();


	delete em ;

}
/*----------------------------------------------------------*/
/**
 *
 * @param em
 * @param tree
 */
void TRACSInterface::DumpToTree( TMeas *em , TTree *tree ) {

	//Time information
	TString sdate; //= TString( line ) ;
	UShort_t  dd, mm, yy , hh , mn, ss;
	yy=GetYear();
	mm=GetMonth();
	dd=GetDay();
	hh=GetHour();
	mn=GetMinute();
	ss=GetSecond();
	TDatime date ;
	date.Set(yy, mm, dd, hh, mn, ss) ;
	em->utc = date ;

	//Total number of Scans

	Int_t NumOfScans = 1 ;
	/*if ( n_vSteps !=0 ) NumOfScans = NumOfScans * (1 + n_vSteps) ;
	if ( n_zSteps !=0 ) NumOfScans = NumOfScans * (1+ n_zSteps) ;
	if ( n_ySteps !=0 ) NumOfScans = NumOfScans * (1+ n_ySteps);*/
	NumOfScans = NumOfScans * voltVector.size();
	NumOfScans = NumOfScans * zVector.size();
	NumOfScans = NumOfScans * yVector.size();


	// NVoltages
	em->NV = voltVector.size();//n_vSteps + 1 ;

	//Create a TMeasHeader object to store info that does not depend on event
	TMeasHeader *emh=new TMeasHeader( em->NV ) ;
	emh->Lambda = waveLength ;
	emh->NV = em->NV;
	emh->comment=TString( "TRACS simulated data" ) ;
	emh->Setup = 5 ;
	emh->Fluence = fluence ;
	emh->Nav  = 1 ;
	emh->Gain = 1.0 ;
	emh->iann = 0. ;
	emh->Illum = -2 ;

	//Cend
	emh->Cend = 0. ;

	//Vbias vector
	for ( int i=0 ; i< em->NV ; i++) emh->vVbias[i] = voltVector[i] ;

	//Nominal power
	emh->Power = 1.0 ;

	//Ax
	emh->Ax = 0. ;
	emh->Nx = 1 ;

	//Ay
	em->Ay  = emh->Ay = deltaY  ;
	em->Ny  = emh->Ny = yVector.size();//n_ySteps + 1 ;

	//Az
	em->Az  = emh->Az = deltaZ  ;
	em->Nz  = emh->Nz = zVector.size(); //n_zSteps + 1;

	em->Temp = temp ;

	int  iRead=0 , iactual = 1 ;
	int  Polarity = 0 ;
	Double_t x0=0., y0=0., z0=0., t0s;
	em->x = 0. ;
	for (Int_t iloop = 0 ; iloop < NumOfScans ; iloop++ ) {

		Int_t tms = iloop * 4000 ;
		em->Vbias = detector->get_vbias();
		em->y = carrierCollection->beamy;
		em->z = carrierCollection->beamy;

		em->Itot = 0. ;

		//Calculate the bin slice along each coordinate
		if (iloop==0) { x0=em->x ;y0=em->y ;z0=em->z ; }
		em->ix = (em->Ax!=0.)? 1 + TMath::Nint((em->x - x0)/(em->Ax)):1 ;
		em->iy = (em->Ay!=0.)? 1 + TMath::Nint((em->y - y0)/(em->Ay)):1 ;
		em->iz = (em->Az!=0.)? 1 + TMath::Nint((em->z - z0)/(em->Az)):1 ;

		//Calculate the time increment wrt t0s
		if (iloop==0) t0s=tms ;
		Double_t Nseconds = 0.001*(tms - t0s) ;
		UShort_t ddi , hhi , mni, ssi , ndd, nhh, nmn , nss  ;
		ndd = TMath::Floor(Nseconds/86400) ;
		ddi = dd + ndd ;
		nhh = TMath::Floor((Nseconds - ndd*86400)/3600.) ;
		hhi = hh+nhh;
		nmn = TMath::Floor((Nseconds - ndd*86400-nhh*3600)/60) ;
		mni = mn+nmn;
		nss = TMath::Floor(Nseconds - ndd*86400-nhh*3600-nmn*60) ;
		ssi = ss + nss;
		if (ssi>60) { mni++ ; ssi=ssi-60; }
		if (mni>60) { hhi++ ; mni=mni-60; }
		if (hhi>24) { ddi++ ; hhi=hhi-24; }
		date.Set(yy, mm, ddi, hhi, mni, ssi) ;
		em->utc = date ;

		em->At = dt*1.e9 ; emh->At = em->At ;
		std::valarray<double> I_tot = vItotals[iloop];//vectorObjetosSimuladosAlmacenadosCorrelativamente.Get_Itot[iloop] ;
		for ( int i=0 ; i< em->Nt ; i++) {
			em->volt[i] = I_tot[i] ;
			em->time[i] = i*em->At ;

		}

		//Estimate polarity
		if ( TMath::Abs(TMath::MaxElement(em->Nt,em->volt)) > TMath::Abs(TMath::MinElement(em->Nt,em->volt)) ) Polarity++ ;
		else Polarity-- ;

		em->event=iRead ;

		//Now postprocess this entry (find out baseline, rtime and so on
		//TWaveform *wvi = new TWaveform( em ) ;
		TWaveform wvi = TWaveform( em ) ;

		if (iRead==0) tree->Branch("proc" , &wvi , 32000 , 0 );


		tree->Fill() ;
		iRead++   ;

		//delete wvi ;

		iactual++ ;

	}

	emh->Polarity = (Polarity>1) ? 1 : -1 ;
	tree->GetUserInfo()->Add( emh ) ;
	em->Ntevent=iRead ; //It does not go into the tree, only in the class!

	//delete emh ;


}


