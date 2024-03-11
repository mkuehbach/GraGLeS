/*
 GraGLeS 2D A grain growth simulation utilizing level set approaches
 Copyright (C) 2015  Christian Miessen, Nikola Velinov

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "ggLS.h"
#include "Settings.h"
#include <sys/time.h>
using namespace voro;
using namespace std;



#include <stdlib.h>
#include <stdio.h>


bool debug_hdf5()
{
	HdfFiveSeqHdl h5w = HdfFiveSeqHdl( Settings::ResultsFileName );
	if( h5w.nexus_create() != MYHDF5_SUCCESS ) { return false; }

	ioAttributes anno = ioAttributes();
	string grpnm = "";
	string dsnm = "";

	unsigned int entry_id = 1;
	grpnm = "/entry" + to_string(entry_id);

	stringstream sstr;
	sstr << "https://github.com/mkuehbach/GraGLeS/<<ADD COMMID ID>>/NXms_gragles_results";
	cout << sstr.str() << "\n";
	anno = ioAttributes();
	anno.add( "version", string(sstr.str()) );
	anno.add( "NX_class", string("NXentry") );
	if ( h5w.nexus_write_group( grpnm, anno ) != MYHDF5_SUCCESS ) { return false; }

	dsnm = grpnm + "/definition";
	string appdef_name = "NXms_gragles_results";
	anno = ioAttributes();
	if ( h5w.nexus_write( dsnm, appdef_name, anno ) != MYHDF5_SUCCESS ) { return false; }

	dsnm = grpnm + "/program";
	string program_name = "twod_obc_solver";
	anno = ioAttributes();
	anno.add( "version", string("<<ADD PARSING OF xstr(GITSHA)>>") );
	anno.add( "preprocessor_date", string(__DATE__) );
	anno.add( "preprocessor_time", string(__TIME__) );
	if ( h5w.nexus_write( dsnm, program_name, anno ) != MYHDF5_SUCCESS ) { return false; }

	dsnm = grpnm + "/analysis_identifier";
	unsigned int identifier = Settings::SimID;
	anno = ioAttributes();
	//anno.add( "unit", string("NX_UNITLESS") );
	if ( h5w.nexus_write( dsnm, identifier, anno ) != MYHDF5_SUCCESS ) { return false; }

	dsnm = grpnm + "/start_time";
	string start_time_stamp = "<<ADD TIME PARSING via utils/src/cxx>>";
	anno = ioAttributes();
	//MK::IF THE VALUE to WRTIE IS A STRING IT MUST NOT BE A CONST STRING
	if ( h5w.nexus_write( dsnm, start_time_stamp, anno ) != MYHDF5_SUCCESS ) { return false; }

	dsnm = grpnm + "/config_filename";
	string config_file_sha256 = "<<ADD SHA256( ConfigShared::ConfigurationFile )>>";
	string config_file_name = "<<ADD ConfigShared::ConfigurationFile>>";
	cout << "config_file_sha256 " << config_file_sha256 << "\n";
	cout << "config_file_name " << config_file_name << "\n";
	anno = ioAttributes();
	anno.add( "version", config_file_sha256 );
	anno.add( "comment", string("SHA256 checksum"));
	if ( h5w.nexus_write( dsnm, config_file_name, anno ) != MYHDF5_SUCCESS ) { return false; }

	grpnm = "/entry" + to_string(entry_id) + "/coordinate_system_set";
	anno = ioAttributes();
	anno.add( "NX_class", string("NXcoordinate_system_set") );
	if ( h5w.nexus_write_group( grpnm, anno ) != MYHDF5_SUCCESS ) { return false; }

	grpnm = "/entry" + to_string(entry_id) + "/coordinate_system_set/gragles";
	anno = ioAttributes();
	anno.add( "NX_class", string("NXtransformations") );
	if ( h5w.nexus_write_group( grpnm, anno ) != MYHDF5_SUCCESS ) { return false; }

	vector<string> axis_name = { "/x", "/y", "/z" };
	for( size_t i = 0; i < 3; i++ ) {
		dsnm = grpnm + axis_name[i];
		vector<double> real = vector<double>( 3, 0. );
		real[i] = 1.;
		anno = ioAttributes();
		anno.add( "depends_on", string(".") );
		anno.add( "offset", string("{0., 0., 0.}, storing 1d array as an attribute value not yet implemented") ); //##MK::not yet implemented
		anno.add( "offset_units", string("nm") );
		if ( h5w.nexus_write(
				dsnm,
				io_info({3, 1},
						{3, 1}),
				real,
				anno) != MYHDF5_SUCCESS ) { return false; }
	}
	//##MK::add dim_alias names
	return true;
}



int main(int argc, char *argv[]) {

	Settings::ResultsFileName = "Twod.Obc.Solver.Results.SimID." + to_string(Settings::SimID) + ".nxs";
	return (int) debug_hdf5();

	if (argc > 1)
		Settings::initializeParameters(argv[1]);
	else
		Settings::initializeParameters();

	grainhdl* my_sim = NULL;

	my_sim = new grainhdl();

	my_sim->setResearchAdjustments(Settings::ResearchProject);
	timeval time_start, time_end;
	gettimeofday(&time_start, NULL);
	my_sim->setSimulationParameter();
	gettimeofday(&time_end, NULL);
	double elapsed_secs = (time_end.tv_sec - time_start.tv_sec)
			+ (time_end.tv_usec - time_start.tv_usec) / 1000000.0;
	cout << "elapsed secs for Initializing network (Read):" << elapsed_secs << endl << endl;

	if (Settings::MicrostructureGenMode == E_GENERATE_WITH_VORONOY )
		my_sim->save_Full_Microstructure_for_Restart();

	gettimeofday(&time_start, NULL);
	my_sim->run_sim();

	gettimeofday(&time_end, NULL);
	elapsed_secs = (time_end.tv_sec - time_start.tv_sec) + (time_end.tv_usec
			- time_start.tv_usec) / 1000000.0;
	cout << "elapsed secs for main loop:" << elapsed_secs << endl;

	my_sim->save_NrGrainsStats();

	my_sim->clear_mem();

	delete my_sim;
	return 0;
}
