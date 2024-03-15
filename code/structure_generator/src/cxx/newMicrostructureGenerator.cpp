/*
	IMM_MicrostructureGenerator
	A program to instantiate two-staged Poisson-Voronoi tessellation microstructures of 
	parent grains and their sub-grains with adjustable properties such as orientation, and dislocation density
	Copyright (C) 2016
	Christian Miessen (data structures), Markus K�hbach (physical metallurgy functionalities, PRNGs), 
	Nikola Velinov (data structures), Luis Antonio Barrales-Mora (PRNGs, Math), Jonathan Nierengarten

	The work was funded by the DFG Reinhart-Koselleck project GO 335/44-1

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

#include <iostream>
#include "microStructureHdl.h"
#include "Settings.h"
using namespace std;

int main(int argc, char *argv[]) {
	if ( argc != 3 ) {
		cout << "Command line call has to be <<app>> <<simid>> <<configfile>>" << "\n";
		return 0;
	}
	double gtic = omp_get_wtime();
	Settings::SimulationId = std::stoul(argv[1]);
	Settings::ConfigFileName = argv[2];
	Settings::ResultsFileName = "StructureGenerator.Results.SimID." 
		+ to_string(Settings::SimulationId) + ".nxs";
	
	Settings::ReadXmlConfig(Settings::ConfigFileName);
	if (Settings::StatusHealthy == false ) {
		cout << "Loading configuration " << Settings::ConfigFileName << " failed or is invalid!" << "\n";
		return 0;
	}

	HdfFiveSeqHdl h5w = HdfFiveSeqHdl( Settings::ResultsFileName );
	if( h5w.nexus_create() != MYHDF5_SUCCESS ) { 
		cout << "Creating results NeXus/HDF5 file " << Settings::ResultsFileName << " failed!" << "\n";
		return 0;
	}

	microStructureHdl myHdl = microStructureHdl();
	myHdl.InitEnvironment();
	myHdl.ReadAdditionalInputFiles();
	myHdl.ReportConfig();
	if (Settings::StatusHealthy == false) {
		cout << "Configuring the tool failed!" << "\n"; return 0;
	}
	myHdl.GenerateGrainStructureAsParentGrains();
	myHdl.DistributeGrainOriAndSee();
	myHdl.GenerateSubgrainStructureInEachGrain();
	myHdl.DistributeSubgrainOrientations();
	myHdl.DistributeSubgrainSee();
	myHdl.RehashGrainIds();
	//myHdl.BreakPeriodicity();
	//from now on all sub-grains and grains have a contiguous numbering
	myHdl.SaveNeXus();

	myHdl.ReportProfile();
	cout << "structure_generator finished successfully in " << (omp_get_wtime() - gtic) << " seconds" << endl;
	return 0;
}

