/*
	IMM_MicrostructureGenerator
	A program to instantiate two-staged Poisson-Voronoi tessellation microstructures of 
	parent grains and their sub-grains with adjustable properties such as orientation, and dislocation density
	Copyright (C) 2016
	Christian Miessen (data structures), Markus Kühbach (physical metallurgy functionalities, PRNGs), 
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
	Settings::StatusHealthy = true;
	try {
		if (argc > 1)
			Settings::readXML(argv[1]);
		else
			Settings::readXML();
	} catch (exception& e) { 
		cout << "Unable to parse parameters file! Details:" << endl << e.what() << endl << "Simulation will now halt!" << endl; return 0;
	}
	if (Settings::StatusHealthy == false ) { cout << "Parameter set inconsistency!" << endl << "Simulation will now halt!" << endl; return 0; }

	microStructureHdl myHdl;

	//myHdl.DebugHDF5XDMF();
	//return 0;

	//myHdl.testprng( 100000, 1.0e+14, 5.0e+13 );
	myHdl.initEnvironment();
	myHdl.ReadAdditionalInputFiles();
	if (Settings::StatusHealthy == false) { cout << "Simulation will now halt!" << endl; return 0; }

	myHdl.GeneratePolycrystallineStructureOfGrains();
	myHdl.DistributeGrainOriAndSEE();
	myHdl.GenerateSubgrainStructureInEachGrain();
	myHdl.DistributeSubgrainOrientations();
	myHdl.DistributeSubgrainSEE();

	myHdl.RehashGrainIDs();
	myHdl.BreakPeriodicity();
	//from now on all sub-grains and grains have a contiguous numbering

	myHdl.SaveDataGraGeLeS();
	myHdl.SaveDataDAMASK();
	myHdl.SaveDetailedDiagnosticsASCII();
	/*	deprecated I/O functions
		myHdl.SaveDetailedDiagnosticsBINARY();
		myHdl.SaveParenthood();
		myHdl.DebugHDF5();
	*/
	myHdl.Plot3DVolume(); //#output too small?
//##MK	myHdl.SaveHDF5();
	myHdl.PlotIPF2DSection();
	
	myHdl.ReportProfile();
	cout << "Program finished successfully!";
	return 0;
}

