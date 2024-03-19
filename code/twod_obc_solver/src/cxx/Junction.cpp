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
#include "Junction.h"
#include "LsBox.h"

grainhdl* GrainJunction::handler = NULL;

double GrainJunction::getWeight(LSbox* me)
{
	unsigned int neighbours[3] = { 0xFFFFFF, 0xFFFFFF, 0xFFFFFF };
	for (int i = 0, j = 0; i < junction_type; i++) {
		if (me->getID() != grains[i]) {
			neighbours[j] = grains[i];
			j++;
		}
	}
	double averageMobility = 0.;
	double sigma = 1.; //TODO:: was undefined
	double gamma[3] = {0., 0., 0.};
	double gamma_hagb = 1.;
	double theta_ref = DEG2RAD(15.);
	double theta_mis = 1.; //TODO::was undefined

	if ( Settings::MicrostructureGenMode == E_READ_VOXELIZED_MICROSTRUCTURE ) {
		//Settings::ResearchMode == 0)
		theta_mis = me->computeMisorientation(neighbours[0]);
		//characteristics& myNeighbor = me->m_grainBoundary.getDirectNeighbourCaracteristic(handler->getGrainByID(neighbours[0]));
		//averageMobility += myNeighbor.mobility;
		averageMobility += me->GbMobilityModel(theta_mis, handler->getGrainByID(neighbours[0]));
		gamma[0] = me->GbEnergyModel(theta_mis, handler->getGrainByID(neighbours[0]));
		
		if (Settings::IdentifyTwins) {
			if (me->MisoriToTwinBoundary(handler->getGrainByID(neighbours[0])) < DEG2RAD(8.66025)) {
				//if the GB is TWIN Boundary set the Energy to the lowest one allowed : theta_mis = 1 degree
				gamma[0] = me->GbEnergyModel(DEG2RAD(1.), handler->getGrainByID(neighbours[0]));
			}
		}
		theta_mis = handler->getGrainByID(neighbours[0])->computeMisorientation(neighbours[1]);
		characteristics & myNeighbor = handler->getGrainByID(
				neighbours[0])->m_grainBoundary.getDirectNeighbourCaracteristic(handler->getGrainByID(neighbours[1]));
		averageMobility += myNeighbor.mobility;
		gamma[1] = handler->getGrainByID(neighbours[0])->GbEnergyModel(theta_mis, handler->getGrainByID(neighbours[1]));
		if (Settings::IdentifyTwins) {
			if (handler->getGrainByID(neighbours[0])->MisoriToTwinBoundary(handler->getGrainByID(neighbours[1])) < DEG2RAD(8.66025)) {
				//if the GB is TWIN Boundary set the Energy to the lowest one allowed : theta_mis = 1 degree
				gamma[1] = handler->getGrainByID(neighbours[0])->GbEnergyModel(DEG2RAD(1.), handler->getGrainByID(neighbours[1]));
			}
		//myNeighbor = me->m_grainBoundary.getDirectNeighbourCaracteristic(handler->getGrainByID(neighbours[1]));
		//averageMobility += myNeighbor.mobility;
		theta_mis = me->computeMisorientation(neighbours[1]);
		averageMobility += me->GbMobilityModel(theta_mis, handler->getGrainByID(neighbours[1]));
		gamma[2] = me->GbEnergyModel(theta_mis, handler->getGrainByID(neighbours[1]));
		if (Settings::IdentifyTwins) {
			if (me->MisoriToTwinBoundary(handler->getGrainByID(neighbours[1])) < DEG2RAD(8.66025)) {
				//if the GB is TWIN Boundary set the Energy to the lowest one allowed : theta_mis = 1 degree
				gamma[2] = me->GbEnergyModel(DEG2RAD(1.), handler->getGrainByID(neighbours[1]));
			}
		}
	}
	else {
		cout << "error in junction.cpp - no rule found" << endl;
		exit(2);
	}

	// find the asociated weight
	sigma = gamma[0] - gamma[1] + gamma[2];
	/*
	if (Settings::IsIsotropicNetwork == 0) {
		if (Settings::UseMobilityModel > 0 && Settings::TriplePointDrag > 0) {
			averageMobility /= 3;
			double ds = handler->get_ds();
			double drag = 1 / ((1 / (ds * Settings::TriplePointDrag)) + 1
					/ averageMobility);
			sigma *= drag;
		}
	}
	*/

	if (sigma < 0.01) {
		//cout << "negative sigma " << endl;
		sigma = 0.01;
	}
	return sigma;
}
