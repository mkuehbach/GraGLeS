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
#ifndef __CURVAPPRX__
#define	 __CURVAPPRX__
#include <cmath>
/*!
 * \struct wghtd_imcurv_res
 * \brief Structure used to represent a two dimensional point
 *
 * The point represented by this structure has coordinates of type double. No operators
 * are overloaded for this structure.
 */

struct wghtd_imcurv_res {
	/*
	double imcurv_with_tjp;					//integral mean curvature including vertices at triple junctions
	double imcurv_without_tjp;				//integral mean curvature excluding vertices at triple junctions
	double effcapindspeed_with_tjp;			//weighted integral mean curvature with tjp
	double effcapindspeed_without_tjp;		//and without
	double effcapdrvforce_with_tjp;			//weighted integral mean curvature with tjp
	double effcapdrvforce_without_tjp;		//and without
	unsigned int cp_with_tjp;				//total number of vertices
	unsigned int cp_without_tjp;			//only virtual vertices not immediately adjacent to triple junctions
	unsigned int gID;
	bool status;
	bool pad1;
	bool pad2;
	bool pad3;
	//three padding bytes remain
	wghtd_imcurv_res(double _i1t, double _i0t, double _wi1t, double _wi0t, double _xi1t, double _xi0t, unsigned int _n1t, unsigned int _n0t, unsigned int _gid, bool _stat) : 
		imcurv_with_tjp(_i1t), imcurv_without_tjp(_i0t), effcapindspeed_with_tjp(_wi1t), effcapindspeed_without_tjp(_wi0t), 
			effcapdrvforce_with_tjp(_xi1t), effcapdrvforce_without_tjp(_xi0t), cp_with_tjp(_n1t), cp_without_tjp(_n0t), gID(_gid), 
			status(_stat), pad1(false), pad2(false), pad3(false) {}
	wghtd_imcurv_res() : imcurv_with_tjp(0.0), imcurv_without_tjp(0.0), effcapindspeed_with_tjp(0.0), effcapindspeed_without_tjp(0.0), 
		effcapdrvforce_with_tjp(0.0), effcapdrvforce_without_tjp(0.0), cp_with_tjp(0), cp_without_tjp(0), gID(0), 
		status(false), pad1(false), pad2(false), pad3(false) {}
	*/
	double imcurv_without_tjp;				//\kappa(s)
	double p_without_tjp;					//\gamma*\kappa(s)
	double v_without_tjp;					//m*\gamma*\kappa(s)
	unsigned int cp_with_tjp;				//total contour supporting points
	unsigned int cp_without_tjp;			//how many of these between triple junction sectors?
	unsigned int gID;
	unsigned int padding;
	wghtd_imcurv_res(double _i0t, double _p, double _v, unsigned int _nt, unsigned int _n0t, unsigned int _gid ) :
		imcurv_without_tjp(_i0t), p_without_tjp(_p), v_without_tjp(_v), cp_with_tjp(_nt),
			cp_without_tjp(_n0t), gID(_gid), padding(0) {};
	wghtd_imcurv_res() : 
		imcurv_without_tjp(std::numeric_limits<double>::max()), p_without_tjp(std::numeric_limits<double>::max()),
			v_without_tjp(std::numeric_limits<double>::max()), cp_with_tjp(0), cp_without_tjp(0), gID(0), padding(0) {};
	~wghtd_imcurv_res(){};
};

#endif	//__CURVAPPRX__
