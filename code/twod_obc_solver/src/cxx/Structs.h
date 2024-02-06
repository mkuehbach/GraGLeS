struct VolEvolution {
	double dA;
	int nVertex;

	VolEvolution(double da, int nr) :
			dA(da), nVertex(nr) {
	}
};

struct Face {
	double length;
	unsigned int grainA;
	unsigned int grainB;
	Face(double _length, unsigned int _grainA, unsigned int _grainB) :
			length(_length), grainA(_grainA), grainB(_grainB) {
	}
};

struct TextureData {
	double volume;
	double perimeter;
	double GBEnergy;
	double BulkEnergy;
	double phi1;
	double PHI;
	double phi2;
	double x;
	double y;
	unsigned int id;
	unsigned int NeighbourCount;
	unsigned int intersectsBoundaryGrain;
	TextureData(unsigned int _id, unsigned int _NeighbourCount,
			bool _intersectsBoundaryGrain, double _volume, double _surfaceArea,
			double _GBEnergy, double _BulkEnergy, Quaternion *ori, double _x,
			double _y) :
			id(_id), NeighbourCount(_NeighbourCount), intersectsBoundaryGrain(
					_intersectsBoundaryGrain), volume(_volume), perimeter(
					_surfaceArea), GBEnergy(_GBEnergy), BulkEnergy(_BulkEnergy), x(
					_x), y(_y) {
		double* newori = ori->quaternion2Euler();
		phi1 = newori[0];
		PHI = newori[1];
		phi2 = newori[2];
		delete[] newori;

	}
	~TextureData() {
	}
};

//SI unit scaling
#define METER2MICRON(x)				((x)*(1.0e6))
#define MICRON2METER(x)				((x)*(1.0e-6))

/*
struct wghtd_imcurv_res {
	double imcurv_without_tjp;
	unsigned int cp_without_tjp;
	unsigned int gID;
	wghtd_imcurv_res(double _i0t, unsigned int _n0t, unsigned int _gid ) :
		imcurv_without_tjp(_i0t), cp_without_tjp(_n0t), gID(_gid) {};
	wghtd_imcurv_res() : 
		imcurv_without_tjp(std::numeric_limits<double>::max()), cp_without_tjp(0), gID(0) {};
	~wghtd_imcurv_res(){};
};
*/

struct GBContourPoint {
	double x;
	double y;
	double energy;
	double mobility;
	unsigned int myid;
	unsigned int nborid;
	GBContourPoint( double _x, double _y, double _en, double _mob,
		unsigned int _myid, unsigned int _nborid ) : 
		x(_x), y(_y), energy(_en), mobility(_mob), myid(_myid), nborid(_nborid) {};
	~GBContourPoint() {}; //
};

struct GBJunctionPoint {
	double x;
	double y;
	unsigned int myid;
	unsigned int jtype;
	GBJunctionPoint( double _x, double _y, unsigned int _myid, unsigned int _jtype ) : 
		x(_x), y(_y), myid(_myid), jtype(_jtype) {};
	~GBJunctionPoint() {};
};
