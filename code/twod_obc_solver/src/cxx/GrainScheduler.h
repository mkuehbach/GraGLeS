#ifndef __GRAIN_SCHEDULER__
#define	__GRAIN_SCHEDULER__

#include <vector>
#include "Spoint.h"


using namespace std;

struct SPoint;

class GrainScheduler
{
public:
	GrainScheduler(){}
	virtual ~GrainScheduler(){}
	virtual void buildGrainWorkloads(vector<vector<SPoint>>&, int) = 0;
	virtual std::vector<unsigned int>&	getThreadWorkload(int threadID) = 0;
};

#endif
