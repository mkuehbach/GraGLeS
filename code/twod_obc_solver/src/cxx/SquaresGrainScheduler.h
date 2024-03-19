#ifndef __SQUARESGRAINSCHEDULER_H__
#define __SQUARESGRAINSCHEDULER_H__

#include "GrainScheduler.h"
using namespace std;

struct SPoint;

class SquaresGrainScheduler: public GrainScheduler {
public:
	SquaresGrainScheduler(int numberOfThreads, int totalNumberOfGrains, int startIndex = 1);
	~SquaresGrainScheduler();

	void buildGrainWorkloads(vector<vector<SPoint>>& contours, int n_gridpoints);
	void buildSubWorkloads(int offsetThreadID, vector<unsigned int>& grainlist,
			int NumberOfThreads);
	SPoint find_center(vector<SPoint>& contour, int n_gridpoints);
	vector<unsigned int>& getThreadWorkload(int threadID);

private:
	int m_totalNumberOfThreads;
	int m_totalNumberOfGrains;
	int m_startIndex;
	vector<vector<unsigned int> > m_threadWorkVectors;
};

#endif
