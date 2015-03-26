#include <vector>
#include "ofVec3f.h"

class crvPt;

class polyLine {
public:
	vector<crvPt> pts;
	float length;
	bool isClosed;

	float getLength();
};

class crvPt : public ofVec3f {
public:
	float t;
	polyLine * crv;
	vector<crvPt *> parents;

	crvPt() {

	}

	crvPt(const ofVec3f & v) {
		x = v.x;
		y = v.y;
		z = v.z;
	}

	void operator=(const ofVec3f &v) {
		x = v.x;
		y = v.y;
		z = v.z;
	}

};

