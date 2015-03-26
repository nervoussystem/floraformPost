#include "polyLine.h"

float polyLine::getLength() {
	length = 0;
	for(int i=0;i<pts.size()-1;++i) {
		pts[i].t = length;
		length += pts[i].distance(pts[i+1]);
	}
	if(isClosed) {
		length += pts[0].distance(pts[pts.size()-1]);
		pts.back().t = length;
	}
	return length;
}