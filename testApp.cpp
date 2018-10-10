#include "testApp.h"

float dt = 1;
float diffusionRateA = .05;
float diffusionRateB = .025;
float reactionRate = 2.0;
float F = 0.029;//.03
float k = 0.06;//.063
float offset = 3.0;//.5;
float threshold = .45;
//it doesn't use this, it uses the settings.xml file
string file = "ringAdjust3_doneJesseEdit.obj";
bool paused = true;

extern vector<float> distances2;

float maximumDist = 0; //10
float geodesicSpacing = 1.25;//2
float holeSpacing = 2.5;//2.5

float colorDist = 30;
float fallOff = 5;

float thickness = 1.3;

float edgeThickness = 0.5;

float centerThickness = 0.8;
float maxThickness = 0.9;
//only these ones are used
float centerThicknessThin = 0.8;//thickness close to the edge
float maxThicknessThin = 1.5;//thickness far from the edge
float centerThicknessThick = centerThicknessThin;
float maxThicknessThick = maxThicknessThin;
/*
float rimT[2] = {.5,.5};
float edgeT[2] = {.85,.95};
float midBodyT[2] = {0.9,1.7};
float bodyT[2] = {1.0,1.8};
float midBodyDist = 2.5;
float edgeDist = 1.5;
float maxY = 4.5;
float minY = 7;
*/
//florescene engagement ring settings
float rimT[2] = { .5,1.2 };
float edgeT[2] = { .75,1.5 };
float midBodyT[2] = { .8,2. };
float bodyT[2] = { .8,2. };
float midBodyDist = 2;
float edgeDist = .7;
float minY = 5;
float maxY = 0;

float rimThick = rimT[0];
float edgeThick = edgeT[0];
float bodyThick = bodyT[0];
float midBodyThick = midBodyT[0];





bool limitEdge = false;

bool doHoles = false;
//cuff ellipse
//necklace? .887, 30.23
//hoop 0,3.887
ofVec2f centerPt(0.0,0.0);//(0,3.887);//(-.69,37.66);//-13.78);
float radX = 9.186;// 32;//8.6;//9.5;//13;//67.157;//62.566;//28.1;//35.485;//9.5;
float radY = 9.186;// 24.5;//8.6;//9.5;//13;//83.66;//98.76;//35.485;//9.5;

hemesh hmesh;
vector<float> a;
vector<float> b;
ofVboMesh mesh;

SparseMatrix A;
SparseMatrix laplacianA;
SparseMatrix laplacianB;
SparseSolver solverA;
//Eigen::BiCGSTAB<SparseMatrix>
SparseSolver solverB;

bool doBorder = true;

typedef pair<ofVec3f, ofVec3f> segment;

class crvVec3 : public ofVec3f {
public:
	float t;

	void operator=(const ofVec3f &v) {
		x = v.x;
		y = v.y;
		z = v.z;
	}
};

vector< vector<polyLine *> > geodesics;

vector< vector< vector< crvPt > * > * > holePts;
vector< vector< crvPt> > holePts2;
vector<ofVec3f> sLines;

unsigned int stepCount = 0;

bool isMovie = false;

void saveLines();

bool processing = false;
int fileIndex = 0;
vector<ofFile> files;

//--------------------------------------------------------------
void testApp::setup(){
	loadSettings();
	loadMesh(file);
	hmesh = hemeshFromOfMesh(mesh);

	for(int i=0;i<hmesh.vertices.size();++i) {
		hevertex * v = hmesh.vertices[i];
		if(v->boundary) {

			ofVec3f pos = mesh.getVertex(i);
			float xTerm= pow( (pos.x - centerPt.x) , 2 ) / (radX*radX);
			float yTerm = pow( (pos.y - centerPt.y) , 2 ) / (radY*radY);
			v->tag = true;
			if(limitEdge && xTerm + yTerm <= 1){
				v->tag = false;
			} else {
				v->tag = true;
			}
		}
	}
	cout << "loaded" << endl;
	setupSolver(hmesh);
	cout << "solver setup" << endl;
	distances.resize(hmesh.vertices.size());
	geodesicDistance(hmesh,distances);
	distances2.resize(hmesh.vertices.size());
	geodesicDistance(hmesh,distances2, true);
	cout << "distance calculated" << endl;
	//exportDistObj();
	maximumDist = 0;
	for(int i=0;i<distances.size();++i) {
		maximumDist  = max(maximumDist, distances[i]);
	}
	cout << "max dist " << maximumDist << endl;
	//adaptiveSubdivision();
	//setupSolver(hmesh);
	//computeSecondNeighbors();
	//a.resize(hmesh.vertices.size());
	//b.resize(hmesh.vertices.size());
	//distances.resize(hmesh.vertices.size());
	//geodesicDistance(hmesh,distances);
	//assembleMatrix();
	
	//float maxDist = computeDistance(mesh,edges,distances);
	//initializeFMM(mesh,edges);
	//initValues();
	//computeDistance(mesh,edges,distances);
	//enable color

	//mesh.enableColors();
	
	ofVec3f minPt,maxPt;
	maxPt = minPt = mesh.getVertex(0);
	for(int i=1;i<mesh.getNumVertices();++i) {
		ofVec3f pt = mesh.getVertex(i);
		minPt.set(min(minPt.x,pt.x),min(minPt.y,pt.y),min(minPt.z,pt.z));
		maxPt.set(max(maxPt.x,pt.x),max(maxPt.y,pt.y),max(maxPt.z,pt.z));
	}
	mCam.setupPerspective(true,60.0,0.001,1000.0);
	mCam.setTarget((minPt+maxPt)*0.5);
	mCam.setDistance(3.0*max(max(maxPt.x-minPt.x,maxPt.y-minPt.y),maxPt.z-minPt.z));
	cout << "normals" << endl;
	updateMeshNormals(mesh);
	
	//calculate holes
	if(doHoles){
		getGeodesics();
		getHolePts();
	}
	cout << "setup done" << endl;
	//markGrowthPts();
	//fixColors();
	//end calculate holes
	
//	saveLines();
//	optimizePts();
//	optimizePts();
//	exportDistObj();
//	exportColor();
//	updateSolver(dt);
//	geodesicDistance(distances);
}

void testApp::markGrowthPts() {
	//check for isolated white spots
	for(int i=0;i<hmesh.vertices.size();++i) {
		ofColor c = mesh.getColor(i);
		float dist = distances[i];
		hevertex * v = hmesh.vertices[i];
		if(c == ofColor::white) {
			bool isSub = true;
			if(dist > maximumDist*0.56) {
				//check neighbors
				hedge * e = v->edge;
				hedge * start = e;
				ofVec3f newColor;
				int neighs = 0;
				do {
					e = e->next;
					int index = e->vertex->index;
					ofColor c2 = mesh.getColor(index);
					if(c2 == ofColor::white) {
						isSub = false;
					}
					newColor += ofVec3f(c2.r,c2.g,c2.b);
					neighs++;
					e = e->pair;
				} while(e != start);
				if(isSub) {
					newColor /= neighs;
					mesh.setColor(i,ofColor(newColor.x,newColor.y,newColor.z));
				} else {
					mesh.setColor(i,ofColor::red);
				}
			}
		}

	}
	for(int i=0;i<hmesh.vertices.size();++i) {
		ofColor c = mesh.getColor(i);
		float dist = distances[i];
		hevertex * v = hmesh.vertices[i];
		v->tag = false;

		if(dist > maximumDist*0.56) {
			if(c == ofColor::red) {
				v->tag = true;
			}
		}
	}
	geodesicDistance(hmesh,distances);

}

void testApp::process() {
	hmesh = hemeshFromOfMesh(mesh);

	for(int i=0;i<hmesh.vertices.size();++i) {
		hevertex * v = hmesh.vertices[i];
		if(v->boundary) {

			//ofVec3f pos = mesh.getVertex(i);
			//float xTerm= pow( (pos.x - centerPt.x) , 2 ) / (radX*radX);
			//float yTerm = pow( (pos.y - centerPt.y) , 2 ) / (radY*radY);
			//v->tag = true;
			//if(limitEdge && xTerm + yTerm <= 1){
			//	v->tag = false;
			//} else {
				v->tag = true;
				//v->tag = false;
			//}
		}
	}
	setupSolver(hmesh);
	distances.resize(hmesh.vertices.size());
	geodesicDistance(hmesh,distances);
	//distances2.resize(hmesh.vertices.size());
	//geodesicDistance(hmesh,distances2, true);
	//exportDistObj();
	maximumDist = 0;
	for(int i=0;i<distances.size();++i) {
		maximumDist  = max(maximumDist, distances[i]);
	}
	markGrowthPts();
	fixColors();
}

void testApp::fixColors() {
	for(int i=0;i<hmesh.vertices.size();++i) {
		float t = ofClamp(1.0-distances[i]/colorDist,0,1);
		t = pow(t,fallOff);
		mesh.setColor(i,ofColor(255,255-t*255.0,255- t*255.0));
	}
}

float crvLen(vector<crvVec3> & crv) {
	return crv[crv.size()-1].t;
}

float crvDistance(crvVec3 & v1, crvVec3 & v2, float len) {
	float d1 = v2.t-v1.t;
	float d2 = len-abs(d1);
	if(d1 > 0) {
		d2 *= -1;
	}
	if(abs(d1) < abs(d2)) {
		return d1;
	} else {
		return d2;
	}
	//return abs(v1.t-v2.t);
}

float crvDistance(crvPt & v1, crvPt & v2) {
	if(v1.crv->isClosed) {
		float d1 = v2.t-v1.t;
		float d2 = v1.crv->length-abs(d1);
		if(d1 > 0) {
			d2 *= -1;
		}
		if(abs(d1) < abs(d2)) {
			return d1;
		} else {
			return d2;
		}
	} else {
		return v2.t-v1.t;
	}
	//return abs(v1.t-v2.t);
}

crvVec3 findPt(vector<crvVec3> &crv, float t, int startI, int endI) {
	int midI = (endI+startI)/2;

	crvVec3 pt = crv[midI*2];
	crvVec3 pt2 = crv[midI*2+1];
	if(t < pt.t) {
		return findPt(crv, t, startI, midI);
	} else if(t > pt2.t) {
		return findPt(crv, t, midI, endI);
	} else {
		crvVec3 out = pt;
		out += (pt2-pt)*(t-pt.t);
		out.t = t;
		return out;
	}
}

crvVec3 getPt(vector<crvVec3> &crv, float t) {
	float cLen = crvLen(crv);
	if(t < 0) t += cLen;
	else if(t >= cLen) t -= cLen;
	int startI = 0;
	int endI = crv.size()/2;
	int midI = (endI+startI)/2;

	return findPt(crv,t,startI,endI);
}

crvPt findPt(polyLine * crv, float t, int startI, int endI) {
	if(endI-startI == 1) {
		crvPt pt = crv->pts[startI];
		crvPt pt2 = crv->pts[(endI)%crv->pts.size()];

		crvPt out = pt;
		out += (pt2-pt)*(t-pt.t);
		out.t = t;
		out.crv = crv;
		return out;
	}
	int midI = (endI+startI)/2;


	crvPt pt = crv->pts[midI];
	//crvPt pt2 = crv->pts[(midI+1)%crv->pts.size()];
	if(t < pt.t) {
		return findPt(crv, t, startI, midI);
	} else if(t > pt.t) { //pt2
		return findPt(crv, t, midI, endI);
	} else {
		return pt;
	}
	/*else {
		crvPt out = pt;
		out += (pt2-pt)*(t-pt.t);
		out.t = t;
		out.crv = crv;
		return out;
	}
	*/
}

crvPt getPt(polyLine * crv, float t) {
	float cLen = crv->length;
	if(crv->isClosed) {
		if(t < 0) t += cLen;
		else if(t >= cLen) t -= cLen;
	} else {
		if(t < 0) t = 0;
		else if(t > cLen) t = cLen;
	}
	int startI = 0;
	int endI = crv->pts.size()-1;
	if(crv->isClosed) {
		endI++;
	}
	int midI = (endI+startI)/2;

	return findPt(crv,t,startI,endI);
}

/*
float projectPtToLine(const ofVec3f &pt, const crvVec3 &l1, const crvVec3 &l2, crvVec3 & out) {
	ofVec3f dir = l2-l1;
	ofVec3f dir2 = pt-l1;
	float len = dir.length();
	float t = dir.dot(dir2)/len;
	t = max(0.f,min(t,1.0f));
	out = l1+dir*t;
	out.t = l1.t+t*len;
	return out.distanceSquared(pt);
}
*/

float projectPtToLine(const ofVec3f &pt, const crvPt &l1, const crvPt &l2, crvPt & out) {
	ofVec3f dir = l2-l1;
	ofVec3f dir2 = pt-l1;
	float len = dir.length();
	float t = dir.dot(dir2)/len;
	t = ofClamp(t,0.0f,len);
	out = l1+dir*t/len;
	out.t = l1.t+t;
	return out.distanceSquared(pt);
}

void intersect_parents(	vector<crvPt *> & parents1, vector<crvPt *> & parents2, vector<crvPt *> & out) {
	for(auto p : parents1) {
		if(find(parents2.begin(), parents2.end(), p) != parents2.end()) {
			out.push_back(p);
		}
	}
}

/*
crvVec3 projectPt(const ofVec3f &pt, vector<crvVec3> & crv) {
	float minDist = 9e9;
	crvVec3 closestPt;
	crvVec3 outPt;
	//find closest pt
	for(int i=0;i<crv.size();i+=2) {
		float dist = projectPtToLine(pt, crv[i], crv[i+1], outPt);
		if(dist < minDist) {
			minDist = dist;
			closestPt = outPt;
		}
	}
	return closestPt;
}
*/

crvPt projectPt(const ofVec3f &pt, polyLine * crv) {
	float minDist = 9e9;
	crvPt closestPt;
	crvPt outPt;
	//find closest pt
	int endI = crv->pts.size();
	if(!crv->isClosed) {
		endI--;
	}
	for(int i=0;i<endI;i++) {
		float dist = projectPtToLine(pt, crv->pts[i], crv->pts[(i+1)%crv->pts.size()], outPt);
		if(dist < minDist) {
			minDist = dist;
			closestPt = outPt;
		}
	}
	closestPt.parents.clear();
	closestPt.crv = crv;
	return closestPt;
}

crvPt projectPt(const ofVec3f &pt, vector<polyLine *> crvs) {
	float minDist = 9e9;
	crvPt closestPt;
	
	for(int i=0;i<crvs.size();++i) {
		crvPt p = projectPt(pt,crvs[i]);
		float d = p.squareDistance(pt);
		
		if( d < minDist ) {
			
			minDist = d;
			closestPt = p;
		}
	}
	return closestPt;
}

crvPt mergePts(crvPt & pt1, crvPt & pt2) {
	polyLine * crv = pt1.crv;

	float t = pt1.t;
	t += crvDistance(pt1,pt2)*0.5;
	crvPt newPt = getPt(crv, t);
	newPt.parents.insert(newPt.parents.end(), pt1.parents.begin(), pt1.parents.end());
	newPt.parents.insert(newPt.parents.end(), pt2.parents.begin(), pt2.parents.end());
	return newPt;

}

bool crvSortComp(crvPt &p1, crvPt &p2) {
	return p1.t < p2.t;
}

void testApp::getHolePts() {
	vector<polyLine *> crvs = geodesics[0];
	
	vector< vector< crvPt > * > * container;
	container = new vector< vector<crvPt> * >();

	for(auto outer : crvs) {
		vector<crvPt> * pts = new vector<crvPt>();
		//polyLine * outer = crvs[0];

		int numPts = (int) (outer->length/holeSpacing);
		float spacing = outer->length/numPts;
		int currCrv = 0;
		crvPt currPt = outer->pts[0];
		currPt.t = 0;
		currPt.crv = outer;
		pts->push_back(currPt);
		float currDist = 0;
		for(int i=0;i<numPts-1;++i) {
			while(currDist < spacing) {
				float lineDist = currPt.distance(outer->pts[currCrv+1]);
				if(currDist + lineDist > spacing) {
					ofVec3f dir = outer->pts[currCrv+1]-currPt;
					dir *= (spacing-currDist)/lineDist;
					currPt += dir;
					currPt.t = (i+1)*spacing;
					pts->push_back(currPt);
					currDist = 0;
					break;
				}
				currDist += lineDist;
				currPt = outer->pts[currCrv+1];
				currCrv += 2;
			}
		}
		if(!outer->isClosed) {
			currPt = outer->pts.back();
			currPt.crv = outer;
			currPt.t = outer->length;
			pts->push_back(currPt);
		}
		container->push_back(pts);
	}
	holePts.push_back(container);

	cout << "first done " << container->front()->size() << endl;
	crvPt newPt, newPt2;
	for(int i=1;i<geodesics.size()-1;i+=2) {
		vector<polyLine *> &crvs = geodesics[i];
		vector<polyLine *> &crvs2 = geodesics[i+1];
		if(crvs2.size() == 0) break;
		vector<crvPt> * pts2 = new vector<crvPt>();
		vector<crvPt> * pts = new vector<crvPt>();
		for(auto crv : *container) {
			pts->insert(pts->end(), crv->begin(), crv->end());
		}
		//container.clear();
		container = new vector<vector<crvPt> *>();
		container->resize(crvs2.size());
		for(int j=0;j<container->size();++j) {
			container->at(j) = new vector<crvPt>();
		}
		cout << "projecting points" << endl;
		for(int j=0;j<pts->size();++j) {
			newPt = projectPt(pts->at(j), crvs);
			if(newPt.distanceSquared(pts->at(j)) < geodesicSpacing*geodesicSpacing*2) {
				newPt2 = projectPt(newPt, crvs2);
				if(newPt2.distanceSquared(newPt) < geodesicSpacing*geodesicSpacing*2) {
					newPt2.parents.push_back(new crvPt(pts->at(j)));
					//pts2.push_back(newPt2);
					for(int k=0;k<crvs2.size();++k) {
						if(newPt2.crv == crvs2[k]) {
							container->at(k)->push_back(newPt2);
						}
					}
				}
			}
		}
		cout << "cont " << container->size() << " " << container->front()->size() << endl;

		//warning memory leak for pts
		//sort
		for(auto &crv : *container) {
			//vector<crvPt> &crv =  container[k];
			if(crv->size() > 2) {
				sort(crv->begin(),crv->end(),crvSortComp);
				pts2->clear();
				//get dual

				//need to check for closed curve
				polyLine * line = crv->at(0).crv;
				int end = crv->size();
				if(!line->isClosed) {
					end--;
				}
				for(int j=0;j<end;++j) {
					crvPt mPt = mergePts(crv->at(j),crv->at((j+1)%crv->size()));
					pts2->push_back(mPt);
				}
			
				//merge
				//pts = new vector<crvPt>();
				pts->clear();

				if(!line->isClosed) {
					crvPt mPt = getPt(line,0);
					//what about parent?
					pts->push_back(mPt);
				}
				bool merged = false;
				crvPt firstMerge;

				end = pts2->size();
				if(!line->isClosed) {
					end--;
				}
				
				for(int j=0;j<end;++j) {
					crvPt p = pts2->at(j);
					crvPt p2 = pts2->at((j+1) % pts2->size());
					float crvDist = abs(crvDistance(p,p2));
					if(crvDist < holeSpacing * .5) {
						if(merged) {
							p.t = firstMerge.t;
							crvPt mPt = mergePts(p, p2);
							pts2->at((j+1) % pts2->size()) = mPt;
						} else {
							crvPt mPt = mergePts(p,p2);
							firstMerge = p;
							pts2->at((j+1) % pts2->size()) = mPt;
							merged = true;
						}
					} else if(crvDist > holeSpacing*1.5) {//split
						pts->push_back(p);
						crvPt mPt = getPt(line,p.t+crvDist*0.5);
						intersect_parents(p.parents,p2.parents, mPt.parents);
						pts->push_back(mPt);
						merged = false;
					} else {
						pts->push_back(p);
						merged = false;
					}

				}



				if(!line->isClosed) {
					pts->push_back(pts2->back());
					crvPt mPt = getPt(line,line->length);
					pts->push_back(mPt);
				}

				//need to split points

				cout << "pts2 " << pts2->size() << endl;
				*crv = *pts;
			}
		}

		holePts.push_back(container);
	}
}
	
void getHolePts2() {
	
	for(int k=0;k<geodesics.size();++k) {
		vector<polyLine *> &crvs = geodesics[k];
		vector<crvPt> pts;
		for(int j=0;j<crvs.size();++j) {
			polyLine * pline = crvs[j];
			int numPts = (int) (pline->length/holeSpacing);
			float spacing = pline->length/numPts;
			int currCrv = 0;
			crvPt currPt = pline->pts[0];
			currPt.t = 0;
			currPt.crv = pline;
			pts.push_back(currPt);
			float currDist = 0;
			for(int i=0;i<numPts-1;++i) {
				while(currDist < spacing) {
					float lineDist = currPt.distance(pline->pts[currCrv+1]);
					if(currDist + lineDist > spacing) {
						ofVec3f dir = pline->pts[currCrv+1]-currPt;
						dir *= (spacing-currDist)/lineDist;
						currPt += dir;
						currPt.t = (i+1)*spacing;
						pts.push_back(currPt);
						currDist = 0;
						break;
					}
					currDist += lineDist;
					currPt = pline->pts[currCrv+1];
					currCrv += 2;
				}
			}
			if(!pline->isClosed) {
				currPt = pline->pts.back();
				currPt.crv = pline;
				currPt.t = pline->length;
				pts.push_back(currPt);
			}
		}
		holePts2.push_back(pts);
	}
}

void testApp::optimizePts() {
	/*
	crvVec3 pt;
	for(int i=0;i<geodesics.size()-1;++i) {
		vector<crvVec3> & crv = geodesics[i+1];
		vector<crvVec3> newPts;
		vector<crvVec3> & pts1 = holePts[i];
		vector<crvVec3> & pts2 = holePts[i+1];
		vector<crvVec3> projPts;
		float cLen = crvLen(crv);
		for(int j=0;j<pts1.size();++j) {
			projPts.push_back(projectPt(pts1[j],crv));
		}

		for(int j=0;j<pts2.size();++j) {
			pt = pts2[j];

			crvVec3 leftPt = pts2[(j-1+pts2.size())%pts2.size()];
			crvVec3 rightPt = pts2[(j+1)%pts2.size()];

			//get projected pts
			crvVec3 closestL, closestR;
			float minL = -99999;
			float minR = 999999;
			for(int k=0;k<projPts.size();++k) {
				float d = crvDistance(pt,projPts[k], cLen);
				if(d < 0) {
					if(d > minL) {
						minL = d;
						closestL = projPts[k];
					}
				} else {
					if(d < minR) {
						minR = d;
						closestR = projPts[k];
					}
				}
			}
			float neighDist = 0;
			neighDist += crvDistance(pt,leftPt, cLen);
			neighDist += crvDistance(pt,rightPt, cLen);
			neighDist += .5*minL;
			neighDist += .5*minR;

			neighDist /= 3.0;
			float newPos = pt.t+neighDist;

			newPts.push_back(getPt(crv,newPos));
		}

		holePts[i+1] = newPts;
	}
	*/
}

void testApp::exportDistObj() {
	ofstream out;
	out.open(ofToDataPath("dist.obj"));
	
	for(int i=0;i<hmesh.vertices.size();++i) {
		ofVec3f v = mesh.getVertex(i);
		out << "v " << v.x << " " << v.y << " " << v.z << " " << distances[i] << endl;
	}

	for(int i=0;i<mesh.getNumIndices();) {
		out << "f " << mesh.getIndex(i++)+1;
		out << " " << mesh.getIndex(i++)+1;
		out << " " << mesh.getIndex(i++)+1 << endl;
	}

	out.flush();
	out.close();
}

void testApp::loadSettings() {
	settings.load("settings.xml");
	dt = settings.getValue("dt",dt);
	diffusionRateA = settings.getValue("diffusionRateA",diffusionRateA);
	diffusionRateB = settings.getValue("diffusionRateB",diffusionRateB);
	F = settings.getValue("F",F);
	k = settings.getValue("k",k);
	offset = settings.getValue("offset",offset);
	threshold = settings.getValue("threshold",threshold);
	file = settings.getValue("file",file);
	geodesicSpacing = settings.getValue("geodesicSpacing",geodesicSpacing);
	holeSpacing = settings.getValue("holeSpacing",holeSpacing);
	doBorder = settings.getValue("doBorder",1)==1;
	cout << doBorder << endl;
}

void testApp::exit() {
	//settings.save("settings.xml");
}

void testApp::initValues() {
	for(int i=0;i<hmesh.vertices.size();++i) {
		a[i] = 1+ofRandom(-0.01,0);
		b[i] = 0+ofRandom(0,.01);
	}
	int numSamples = 0;//900;
	float radius = 2.0;
	float radiusSq = radius*radius;
	int nv = mesh.getNumVertices();
	for(int i=0;i<numSamples;++i) {
		//get pt
		ofVec3f pt = mesh.getVertex((int) ofRandom(nv));
		for(int j=0;j<nv;++j) {
			ofVec3f pt2 = mesh.getVertex(j);
			if(pt2.squareDistance(pt) < radiusSq) {
				a[j] = 0.5+ofRandom(-.01,0.01);
				b[j] = 0.25+ofRandom(-.01,0.01);
			}
		}
	}
}

void updateMeshNormals(ofVboMesh &mesh) {
	int numNormals = mesh.getNumNormals();
	ofVec3f * normals  = mesh.getNormalsPointer();
	for(int i=0;i<numNormals;++i) {
		normals[i] = ofVec3f(0,0,0);
	}
	ofVec3f p1,p2,p3;
	for(int i=0;i<mesh.getNumIndices();) {
		int i1 = mesh.getIndex(i++);
		int i2 = mesh.getIndex(i++);
		int i3 = mesh.getIndex(i++);

		p1 = mesh.getVertex(i1);
		p2 = mesh.getVertex(i2);
		p3 = mesh.getVertex(i3);

		p2 -= p1;
		p3 -= p1;
		p2.cross(p3);
		p2.normalize();
		normals[i1] -= p2;
		normals[i2] -= p2;
		normals[i3] -= p2;
	}
	for(int i=0;i<numNormals;++i) {
		normals[i].normalize();
	}
	mesh.getVbo().updateNormalData((float *)normals, numNormals);
}
//--------------------------------------------------------------
void testApp::update() {
	if(processing) {
		if(fileIndex >= files.size()) {
			processing = false;
		} else {
			mesh.load(files[fileIndex].getAbsolutePath());
			cout << files[fileIndex].getAbsolutePath() << endl;
			//process();
			mesh.save(files[fileIndex].getAbsolutePath());
			fileIndex++;
		}
	}
}

float reactionA(float aVal, float bVal, float reactionRate = 1.0f) {
  return reactionRate*(-aVal*bVal*bVal+F*(1-aVal));
}

float reactionB(float aVal, float bVal, float reactionRate = 1.0f) {
  return reactionRate*(aVal*bVal*bVal-(F+k)*bVal);
}


void testApp::reactionStep() {
	for(int i=0;i<a.size();++i) {
		float rate = (1.0-distances[i]/maximumDist)*28.0+0.5;
		float aVal = a[i];
		float bVal = b[i];
		
		float k1A = dt*reactionA(aVal,bVal,rate);
		float k1B = dt*reactionB(aVal,bVal,rate);
		float k2A = dt*reactionA(aVal+.5*k1A,bVal+.5*k1B,rate);
		float k2B = dt*reactionB(aVal+.5*k1A,bVal+.5*k1B,rate);
		float k3A = dt*reactionA(aVal+.5*k2A,bVal+.5*k2B,rate);
		float k3B = dt*reactionB(aVal+.5*k2A,bVal+.5*k2B,rate);
		float k4A = dt*reactionA(aVal+k3A,bVal+k3B);
		float k4B = dt*reactionB(aVal+k3A,bVal+k3B);		

		a[i] += 1.0/6.0*(k1A+2*k2A+2*k3A+k4A);
		b[i] += 1.0/6.0*(k1B+2*k2B+2*k3B+k4B);
		
		//a[i] += dt*reactionA(aVal,bVal);
		//b[i] += dt*reactionB(aVal,bVal);
		a[i] = min(1.0f,max(a[i],0.0f));
		b[i] = min(1.0f,max(b[i],0.0f));
	}
}

void testApp::updateSolver(float dt) {
	//Eigen::Map<Vectorf> aV(a.data(),hmesh.vertices.size());
	//Eigen::Map<Vectorf> bV(b.data(),hmesh.vertices.size());

	//aV = solverA.solveWithGuess(A*aV,aV);
	//bV = solverB.solveWithGuess(A*bV,bV);
		
	//reactionStep();

	cout << solverA.error() << " i " << solverA.iterations() << endl;
}

void testApp::adaptiveSubdivision() {
	float minLength = .2;
	float maxLength = 1.5;
	hmesh.facesToUpdate.clear();
	for(int i=0;i<15;++i) {
		int numEdges = hmesh.edges.size();
		for(int j=0;j<numEdges;++j) {
			heedge * e = hmesh.edges[j];
			hedge * he = e->edge;
			hedge * hep = he->pair;
			hevertex * v1 = he->vertex;
			hevertex * v2 = hep->vertex;

			ofVec3f p1 = mesh.getVertex(v1->index);
			ofVec3f p2 = mesh.getVertex(v2->index);

			float d1 = distances[v1->index];
			float d2 = distances[v2->index];
			float tLen = (d1+d2)*0.5/maximumDist*(maxLength-minLength)+minLength;
			float len = e->length;
			if(len > tLen) {
				//subdivide
				bool largestEdge = true;
				if(he->face != NULL) {
					if(he->next->edge->length > len) largestEdge = false;
					if(he->next->next->edge->length > len) largestEdge = false;
				}
				if(largestEdge && hep->face != NULL) {
					if(hep->next->edge->length > len) largestEdge = false;
					if(hep->next->next->edge->length > len) largestEdge = false;
				}
				if(largestEdge) {
					hmesh.subdivideEdge(e->edge);
				}
			}
		}
	}
	while(hmesh.mesh->getNumIndices() < hmesh.faces.size()*3) {
		hmesh.mesh->addIndex(0);
	}
	for(auto f : hmesh.facesToUpdate) {
		//indices
		int index = f->index*3;
		hedge * he = f->edge;
		hmesh.mesh->setIndex(index,he->vertex->index);
		he = he->next;
		hmesh.mesh->setIndex(index+1,he->vertex->index);
		he = he->next;
		hmesh.mesh->setIndex(index+2,he->vertex->index);
	}
}

void testApp::assembleMatrix() {
    // DEC
    SparseMatrix star0;
    HodgeStar0Form::build( hmesh, star0 );
    
    SparseMatrix  star1;
    HodgeStar1Form::build( hmesh, star1 );
         
    SparseMatrix d0;
    ExteriorDerivative0Form::build( hmesh, d0 );
    
	int nV = hmesh.vertices.size();
	SparseMatrix id(nV,nV);
	id.setIdentity();
    // zero Neumann boundary condition

	SparseMatrix L =  d0.transpose() * star1 * d0;
    
    //make L positive-definite
	//changed from 1.0e-8 to 1.0e-6 and it magically worked
    //L += 1.0e-6*star0;
	//A = star0;
    // heat flow for short interval
    //laplacianA = A + dt * diffusionRateA * L;
    //laplacianB = A + dt * diffusionRateB * L;
	
	//solverA.compute(laplacianA);
	//solverB.compute(laplacianB);
	//laplacian =  star0 * d0.transpose() * star1 * d0;
    
    
    // heat flow for short interval
	//laplacian = id-dt*laplacian;

}

float EPS = .0000001;

bool addSegment(polyLine * line, list<segment> & segs) {
	ofVec3f endPt = line->pts.back();
	for(list<segment>::iterator it = segs.begin();it != segs.end();it++) {
		if(endPt.squareDistance(it->first) < EPS) {
			crvPt newPt;
			newPt = it->second;
			newPt.crv = line;
			line->pts.push_back(newPt);
			//remove found seg from segs;
			segs.erase(it);
			return true;
		}
	}
	return false;
}

bool addSegmentReverse(polyLine * line, list<segment> & segs) {
	ofVec3f endPt = line->pts.back();
	for(list<segment>::iterator it = segs.begin();it != segs.end();it++) {
		if(endPt.squareDistance(it->second) < EPS) {
			crvPt newPt;
			newPt.crv = line;
			newPt = it->first;
			line->pts.push_back(newPt);
			//remove found seg from segs;
			segs.erase(it);
			return true;
		}
	}
	return false;
}

void testApp::getGeodesics() {
	list<segment> crv;
	float dist = geodesicSpacing;
	float halfDistSq = dist*dist/16;
	float bDistSq = 6*6;
	for(int k=0;k<80;++k) {
		vector<polyLine *> lines;
		crv.clear();
		dist = (k+.1)*geodesicSpacing;
		for(int i=0;i<mesh.getNumIndices();) {
			int i1 = mesh.getIndex(i++);
			int i2 = mesh.getIndex(i++);
			int i3 = mesh.getIndex(i++);

			float d1 = distances[i1];
			float d2 = distances[i2];
			float d3 = distances[i3];
		
			int below = 0, above = 0;
			if(d1 < dist) below++;
			else if(d1>dist) above++;
			if(d2 < dist) below++;
			else if(d2>dist) above++;
			if(d3 < dist) below++;
			else if(d3>dist) above++;
			if(below > 0 && above > 0) {
				ofVec3f p1 = mesh.getVertex(i1);
				ofVec3f p2 = mesh.getVertex(i2);
				ofVec3f p3 = mesh.getVertex(i3);

				ofVec3f n1 = mesh.getNormal(i1);
				ofVec3f n2 = mesh.getNormal(i2);
				ofVec3f n3 = mesh.getNormal(i3);
				
				crvVec3 int1,int2;
				crvVec3 inorm1,inorm2;
				if(below == 1) {
					if(d1 < dist) {
						int2 = abs(d2-d1) < .000001 ? p1 : (dist-d1)/(d2-d1)*(p2-p1)+p1;
						int1 = abs(d3-d1) < .000001 ? p1 : (dist-d1)/(d3-d1)*(p3-p1)+p1;
					} else if(d2 < dist) {
						int1 = abs(d1-d2) < .000001 ? p2 : (dist-d2)/(d1-d2)*(p1-p2)+p2;
						int2 = abs(d3-d2) < .000001 ? p2 : (dist-d2)/(d3-d2)*(p3-p2)+p2;
					} else {
						int1 = abs(d2-d3) < .000001 ? p3 : (dist-d3)/(d2-d3)*(p2-p3)+p3;
						int2 = abs(d1-d3) < .000001 ? p3 : (dist-d3)/(d1-d3)*(p1-p3)+p3;
					}
				} else {
					if(d1 > dist) {
						int1 = abs(d2-d1) < .000001 ? p1 : (dist-d1)/(d2-d1)*(p2-p1)+p1;
						int2 = abs(d3-d1) < .000001 ? p1 : (dist-d1)/(d3-d1)*(p3-p1)+p1;
					} else if(d2 > dist) {
						int2 = abs(d1-d2) < .000001 ? p2 : (dist-d2)/(d1-d2)*(p1-p2)+p2;
						int1 = abs(d3-d2) < .000001 ? p2 : (dist-d2)/(d3-d2)*(p3-p2)+p2;
					} else {
						int2 = abs(d2-d3) < .000001 ? p3 : (dist-d3)/(d2-d3)*(p2-p3)+p3;
						int1 = abs(d1-d3) < .000001 ? p3 : (dist-d3)/(d1-d3)*(p1-p3)+p3;
					}
				}
				if(int1.distanceSquared(int2) >= EPS) {
					crv.push_back(segment(int1,int2));
					//crv.push_back(int1);
					//crv.push_back(int2);
				}

			}
		}
		if(crv.size() > 0) {
			while(crv.size() > 0) {
				polyLine * newLine = new polyLine();
				newLine->pts.push_back(crv.front().first);
				while(addSegment(newLine,crv)) {

				}
				if(newLine->pts.front().squareDistance(newLine->pts.back()) < EPS) {
					newLine->pts.pop_back();
					newLine->isClosed = true;
				} else {
					reverse(newLine->pts.begin(), newLine->pts.end());
					while(addSegmentReverse(newLine,crv)) {

					}
					reverse(newLine->pts.begin(), newLine->pts.end());
					newLine->isClosed = false;
				}
				newLine->getLength();
				lines.push_back(newLine);
			}
			geodesics.push_back(lines);
		}
	}
}


//--------------------------------------------------------------
void testApp::draw(){
	//exitApp();
	//ofDisableLighting();

	mCam.begin();
	//ofLight light0;
	//light0.enable();
	//light0.setPosition(-15,9.3,0);
	//light0.setDiffuseColor(ofFloatColor(.8,.8,.8));
	//ofLight light1;
	//light1.enable();
	//light1.setPosition(26.879,-3.8,0);
	//ofSetColor(110);
	//ofSetColor(255);
	//glEnable(GL_COLOR_MATERIAL);
	//glColorMaterial(GL_BACK, GL_AMBIENT_AND_DIFFUSE);
	
	//glColor4f(44.0/255.0,110.0/255.0,0.0,1.0);
	//glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	
	//glColor4f(1.0,251.0/255.0,0.0,1.0);
	glEnable(GL_DEPTH_TEST);
	
	mesh.enableColors();
	//float frontColor[4] = {44.0/255.0,110.0/255.0,0.0,1.0};
	//glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, frontColor);
	//float backColor[4] = {1.0,251.0/255.0,0.0,1.0};
	//glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, backColor);
	for(int i=0;i<hmesh.vertices.size();++i) {
		//mesh.setColor(i,ofColor(a[i]*255));
		//mesh.setColor(i,ofColor::fromHsb(distances[i]/maximumDist*255,200.0,a[i]*255));
		mesh.setColor(i,ofColor::fromHsb(distances[i]/maximumDist*255,200.0,255));
	}
	mesh.draw();
}

void saveLines() {
	ostringstream filename;
	filename<< file.substr(0,file.length()-4)<<geodesicSpacing<<"_"<<holeSpacing<<".csv";

	ofstream out(ofToDataPath(filename.str()));
	//ofstream out(ofToDataPath(file + "_holes.csv"));
	for(int i=1;i<holePts.size();++i) {
		for(auto pts : *holePts[i]) {
			//vector<crvPt> & pts = holePts[i];
			for(int j=0;j<pts->size();++j) {
				crvPt p1 = pts->at(j);
				for(int k=0;k<p1.parents.size();++k) {
					crvPt p2 = *p1.parents[k];
					out << p1.x << "," << p1.y << "," << p1.z << "," << p2.x << "," << p2.y << "," << p2.z << "," << thickness*0.5 << endl; 
				}
			}
		}
	}
	if(doBorder) {
		for(int j=0;j<geodesics[0].size();++j) {
			polyLine * crv = geodesics[0][j];
			for(int i=0;i<crv->pts.size();++i) {
				crvPt p1 = crv->pts[i];
				crvPt p2 = crv->pts[(i+1)%crv->pts.size()];
				out << p1.x << "," << p1.y << "," << p1.z << "," << p2.x << "," << p2.y << "," << p2.z << "," << thickness*0.5 << endl;
			}
		}
	}
	/*
	for(int i=0;i<geodesics.size()-1;++i) {
		
		vector<polyLine *> crvs = geodesics[i];
		vector<crvPt> &crv = holePts[i];
		vector<crvPt> &crv2 = holePts[i+1];

		vector<crvPt> projPts;
		for(int j=0;j<crv2.size();++j) {
			projPts.push_back(projectPt(crv2[j],crvs));
		}

		for(int j=0;j<crv.size();++j) {
			crvPt p1 = crv[j];
			for(int k=0;k<crv2.size();++k) {
				crvPt p2 = projPts[k];
				if(p1.crv == p2.crv) {
					if(abs(crvDistance(projPts[k], p1)) < holeSpacing) {
						p2 = crv2[k];
						out << p1.x << "," << p1.y << "," << p1.z << "," << p2.x << "," << p2.y << "," << p2.z << "," << thickness*0.5 << endl; 
					}
				}
				*/
					/*
				if(d < 0) {
					if(d > minL) {
						minL = d;
						closestL = crv[k];
					}
				} else {
					if(d < minR) {
						minR = d;
						closestR = crv[k];
					}
				}
				
			}
			//ofLine(p1,closestL);
			//ofLine(p1,closestR);
			//if(p1.squareDistance(closestL) < 100) {
			//	out << p1.x << "," << p1.y << "," << p1.z << "," << closestL.x << "," << closestL.y << "," << closestL.z << "," << thickness*0.5 << endl;
			//}
			//if(p1.squareDistance(closestR) < 100) {
			//	out << p1.x << "," << p1.y << "," << p1.z << "," << closestR.x << "," << closestR.y << "," << closestR.z << "," << thickness*0.5 << endl; 
			//}
			/*
			for(int k=0;k<crv2.size();++k) {
				ofVec3f p2 = crv2[k];
				//if(p1.distanceSquared(projPts[k]) < 4) {
				if(abs(crvDistance(projPts[k], p1,crv1Len)) < 2.0f) {
					ofLine(p1,p2);
				}
			}
		
			
		}
	}
	*/
	out.flush();
	out.close();
}

void testApp::exportHoles() {
	float rThreshold = 1.0-threshold;
	ofMesh holeMesh;
	//should do the dual but for now we won't
	for(int i=0;i<mesh.getNumVertices();++i) {
		float d = distances[i];
		if(d < 2) {
			a[i] -= 2-d;
			a[i] = max(a[i],0.0f);
		}
		float aVal = a[i];
		bool in = aVal < threshold;

		ofVec3f norm = mesh.getNormal(i);
		ofVec3f pos = mesh.getVertex(i);

		aVal = max(aVal,.3f);
		float normVal = (threshold-aVal)/(threshold-.3);

		holeMesh.addVertex(pos+norm*offset*normVal);
		holeMesh.addVertex(pos-norm*offset*normVal);
	}
	for(int i=0;i<mesh.getNumIndices();) {
		int i1 = mesh.getIndex(i++);
		int i2 = mesh.getIndex(i++);
		int i3 = mesh.getIndex(i++);

		bool in1 = a[i1] < threshold;
		bool in2 = a[i2] < threshold;
		bool in3 = a[i3] < threshold;
		bool b1 = hmesh.vertices[i1]->boundary;
		bool b2 = hmesh.vertices[i2]->boundary;
		bool b3 = hmesh.vertices[i3]->boundary;


		if(in1 && in2 && in3) {
			holeMesh.addIndex(i1*2);
			holeMesh.addIndex(i2*2);
			holeMesh.addIndex(i3*2);
			holeMesh.addIndex(i1*2+1);
			holeMesh.addIndex(i3*2+1);
			holeMesh.addIndex(i2*2+1);
		} else if(in1 && in2) {
			holeMesh.addIndex(i1*2);
			holeMesh.addIndex(i2*2+1);
			holeMesh.addIndex(i1*2+1);
			holeMesh.addIndex(i1*2);
			holeMesh.addIndex(i2*2);
			holeMesh.addIndex(i2*2+1);
		} else if(in2 && in3) {
			holeMesh.addIndex(i2*2);
			holeMesh.addIndex(i3*2+1);
			holeMesh.addIndex(i2*2+1);
			holeMesh.addIndex(i2*2);
			holeMesh.addIndex(i3*2);
			holeMesh.addIndex(i3*2+1);
		} else if(in1 && in3) {
			holeMesh.addIndex(i3*2);
			holeMesh.addIndex(i1*2+1);
			holeMesh.addIndex(i3*2+1);
			holeMesh.addIndex(i3*2);
			holeMesh.addIndex(i1*2);
			holeMesh.addIndex(i1*2+1);
		}
		if(b1 && b2) {
			holeMesh.addIndex(i1*2);
			holeMesh.addIndex(i1*2+1);
			holeMesh.addIndex(i2*2+1);
			holeMesh.addIndex(i1*2);
			holeMesh.addIndex(i2*2+1);
			holeMesh.addIndex(i2*2);
		} else if(b2 && b3) {
			holeMesh.addIndex(i2*2);
			holeMesh.addIndex(i2*2+1);
			holeMesh.addIndex(i3*2+1);
			holeMesh.addIndex(i2*2);
			holeMesh.addIndex(i3*2+1);
			holeMesh.addIndex(i3*2);
		} else if(b1 && b3) {
			holeMesh.addIndex(i3*2);
			holeMesh.addIndex(i3*2+1);
			holeMesh.addIndex(i1*2+1);
			holeMesh.addIndex(i3*2);
			holeMesh.addIndex(i1*2+1);
			holeMesh.addIndex(i1*2);
		}
	}
	holeMesh.save("holes.ply");
}


void testApp::exportGeo() {
	ofMesh holeMesh;
	//sort geodesics
	vector<ofVec3f> newCrv;
	ofVec3f startPt, pt;
	/*
	for(int i=0;i<geodesics.size();++i) {
		vector<ofVec3f> & crv = geodesics[i];
		newCrv.clear();
		newCrv.push_back(crv[0]);
		newCrv.push_back(crv[1]);
		startPt.set(newCrv[0]);
		pt.set(newCrv[1]);
		while(pt.distanceSquared( startPt) >= EPS) {
			//get next pt
			bool found = false;
			for(int j=2;j<crv.size();j+=2) {
				if(crv[j].distanceSquared(pt) < EPS) {
					pt = crv[j+1];
					found = true;
					break;
				}
			}
			if(!found) {
				break;
			}
			newCrv.push_back(pt);
		}
		cout << crv.size() << " " << newCrv.size() << endl;

		crv = newCrv;

	}
	*/

	holeMesh.save("ribs.ply");
}

/*
void testApp::exportGeo() {
	ofMesh holeMesh;
	float maxLenSq = 1.2*1.2;
	//subdivide
	
	for(int i=0;i<mesh.getNumIndices();) {
		int i1 = mesh.getIndex(i++);
		int i2 = mesh.getIndex(i++);
		int i3 = mesh.getIndex(i++);

		float d1 = distances[i1];
		float d2 = distances[i2];
		float d3 = distances[i3];

		//check distances
		float val1 = cos((d1+5)/10.0*TWO_PI);
		float val2 = cos((d2+5)/10.0*TWO_PI);
		float val3 = cos((d3+5)/10.0*TWO_PI);

		int inCount = (val1 > .5) + (val2 > .5) +(val3 > .5);
		if(inCount > 0) {
			//check lengths
			ofVec3f p1 = mesh.getVertex(i1);
			ofVec3f p2 = mesh.getVertex(i2);
			ofVec3f p3 = mesh.getVertex(i3);

			if(p1.distanceSquared(p2) > maxLenSq || p1.distanceSquared(p3) > maxLenSq || p2.distanceSquared(p3) > maxLenSq) {
				//subdivide
				ofVec3f mid1 = (p2+p3)*0.5;
				ofVec3f mid2 = (p1+p3)*0.5;
				ofVec3f mid3 = (p1+p2)*0.5;

				ofVec3f norm1 = mesh.getNormal(i1);
				ofVec3f norm2 = mesh.getNormal(i2);
				ofVec3f norm3 = mesh.getNormal(i3);

				ofVec3f nNorm1 = (norm2+norm3).getNormalized();
				ofVec3f nNorm2 = (norm1+norm3).getNormalized();
				ofVec3f nNorm3 = (norm2+norm1).getNormalized();

				float nd1 = (d2+d3)*0.5;
				float nd2 = (d1+d3)*0.5;
				float nd3 = (d2+d1)*0.5;

				int ni1 = mesh.getNumVertices();
				mesh.addVertex(mid1);
				mesh.addNormal(nNorm1);
				int ni2 = mesh.getNumVertices();
				mesh.addVertex(mid2);
				mesh.addNormal(nNorm2);
				int ni3 = mesh.getNumVertices();
				mesh.addVertex(mid3);
				mesh.addNormal(nNorm3);

				distances.push_back(nd1);
				distances.push_back(nd2);
				distances.push_back(nd3);

				//indices

				mesh.addIndex(ni1);mesh.addIndex(ni2);mesh.addIndex(ni3);
				mesh.addIndex(i2);mesh.addIndex(ni1);mesh.addIndex(ni3);
				mesh.addIndex(i3);mesh.addIndex(ni2);mesh.addIndex(ni1);
				mesh.setIndex(i-3, i1); mesh.setIndex(i-2,ni3);	mesh.setIndex(i-1,ni2);


				i-=3;
			}
		}
	}
	
	for(int i=0;i<mesh.getNumVertices();++i) {
		float d = distances[i]+5;
		d = max(d,0.0f);

		ofVec3f norm = mesh.getNormal(i);
		ofVec3f pos = mesh.getVertex(i);

		float normVal = cos(d/10.0*TWO_PI);
		if(normVal>.5) {
			normVal = (cos(d/10*3*TWO_PI)+1)*0.5;
		} else {
			normVal = 0;
		}
		pos += norm*(normVal*1.5+.6);

		holeMesh.addVertex(pos);
		pos -= norm*(normVal*1.5+1.6);
		holeMesh.addVertex(pos);
	}
	for(int i=0;i<mesh.getNumIndices();) {
		int i1 = mesh.getIndex(i++);
		int i2 = mesh.getIndex(i++);
		int i3 = mesh.getIndex(i++);

		holeMesh.addIndex(i1*2);
		holeMesh.addIndex(i2*2);
		holeMesh.addIndex(i3*2);
		holeMesh.addIndex(i1*2+1);
		holeMesh.addIndex(i3*2+1);
		holeMesh.addIndex(i2*2+1);

	}

	holeMesh.save("ribs.ply");
}
*/

int numColors = 8;
int numColors2 = 5;
//float colorD[8] = {0,8,12,16,25,35,60,95};
float colorD[8] = {0,8*2,12*2,16*2,25*2,35*2,60*2,95*2};
float colorD2[6] = {0,12,35,60,90};

ofVec3f colors[8] = {
ofVec3f(254,254,1),
ofVec3f(254,154,1),
ofVec3f(254,91,1),
ofVec3f(189,16,1),
ofVec3f(208,65,141),
ofVec3f(147,52,149),
ofVec3f(92,78,152),
ofVec3f(91,135,186)
};

ofVec3f colors2[5] = {
ofVec3f(62,114,183),
ofVec3f(27,188,239),
ofVec3f(7,46,109),
ofVec3f(119,204,204),
ofVec3f(16,48,34)
};



float cubicInterp(float d, float x0, float x1, float x2, float x3,float d0, float d1, float d2,float d3) {
	float t = (d-d1)/(d2-d1);
	float m0 = (x2-x0)/sqrt(d2-d0);
	float m1 = (x3-x1)/sqrt(d3-d1);
   //float a0,a1,a2,a3;
   //a0 = -0.5*x0 + 2*x1 - 2.5*x2 + x3;
   //a1 = x0 - 2*x1 + 3*x2 - x3;
   //a2 = -x0 +x1;
   //a3 = x1;
   //return (a0*t*t*t+a1*t*t+a2*t+a3);
	float t2 = t*t;
	float t3 = t2*t;
	return (2*t3-3*t2+1)*x1+
			(t3-2*t2+t)*m0+
			(-2*t3+3*t2)*x2+
			(t3-t2)*m1;
}

ofColor cubicInterpColor(float d, ofVec3f c0, ofVec3f c1, ofVec3f c2, ofVec3f c3,float d0, float d1, float d2,float d3) {
	float t = (d-d1)/(d2-d1);
	ofVec3f m0 = (c2-c0)/sqrt(d2-d0);
	ofVec3f m1 = (c3-c1)/sqrt(d3-d1);
	float t2 = t*t;
	float t3 = t2*t;
	ofVec3f c = (2*t3-3*t2+1)*c1+
			(t3-2*t2+t)*m0+
			(-2*t3+3*t2)*c2+
			(t3-t2)*m1;
	return ofColor(ofClamp(c.x,0,255),ofClamp(c.y,0,255),ofClamp(c.z,0,255)).getClamped();
}

ofColor getColor1(float d) {
	for(int i=0;i<numColors;++i) { //7
		if(d<colorD[i]) {
			int i0 = max(0,i-2);
			int i1 = max(0,i-1);
			int i2 = min(i,numColors-1);
			int i3 = min(i+1,numColors-1);
			return cubicInterpColor(d,colors[i0],colors[i1],colors[i2],colors[i3],colorD[i0],colorD[i1],colorD[i2],colorD[i3]);
		}
	}
	int i0 = max(0,numColors-3);
	int i1 = max(0,numColors-2);
	int i2 = min(numColors-1,numColors-1);
	int i3 = min(numColors,numColors-1);
	return cubicInterpColor(d,colors[i0],colors[i1],colors[i2],colors[i3],colorD[i0],colorD[i1],colorD[i2],colorD[i3]);
	ofColor c(100,217,255);
	c.lerp(ofColor(0,175,230),min(d/20.0f,1.0f));
	return c;
}

ofColor getColor2(float d) {
	return getColor1(d);
	for(int i=0;i<numColors2;++i) { //7
		if(d<colorD2[i]) {
			int i0 = max(0,i-2);
			int i1 = max(0,i-1);
			int i2 = min(i,numColors2-1);
			int i3 = min(i+1,numColors2-1);
			return cubicInterpColor(d,colors2[i0],colors2[i1],colors2[i2],colors2[i3],colorD2[i0],colorD2[i1],colorD2[i2],colorD2[i3]);
		}
	}
	int i0 = max(0,numColors2-3);
	int i1 = max(0,numColors2-2);
	int i2 = numColors2-1;
	int i3 = numColors2-1;
	return cubicInterpColor(d,colors2[i0],colors2[i1],colors2[i2],colors2[i3],colorD2[i0],colorD2[i1],colorD2[i2],colorD2[i3]);
}

float getThickness(float d) {
	//return d/maximumDist*.75+.75;
	//return .75;
	//if(d < edgeDist) {
		//return cubicInterp(d/2.0,.5, .5, 1, 1.5);
		//return cubicInterp(d,.55,.55,1,1.5,0,0,3,25); 
		//return cubicInterp(d,edgeThick*0.25,edgeThick*0.5,edgeThick*0.5,midBodyThick*0.5,-edgeDist,0,edgeDist,midBodyDist);
	//} else 
	if(d < midBodyDist) {
		//float t = (d-2)/(25.0f-2.0f);
		//return cubicInterp(t,.5, 1, 1.5, 2);
		//return cubicInterp(d,.5,1,1.5,2,0,2,25,maximumDist);
		//return cubicInterp(d,edgeThick*0.5,edgeThick*0.5,midBodyThick*0.5,bodyThick*0.5,0,edgeDist,midBodyDist,maximumDist);
		return cubicInterp(d,edgeThick*0.25,edgeThick*0.5,midBodyThick*0.5,bodyThick*0.5,-edgeDist,0,midBodyDist,maximumDist);
	} else {
		//float t = (d-25.0)/(maximumDist-25.0);
		//return cubicInterp(t,1, 1.5, 2, 2);
		d = min(d,maximumDist);
		//return cubicInterp(d,.75,1.5,2.,2., 2, 25,maximumDist,maximumDist);
		return cubicInterp(d, edgeThick*0.5, midBodyThick*0.5, bodyThick*0.5, bodyThick*0.5, edgeDist, midBodyDist, maximumDist, maximumDist);
	}
}

void testApp::exportColor() {
	ofMesh holeMesh;
	holeMesh.enableColors();
	ofColor bColor(255,255,255);
	int bIndex = 0;
	//should do the dual but for now we won't
	for(int i=0;i<mesh.getNumVertices();++i) {
		float d = distances[i];
		d = max(d,0.0f);

		float d2 = distances2[i];
		ofVec3f norm = mesh.getNormal(i);
		ofVec3f pos = mesh.getVertex(i);

		//collar -17 - 5
		float t_d = ofClamp((pos.y-minY)/(maxY-minY),0,1);
		rimThick = ofLerp(rimT[0], rimT[1], t_d);
		edgeThick = ofLerp(edgeT[0], edgeT[1], t_d);
		midBodyThick = ofLerp(midBodyT[0], midBodyT[1], t_d);
		bodyThick = ofLerp(bodyT[0], bodyT[1], t_d);
		float normVal = getThickness(d);

		if(d2 < 1.1) {
			normVal = ofLerp(normVal,rimThick*.5, pow(1-d2/1.1,2));
		}
		hmesh.vertices[i]->index = holeMesh.getNumVertices();
		//normVal = .75*normVal+.75;
		//if(d <= .0001) {
		//	normVal = .5;
		//}
		//holeMesh.addVertex(pos+norm*offset*normVal);
		//holeMesh.addVertex(pos-norm*offset*normVal);
		holeMesh.addVertex(pos+norm*normVal);
		holeMesh.addColor(getColor1(d));
		holeMesh.addVertex(pos-norm*normVal);	
		holeMesh.addColor(getColor2(d));
		
		bool b1 = hmesh.vertices[i]->boundary;
		if(b1) {
			holeMesh.addVertex(pos+norm*normVal);
			holeMesh.addColor(bColor);
			holeMesh.addVertex(pos-norm*normVal);	
			holeMesh.addColor(bColor);
			bIndex++;
		}
	}
	for(int i=0;i<mesh.getNumIndices();) {
		int i1 = mesh.getIndex(i++);
		int i2 = mesh.getIndex(i++);
		int i3 = mesh.getIndex(i++);

		hevertex * v1 = hmesh.vertices[i1];
		hevertex * v2 = hmesh.vertices[i2];
		hevertex * v3 = hmesh.vertices[i3];
		i1 = v1->index;
		i2 = v3->index;
		i3 = v2->index;
		bool b1 = v1->boundary;
		bool b2 = v3->boundary;
		bool b3 = v2->boundary;

		int bi1 = i1+2;
		int bi2 = i2+2;
		int bi3 = i3+2;

		holeMesh.addIndex(i1);
		holeMesh.addIndex(i2);
		holeMesh.addIndex(i3);
		holeMesh.addIndex(i1+1);
		holeMesh.addIndex(i3+1);
		holeMesh.addIndex(i2+1);
		
		if(b1 && b2) {
			holeMesh.addIndex(bi1);
			holeMesh.addIndex(bi1+1);
			holeMesh.addIndex(bi2+1);
			holeMesh.addIndex(bi1);
			holeMesh.addIndex(bi2+1);
			holeMesh.addIndex(bi2);
		} else if(b2 && b3) {
			holeMesh.addIndex(bi2);
			holeMesh.addIndex(bi2+1);
			holeMesh.addIndex(bi3+1);
			holeMesh.addIndex(bi2);
			holeMesh.addIndex(bi3+1);
			holeMesh.addIndex(bi3);
		} else if(b1 && b3) {
			holeMesh.addIndex(bi3);
			holeMesh.addIndex(bi3+1);
			holeMesh.addIndex(bi1+1);
			holeMesh.addIndex(bi3);
			holeMesh.addIndex(bi1+1);
			holeMesh.addIndex(bi1);
		}
	}
	ostringstream filename;
	filename<< file.substr(0,file.length()-4)<<"thick.ply";
	holeMesh.save(filename.str());
}

void testApp::loadMesh(string path) {
	mesh.load(path);
	//mesh.indices.reserve(mesh.getNumIndices()*20);
	//perturb
	for(int i=0;i<mesh.getNumVertices();++i) {
		ofVec3f pt = mesh.getVertex(i);
		pt += ofVec3f(ofRandom(-1,1),ofRandom(-1,1),ofRandom(-1,1))*0.0001;
		mesh.setVertex(i,pt);
	}
	//edges
	//getMaxEdgeLength();
}

void testApp::dragged(ofDragInfo & info) {
	if(!processing) {
		processing = true;
		ofDirectory dir;
		dir = ofDirectory(info.files[0]);
		dir.allowExt("ply");
		dir.listDir();
		files = dir.getFiles();
		fileIndex = 0;
	}
}

//--------------------------------------------------------------
void testApp::keyPressed(int key){
	switch(key) {
		case 'p':
			paused = !paused;
			break;
		case 'e':
			exportColor();
			break;
		case 's':
			//exportHoles();
			//exportGeo();
			saveLines();
			break;
		case 'o':
			optimizePts();
			optimizePts();
			break;
	}
}

//--------------------------------------------------------------
void testApp::keyReleased(int key){

}

//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y){

}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){ 

}