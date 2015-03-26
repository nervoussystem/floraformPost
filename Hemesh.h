#pragma once

#include "ofMesh.h"
#include "ofVec3f.h"
#include "ofVec3d.h"

struct hedge;
struct heface;
struct hevertex;
struct heedge;

struct hedge {
	hedge * next;
	hedge * pair;
	hevertex * vertex;
	heedge * edge;
	heface * face;
	float cosA;
	ofVec3f norm,dir;
	float invAlt;
	unsigned int index;
	int tag;
	hedge() : next(NULL), pair(NULL), vertex(NULL), face(NULL),tag(0) {}

	float cotan();
};

struct heface {
	hedge * edge;
	double scalar;
	ofVec3f norm;
	ofVec3d vector;
	unsigned int index;
	//ofVec3f edgeNorm[3];
	float area;
	//float cosA[3];
	//float invAltitude[3];

	heface() : edge(NULL) {}
	
};

struct hevertex {
	hedge * edge;
	ofVec3f position;
	ofVec3f normal;
	unsigned int index;
	int flag;
	bool boundary;
	bool tag;
	hevertex(ofVec3f & pos) : edge(NULL), position(pos), boundary(false), tag(false) {}

	float area() {
		hedge * e = edge;
		float area = 0;
		do {
			if(e->face != NULL) {
				area += e->face->area;
			}
			e = e->next->pair;
		} while(e != edge);
		area /= 3.0;
		return area;
	}
};

struct heedge {
	hedge * edge;
	void * f;
	unsigned int index;
	float length;
	heedge() : edge(NULL) {};
};


class hemesh {
public:
	ofMesh * mesh;
	vector<heedge *> edges;
	vector<heface *> faces;
	vector<hevertex *> vertices;

	vector<heedge *> edgesToUpdate;
	vector<heface *> facesToUpdate;

	hevertex * subdivideEdge(hedge * e);

	hevertex * addVertex(ofVec3f & pos) {
		hevertex * v = new hevertex(pos);
		v->index = vertices.size();
		vertices.push_back(v);
		return v;
	}

	heface * addFace() {
		heface * f = new heface();
		f->index = faces.size();
		faces.push_back(f);
		return f;
	}

	heedge * addEdge() {
		heedge * e = new heedge();
		e->index = edges.size();
		edges.push_back(e);
		return e;
	}
	void addHole(heface * f);
};

hedge * getEdgePair(int index1, int index2, vector<vector<hedge * > > & edgeMap);
hemesh hemeshFromOfMesh(ofMesh & mesh);