#pragma once

#include "ofMain.h"
//#include "ofxXmlSettings.h"
#include "solver.h"
#include "Hemesh.h"
#include "DiscreteExteriorCalculus.h"
#include "polyLine.h"

class testApp : public ofBaseApp{
	public:
		void setup();
		void update();
		void draw();
		void loadMesh(string path);

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y);
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
		void exit();

		//ofxXmlSettings settings;
		void loadSettings();
		ofEasyCam mCam;
		hemesh hmesh;
		ofVboMesh mesh;
		void reactionStep();
		void updateSolver(float dt);
		void initValues();
		void assembleMatrix();
		void exportHoles();
		void adaptiveSubdivision();
		void exportDistObj();
		void exportColor();
		void exportGeo();
		void getGeodesics();
		void getHolePts();
		void optimizePts();
};

void updateMeshNormals(ofVboMesh &mesh);
