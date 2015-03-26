#include "solver.h"
#include "Hemesh.h"

hevertex * hemesh::subdivideEdge(hedge * e) {
	hedge * epair = e->pair;
	hevertex *v1 = e->vertex;
	hevertex *v2 = epair->vertex;
	hedge *enext = e->next;
	hedge *epnext = epair->next;

	//hevertex * newVertex = addVertex(0.5*(v1->position+v2->position));
	ofVec3f newPos = 0.5f*(mesh->getVertex(v1->index)+mesh->getVertex(v2->index));
	hevertex * newVertex = addVertex(newPos);
	distances.push_back(0.5*(distances[v1->index]+distances[v2->index]));
	//distances.push_back((distances[v1->index]+distances[v2->index])*0.5);
	mesh->addVertex(newPos);
	mesh->addNormal(newPos);
	mesh->addColor(ofColor());
	//do e first
	hedge * newEdge = new hedge();
	newEdge->vertex = v1;
	v1->edge = newEdge;
	e->vertex = newVertex;
	newEdge->next = enext;
	newEdge->pair = epair;
	epair->pair = newEdge;
	//tag edges so we can track boundaries
	newEdge->tag = e->tag;
    
	hedge * newEdgePair = new hedge();
	newEdgePair->vertex = v2;
	epair->vertex = newVertex;
	newEdgePair->next = epnext;
	newEdgePair->pair = e;
	e->pair = newEdgePair;
	newVertex->edge = e;
	//tag edges so we can track boundaries
	newEdgePair->tag = epair->tag;
	//NEED TO DEAL WITH "EDGES"

	heedge * oldEdge = e->edge;
	//float rLength = ((Force *) (oldEdge->f))->restLength;
	float rLength = oldEdge->length;
	rLength *= 0.5;
	//((Force *) (oldEdge->f))->restLength = rLength;
	e->pair->edge = oldEdge;
	oldEdge->edge = e;
	oldEdge->length = rLength;

	heedge * redge = addEdge();
	Force * newForce = new Force();
	redge->f = newForce;
	newForce->restLength = rLength;
	redge->length = rLength;

	redge->edge = newEdge;
	newEdge->edge = redge;
	newEdge->pair->edge = redge;
    
	//set b to neighboring b, it p1.b should equal p2.b
	//if(e->face == NULL || epair->face == NULL) { newVertex->b = v1->b;}
	edgesToUpdate.push_back(redge);
	edgesToUpdate.push_back(e->edge);

	if(e->face != NULL) {
		//face 1
		heface * newFace = addFace();
		hedge * splitEdge1 = new hedge();
		hedge * splitEdge2 = new hedge();
		splitEdge1->pair = splitEdge2;
		splitEdge2->pair = splitEdge1;
		splitEdge1->vertex = enext->vertex;
		splitEdge2->vertex = newVertex;
      
		//e.f
		e->next = splitEdge1;
		hedge * enext2 = enext->next;
		splitEdge1->next = enext2;
		e->face->edge = e;
		splitEdge1->face = e->face;
		//newFace
		newEdge->face = newFace;
		splitEdge2->face = newFace;
		enext->face = newFace;
		newFace->edge = newEdge;
		enext->next = splitEdge2;
		splitEdge2->next = newEdge;

		heedge * sedge = addEdge();
		sedge->edge = splitEdge1;
		splitEdge1->edge = sedge;
		splitEdge2->edge = sedge;

		Force * newForce = new Force();
		sedge->f = newForce;
		//float rLength1 = ((Force *) enext->edge->f)->restLength;
		//float rLength2 = ((Force *) enext2->edge->f)->restLength;
		float rLength1 = enext->edge->length;
		float rLength2 = enext2->edge->length;
		newForce->restLength = sqrt(0.5*(rLength1*rLength1+rLength2*rLength2)-rLength*rLength);
		sedge->length = newForce->restLength;

		facesToUpdate.push_back(newFace);
		facesToUpdate.push_back(e->face);
		edgesToUpdate.push_back(sedge);
		edgesToUpdate.push_back(enext->edge);
		edgesToUpdate.push_back(enext2->edge);
	} else {
		e->next = newEdge;
	}
    
	if(epair->face != NULL) {
		heface * newFace = addFace();
		hedge * splitEdge1 = new hedge();
		hedge * splitEdge2 = new hedge();
		splitEdge1->pair = splitEdge2;
		splitEdge2->pair = splitEdge1;
		splitEdge1->vertex = epnext->vertex;
		splitEdge2->vertex = newVertex;
      
		//epair.f
		hedge * epnext2 = epnext->next;
		epair->next = splitEdge1;
		splitEdge1->next = epnext2;
		epair->face->edge = epair;
		splitEdge1->face = epair->face;
      
		//newFace
		newEdgePair->face = newFace;
		splitEdge2->face = newFace;
		epnext->face = newFace;
		newFace->edge = newEdgePair;
		epnext->next = splitEdge2;
		splitEdge2->next = newEdgePair;

		heedge * sedge = addEdge();
		sedge->edge = splitEdge1;
  		splitEdge1->edge = sedge;
		splitEdge2->edge = sedge;

		Force * newForce = new Force();
		sedge->f = newForce;
		//float rLength1 = ((Force *) epnext->edge->f)->restLength;
		//float rLength2 = ((Force *) epnext2->edge->f)->restLength;
		float rLength1 = epnext->edge->length;
		float rLength2 = epnext2->edge->length;
		newForce->restLength = sqrt(0.5*(rLength1*rLength1+rLength2*rLength2)-rLength*rLength);
		sedge->length = newForce->restLength;

		//stuff to update
		facesToUpdate.push_back(newFace);
		facesToUpdate.push_back(epair->face);
		edgesToUpdate.push_back(sedge);
		edgesToUpdate.push_back(epnext->edge);
		edgesToUpdate.push_back(epnext2->edge);
	} else {
		epair->next = newEdgePair;
	}
    
	return newVertex;
}

void hemesh::addHole(heface * f) {
	hedge * e = f->edge;
	hedge * pair;
	do {
		e->face = NULL;
		e = e->next;
	} while(e != f->edge);
	faces.erase(find(faces.begin(), faces.end(),f));
	//redo indices
	for(int i=0;i<faces.size();++i) {
		faces[i]->index = i;
	}
}

float hedge::cotan() {
	if(face == NULL) return 0.0;
	//zero area faces are bad
	return cosA/face->area*next->next->edge->length*next->edge->length*0.5;

}

hedge * getEdgePair(int index1, int index2, vector<vector<hedge * > > & edgeMap) {
	vector<hedge *> & map = edgeMap[index1];
	for(int i=0;i<map.size();++i) {
		hedge * e = map[i];
		if(e->vertex->index == index2) {
			return e;
		}
	}
	return NULL;
}

hemesh hemeshFromOfMesh(ofMesh & mesh) {
	hemesh newMesh;
	newMesh.mesh = &mesh;
	int numVertices = mesh.getNumVertices();
	vector<vector<hedge * > > edgeMap;
	edgeMap.resize(numVertices);
	for(int i=0;i<numVertices;++i) {
		newMesh.addVertex(mesh.getVertex(i));
	}
	for(int i=0;i<mesh.getNumIndices();) {
		int i0 = mesh.getIndex(i++);
		int i1 = mesh.getIndex(i++);
		int i2 = mesh.getIndex(i++);

		hevertex * v0 = newMesh.vertices[i0];
		hevertex * v1 = newMesh.vertices[i1];
		hevertex * v2 = newMesh.vertices[i2];

		hedge *e0, *e1, *e2;


		hedge * pair = getEdgePair(i0,i1,edgeMap);
		if(pair != NULL) {
			e0 = pair;
		} else {
			e0 = new hedge();
			pair  = new hedge();
			e0->vertex = v1;
			pair->vertex = v0;
			e0->pair = pair;
			pair->pair = e0;
			v1->edge = e0;

			edgeMap[i0].push_back(e0);
			edgeMap[i1].push_back(pair);
			heedge * e = newMesh.addEdge();
			e->edge = e0;
			e0->edge = e;
			pair->edge = e;
		}

		pair = getEdgePair(i1,i2,edgeMap);
		if(pair != NULL) {
			e1 = pair;
		} else {
			e1 = new hedge();
			pair = new hedge();
			e1->vertex = v2;
			pair->vertex = v1;
			e1->pair = pair;
			pair->pair = e1;
			v2->edge = e1;

			edgeMap[i1].push_back(e1);
			edgeMap[i2].push_back(pair);
			heedge * e = newMesh.addEdge();
			e->edge = e1;
			e1->edge = e;
			pair->edge = e;
		}

		pair = getEdgePair(i2,i0,edgeMap);
		if(pair != NULL) {
			e2 = pair;
		} else {
			e2 = new hedge();
			pair = new hedge();
			e2->vertex = v0;
			pair->vertex = v2;
			e2->pair = pair;
			pair->pair = e2;
			v0->edge = e2;


			edgeMap[i2].push_back(e2);
			edgeMap[i0].push_back(pair);
			heedge * e = newMesh.addEdge();
			e->edge = e2;
			e2->edge = e;
			pair->edge = e;
		}

		e0->next = e1;
		e1->next = e2;
		e2->next = e0;

		heface *f = newMesh.addFace();
		f->edge = e0;

		e2->face = e1->face = e0->face = f;

	}
	//do boundary
	for(auto e : newMesh.edges) {
		//breaks for non manifold meshes
		hedge * he = e->edge;
		if(he->pair->face == NULL) {
			he = he->pair;
		}
		if(he->face == NULL) {
			he->tag = 1;
			he->vertex->boundary = true;
			unsigned int startIndex = he->vertex->index;
			vector<hedge *> & map = edgeMap[startIndex];
			for(int i=0;i<map.size();++i) {
				hedge * he2 = map[i];
				if(he2->face == NULL) {
					he->next = he2;
					break;
				}
			}
			if(he->next == NULL) {
				cout << "ERROR NO NEXT BOUNDARY FOUND" << endl;
			}
		}
	}
	return newMesh;
}