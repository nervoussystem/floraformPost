#include "solver.h"


static SpMat solverMatrix;
static Cg linearSolver;
static Vectorf dVelocities;
static ofVec3f * positions;
static ofVec3f * velocities;
static Vectorf forces;
static Eigen::BiCGSTAB<SpMat> solver;
static SpMat IDENTITY;
static int numVertices = 0;
static vector< Eigen::Triplet<double> > triplets;

float maxEdgeLength;
float maxEdgeLengthSq;
vector<float> distances;
vector<float> distances2;


void setupSolver(hemesh & hmesh) {
	//edges
	for(int i=0;i<hmesh.edges.size();++i) {
		heedge * e = hmesh.edges[i];
		hedge *he = e->edge;
		hedge *hep = he->pair;
		hevertex *p1 = he->vertex;
		hevertex *p2 = hep->vertex;
		ofVec3f v1 = hmesh.mesh->getVertex(p1->index);//p1->position;
		ofVec3f v2 = hmesh.mesh->getVertex(p2->index);//p2->position;
		ofVec3f dir = v1-v2;
		e->length = dir.length();
		//e->invLength = 1.0/e->length;
		float invLength = 1.0/e->length;
		dir *= invLength;
		he->dir = dir;
		hep->dir = -dir;
	}
	//faces
	for(auto f : hmesh.faces) {
		hedge * he = f->edge;
		ofVec3f p0 = hmesh.mesh->getVertex(he->vertex->index);
		hedge * he1 = he->next;
		ofVec3f p1 = hmesh.mesh->getVertex(he1->vertex->index);
		hedge * he2 = he1->next;
		ofVec3f p2 = hmesh.mesh->getVertex(he2->vertex->index);

		ofVec3f v0 = p1-p0;
		ofVec3f v1 = p2-p0;

		f->norm = v0.getCrossed(v1);
		f->area = f->norm.length();
		f->norm /= f->area;
		f->area *= 0.5;
		
		he->cosA = -(he1->dir.dot(he2->dir));
		he1->cosA = -(he->dir.dot(he2->dir));
		he2->cosA = -(he1->dir.dot(he->dir));

		he->norm = he->dir.getCrossed(f->norm);
		he1->norm = he1->dir.getCrossed(f->norm);
		he2->norm = he2->dir.getCrossed(f->norm);

		he->invAlt = he->edge->length*0.5/f->area;
		he1->invAlt = he1->edge->length*0.5/f->area;
		he2->invAlt = he2->edge->length*0.5/f->area;
	}
}

void setDistanceBoundary(hemesh & hmesh, Vectord & u, bool tag) {
	u.setZero();
	/*
	for(auto e : hmesh.edges) {
		hedge *he = e->edge;
		hedge *hep = he->pair;
		
		if(he->face == NULL) {
			if(tag || he->tag) {
				//check edge normal
				ofVec3f v = hmesh.mesh->getVertex(he->vertex->index);

				//if(v.z > 60.3) {
					u[he->vertex->index] = 1.0;
				//}
			}
		} else if(hep->face == NULL) {
			if(tag || hep->tag) {
				ofVec3f v = hmesh.mesh->getVertex(hep->vertex->index);

				//if(v.z > 60.3) {
					u[hep->vertex->index] = 1.0;
				//}
			}
		}
	}
	*/
	for(auto v : hmesh.vertices) {
		if((tag && v->boundary) || v->tag) {
			u[v->index] = 1.0;
		}
	}
	/*
	for(int i=0;i<10;++i) {
		int index = (int)ofRandom(hmesh.vertices.size());
		u[index] = 1.0;
	}*/
}

void computeVectorField(const Vectord & u, hemesh & mesh)
{
	for(auto f : mesh.faces) {
		hedge * hij = f->edge;
		hedge * hjk = hij->next;
		hedge * hki = hjk->next;
		hevertex * vj = hij->vertex;
		hevertex * vk = hjk->vertex;
		hevertex * vi = hki->vertex;

		float ui = u[vi->index];
		float uj = u[vj->index];
		float uk = u[vk->index];

		ofVec3d eij = hij->norm*hij->edge->length;
		ofVec3d ejk = hjk->norm*hjk->edge->length;
		ofVec3d eki = hki->norm*hki->edge->length;

		ofVec3d X = 0.5f * ( ui*ejk + uj*eki + uk*eij ) / f->area;
		X.normalize();
		f->vector = - X;
		//f->vector = f->vector.rotate(85,f->norm);
	}

}
      
void computeDivergence(const hemesh& mesh, Vectord & div)
{
	for(auto v : mesh.vertices) {
		double sum = 0.0;
		hedge * he = v->edge;
		do {
			if(he->face != NULL) {
				hedge * opp = he->next->next;
				ofVec3d n = opp->norm*opp->edge->length;
				ofVec3d v = he->face->vector;
				sum += n.dot(v);
			}
			he = he->next->pair;
		} while(he != v->edge);

		div[v->index] = sum;
	}
}

struct PruneFunc {
	int i;
	bool operator()(int r, int c, float f) const 
	{
		return (r!=i && c!=i) || (r==i &&c ==i);
	}
};

void geodesicDistance(hemesh &hmesh, vector<float> & out, bool tag) {
	//set initialize conditions
	int nV = hmesh.vertices.size();
	out.resize(nV);
	Vectord u0(nV);
	setDistanceBoundary(hmesh, u0,tag);
         
    // DEC
    SparseMatrix star0;
    HodgeStar0Form::build( hmesh, star0 );
         
    SparseMatrix  star1;
    HodgeStar1Form::build( hmesh, star1 );
         
    SparseMatrix d0;
    ExteriorDerivative0Form::build( hmesh, d0 );
         
    // zero Neumann boundary condition
    SparseMatrix L =  d0.transpose() * star1 * d0;
    
    // make L positive-definite
	//changed from 1.0e-8 to 1.0e-6 and it magically worked
    L += 1.0e-8*star0;
  
    // heat flow for short interval
    float step = 2.0*2.f;//0.6;
    SparseMatrix A = star0 + step * L;

    Vectord u;
	SparseSolver solver(A);
	solver.setTolerance(1e-22);
	solver.setMaxIterations(1000);

	u = solver.solve(u0);
    //solvePositiveDefinite(A, u, u0);

    // extract geodesic
    computeVectorField(u, hmesh);
    
    Vectord div(nV);
    computeDivergence(hmesh, div);
	//Eigen::Map<Vectorf> phi(out.data(),nV);
	Vectord phi(nV);
	//boundary conditions
	
	PruneFunc func;
	for(int i=0;i<nV;++i) {
		if(u0[i] > .5) {
			func.i = i;
			L.prune(func);
			L.coeffRef(i,i) = 1.0;
			div[i] = 0.0;
		}
	}
	
	SparseSolver solver2(L);
	solver2.setTolerance(1e-20);
	//solver.compute(L);
	phi = solver2.solve(div);
	cout << "div err " << solver2.error() << " I " << solver2.iterations() << endl;
    //solvePositiveDefinite(L, phi, div);
	//phi = A*phi;
	double minV = 1e20;
	for(int i=0;i<nV;++i) {
		minV = min(minV, phi[i]);
	}
	for(int i=0;i<nV;++i) {
		phi[i] -= minV;
	}

	for(int i=0;i<nV;++i) {
		out[i] = phi[i];
	}

    //setMinToZero(phi);
    //assignDistance(phi, mesh);         
    //return phi.norm();
}