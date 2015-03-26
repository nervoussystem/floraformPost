#include "DiscreteExteriorCalculus.h"

   void HodgeStar0Form :: build( const hemesh& mesh,
                                    SparseMatrix& star0 )
   // builds a diagonal matrix mapping primal discrete 0-forms
   // to dual discrete 2-forms
   {
      int nV = mesh.vertices.size();
	  vector<Triplet> trips;
	  star0.resize(nV,nV);
	  
	
	  for(auto v : mesh.vertices) {
		  int i = v->index;
		  trips.push_back(Triplet(i,i,v->area()));
		  //star0.diagonal[i] = v->area();
	  }
	  star0.setFromTriplets(trips.begin(),trips.end());
   }

   void HodgeStar1Form :: build( const hemesh& mesh,
                                    SparseMatrix& star1 )
   // builds a diagonal matrix mapping primal discrete 1-forms
   // to dual discrete 1-forms
   {
      int nE = mesh.edges.size();

      star1.resize(nE,nE);
	  vector<Triplet> trips;

      for(auto e : mesh.edges )
      {
         // get the cotangents of the two angles opposite this edge
         double cotAlpha = e->edge->cotan();
         double cotBeta  = e->edge->pair->cotan();

         int i = e->index;
		 trips.push_back(Triplet(i,i,(cotAlpha+cotBeta)*0.5));
		 //star1.diagonal[i] = (cotAlpha+cotBeta)*0.5;
      }
	  star1.setFromTriplets(trips.begin(),trips.end());
   }

   void HodgeStar2Form :: build( const hemesh& mesh,
                                    SparseMatrix& star2 )
   // builds a diagonal matrix mapping primal discrete 2-forms
   // to dual discrete 2-forms
   {
      int nF = mesh.faces.size();

	  star2.resize(nF,nF);
	  vector<Triplet> trips;

	  for( auto f : mesh.faces) {
		  int i= f->index;
		  trips.push_back(Triplet(i,i,1.0/f->area));
		  //star2.diagonal[i] = 1.0/f->area;
	  }
	  star2.setFromTriplets(trips.begin(),trips.end());
	  
   }

   void ExteriorDerivative0Form :: build( const hemesh& mesh,
                                             SparseMatrix& d0 )
   {
      int nV = mesh.vertices.size();
      int nE = mesh.edges.size();

	  d0.resize(nE,nV);
	  vector<Triplet> trips;

	  for(auto e : mesh.edges) {
		  int r = e->index;
		  
		  int ci = e->edge->vertex->index;
		  int cj = e->edge->pair->vertex->index;
		  trips.push_back(Triplet(r,ci,1.0));
		  trips.push_back(Triplet(r,cj,-1.0));
	  }
	  d0.setFromTriplets(trips.begin(),trips.end());
   }

   void ExteriorDerivative1Form :: build( const hemesh& mesh,
                                             SparseMatrix& d1 )
   {
      int nE = mesh.edges.size();
      int nF = mesh.faces.size();

	  d1.resize(nF,nE);
	  vector<Triplet> trips;

      // visit each face
	  for( auto f : mesh.faces) {
         int r = f->index;
		 hedge * he = f->edge;
		 do {
			 int c = he->edge->index;
			 double s = he->edge->edge == he ? 1.0 : -1.0;
			 trips.push_back(Triplet(r,c,s));
			 he = he->next;
		 } while(he != f->edge);
	  }
   }


