#pragma once
#include "myFace.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <vector>
#include <string>

class myMesh
{
public:
	std::vector<myVertex *> vertices;
	std::vector<myHalfedge *> halfedges;
	std::vector<myFace *> faces;
	std::string name;

	void checkMesh();

	void checkVerticesNull();
	void checkHalfedgesNull();
	void checkFacesNull();

	void checkCircuitAroundVertex();
	void checkCircuitHalfedge();
	void checkCircuitAroundFace();


	bool readFile(std::string filename);
	void checkNonManifoldEdges();
	void checkDisconnectedComponents();
	void computeNormals();
	void normalize();

	void subdivisionCatmullClark();

	void splitFaceTRIS(myFace *, myPoint3D *);

	void splitEdge(myHalfedge *, myPoint3D *);
	void splitFaceQUADS(myFace *, myPoint3D *);

	void triangulate();
	bool triangulate(myFace *);
	bool triangulateF(myFace *);
	void simplify();
	void simplify(myVertex *);
	void simplifyOLD();
	void logMeshStatistics();

	void clear();

	myMesh(void);
	~myMesh(void);
};

