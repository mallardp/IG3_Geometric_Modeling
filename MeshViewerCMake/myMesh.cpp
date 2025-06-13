#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <GL/glew.h>
#include "myVector3D.h"
#include <algorithm>
#include <queue>
#include <set>

using namespace std;

myMesh::myMesh(void)
{
	/**** TODO ****/
}


myMesh::~myMesh(void)
{
	/**** TODO ****/
}

void myMesh::clear()
{
	for (unsigned int i = 0; i < vertices.size(); i++) if (vertices[i]) delete vertices[i];
	for (unsigned int i = 0; i < halfedges.size(); i++) if (halfedges[i]) delete halfedges[i];
	for (unsigned int i = 0; i < faces.size(); i++) if (faces[i]) delete faces[i];

	vector<myVertex *> empty_vertices;    vertices.swap(empty_vertices);
	vector<myHalfedge *> empty_halfedges; halfedges.swap(empty_halfedges);
	vector<myFace *> empty_faces;         faces.swap(empty_faces);
}

void myMesh::checkMesh()
{	
	//Verify every vertices, every halfedges, every faces.

	cout << "\nChecking Vertices :" << endl;
	checkVerticesNull();
	checkCircuitAroundVertex();
	
	cout << "\nChecking Halfedges :" << endl;
	checkHalfedgesNull();
	checkCircuitHalfedge();

	cout << "\n Checking Faces :" << endl;
	checkFacesNull();
	checkCircuitAroundFace();
}

void myMesh::checkVerticesNull()
{
	vector<myVertex *>::iterator it;
	for (it = vertices.begin(); it != vertices.end(); it++)
	{
		if ((*it)->originof == NULL)
			throw std::runtime_error("Error! Not all vertices have their origin halfedge!");
		if ((*it)->point == NULL)
			throw std::runtime_error("Error! Not all vertices have their point!");
		if ((*it)->normal == NULL)
			throw std::runtime_error("Error! Not all vertices have their normal!");
		
		if ((*it)->originof->source != *it)
			throw std::runtime_error("Error! The halfedge around the vertex doesn't have the same source!");
	}
	cout << "All vertices have their originof halfedge, point and normal!" << endl;
}

void myMesh::checkHalfedgesNull()
{
	vector<myHalfedge *>::iterator it;
	for (it = halfedges.begin(); it != halfedges.end(); it++)
	{
		if ((*it)->twin == NULL)
			throw std::runtime_error("Error! Not all edges have their twins!");
		else if ((*it)->twin->twin != *it)
			throw std::runtime_error("Error! The twins are not reciprocal!");

		if ((*it)->source == NULL)
			throw std::runtime_error("Error! Not all edges have their source vertex!");
		if ((*it)->adjacent_face == NULL)
			throw std::runtime_error("Error! Not all edges have their adjacent face!");
		if ((*it)->next == NULL)
			throw std::runtime_error("Error! Not all edges have their next halfedge!");
		if ((*it)->prev == NULL)
			throw std::runtime_error("Error! Not all edges have their previous halfedge!");
	}
	cout << "All halfedges have twins, source, adjacent faces, next and prev he!" << endl;
}

void myMesh::checkFacesNull()
{
	vector<myFace *>::iterator it;
	for (it = faces.begin(); it != faces.end(); it++)
	{
		if ((*it)->adjacent_halfedge == NULL)
			throw std::runtime_error("Error! Not all faces have their adjacent halfedge!");
		if ((*it)->normal == NULL)
			throw std::runtime_error("Error! Not all faces have their normals!");
	}
	cout << "All faces have their adjacent halfedge and normals!" << endl;
}

void myMesh::checkCircuitAroundVertex()
{	
	vector<myVertex *>::iterator it;
	for (it = vertices.begin(); it != vertices.end(); it++)
	{
		myHalfedge* start = (*it)->originof;
		myHalfedge* heNext = start;
		myHalfedge* hePrev = start;
		int count = 0;

		do {
			count++;
			heNext = heNext->twin->next;

			if (heNext->source != (*it)) {
				throw std::runtime_error("Error! Halfedges around the vertex don't have the same source!");
			}
			if (count > 35) {
				throw std::runtime_error("Error! Too many halfedges around the vertex!");
			}
		} while (heNext != start);

		if (heNext != start) {
			throw std::runtime_error("Error! Halfedges twin->next circuit is not closed!");
		}
	}
	cout << "All halfedges around vertices have the same source vertex!" << endl;
}

void myMesh::checkCircuitHalfedge()
{
	vector<myHalfedge *>::iterator it;
	for (it = halfedges.begin(); it != halfedges.end(); it++)
	{
		myHalfedge* start = *it;
		myHalfedge* heNext = start;
		myFace* f = (*it)->adjacent_face;
		int count = 0;
		
		do{
			count++;
			heNext = heNext->next;
			if (heNext->adjacent_face != f) {
				throw std::runtime_error("Error! Halfedge next circuit is not around the same face!");
			}
			if (count > 15) {
				throw std::runtime_error("Error! Halfedge circuit is not closed (> 15 edges)!");
			}

		}
		while (heNext != start);

		if (heNext != start) {
			throw std::runtime_error("Error! Halfedge next circuit is not closed!");
		}
	}
	cout << "All halfedges are closed circuits!" << endl;
}

void myMesh::checkCircuitAroundFace()
{	
	vector<myFace *>::iterator it;
	for (it = faces.begin(); it != faces.end(); it++)
	{
		myHalfedge* start = (*it)->adjacent_halfedge;
		myHalfedge* he = start;
		int count = 0;

		do {
			count++;
			he = he->next;
			if (count > 25) {
				throw std::runtime_error("Error! Too many halfedges around the face!");
			}
		} while (he != start);

		if (he != start) {
			throw std::runtime_error("Error! Halfedges circuit is not closed!");
		}
	}
	cout << "All faces are closed circuits!" << endl;
}

bool myMesh::readFile(std::string filename)
{
	string s, t, u;
	vector<int> faceids;
	myHalfedge** hedges;

	ifstream fin(filename);
	if (!fin.is_open()) {
		cout << "Unable to open file!\n";
		return false;
	}
	
	name = filename;

	map<pair<int, int>, myHalfedge*> twin_map;
	map<pair<int, int>, myHalfedge*>::iterator it;

	myFace* f = nullptr;

	while (getline(fin, s))
	{
		stringstream myline(s);
		myline >> t;

		if (t == "g") {}
		else if (t == "v")
		{
			float x, y, z;
			myline >> x >> y >> z;
			cout << "v " << x << " " << y << " " << z << endl;
			myVertex* v = new myVertex();
			v->point = new myPoint3D(x, y, z);
			vertices.push_back(v);
		}
		else if (t == "mtllib") {}
		else if (t == "usemtl") {}
		else if (t == "s") {}
		else if (t == "f")
		{
			faceids.clear();
			while (myline >> u) // read indices of vertices from a face into a container - it helps to access them later 
				faceids.push_back(atoi((u.substr(0, u.find("/"))).c_str()) - 1);
			if (faceids.size() < 3) // ignore degenerate faces
				continue;

			hedges = new myHalfedge * [faceids.size()]; // allocate the array for storing pointers to half-edges
			for (unsigned int i = 0; i < faceids.size(); i++)
				hedges[i] = new myHalfedge(); // pre-allocate new half-edges

			f = new myFace(); // allocate the new face
			f->adjacent_halfedge = hedges[0]; // connect the face with incident edge

			for (unsigned int i = 0; i < faceids.size(); i++)
			{
				int iplusone = (i + 1) % faceids.size();
				int iminusone = (i - 1 + faceids.size()) % faceids.size();

				// YOUR CODE COMES HERE!

				// connect prevs, and next
				hedges[i]->next = hedges[iplusone];
				hedges[i]->prev = hedges[iminusone];

				hedges[i]->adjacent_face = f;

				// search for the twins using twin_map
				pair<int, int> edge_key(faceids[i], faceids[iplusone]);
				it = twin_map.find(edge_key);
				if (it != twin_map.end()) {
					hedges[i]->twin = it->second;
					it->second->twin = hedges[i];
				}
				else {
					twin_map[make_pair(faceids[iplusone], faceids[i])] = hedges[i];
				}

				// set originof
				hedges[i]->source = vertices[faceids[i]];
				vertices[faceids[i]]->originof = hedges[i];

				// push edges to halfedges in myMesh
				halfedges.push_back(hedges[i]);
			}
			// push faces to faces in myMesh
			faces.push_back(f);
		}
	}
	

	checkNonManifoldEdges();
	checkDisconnectedComponents();
	checkMesh();
	normalize();

	logMeshStatistics();

	return true;
}

void myMesh::checkNonManifoldEdges() {
    map<pair<int, int>, int> edgeCount;

    for (myHalfedge* he : halfedges) {
        // Find the indices of the source and next vertices
        auto it1 = find(vertices.begin(), vertices.end(), he->source);
        auto it2 = find(vertices.begin(), vertices.end(), he->next->source);

        if (it1 == vertices.end() || it2 == vertices.end()) {
            throw std::runtime_error("Error! Halfedge source or next vertex not found in vertices list!");
        }

        int v1 = distance(vertices.begin(), it1);
        int v2 = distance(vertices.begin(), it2);

        // Create an edge key with sorted vertex indices
        pair<int, int> edgeKey = make_pair(min(v1, v2), max(v1, v2));
        edgeCount[edgeKey]++;
    }

    for (auto& entry : edgeCount) {
        if (entry.second > 2) {
            cout << "Non-manifold edge detected between vertices " << entry.first.first
                 << " and " << entry.first.second << " (shared by " << entry.second << " faces)" << endl;
        }
    }
	cout << "Non-manifold edge check completed." << endl;
}

void myMesh::checkDisconnectedComponents() {
    set<myVertex*> visited;
    queue<myVertex*> toVisit;

    if (!vertices.empty()) {
        toVisit.push(vertices[0]);
        visited.insert(vertices[0]);
    }

    while (!toVisit.empty()) {
        myVertex* v = toVisit.front();
        toVisit.pop();

        myHalfedge* start = v->originof;
        myHalfedge* he = start;

        do {
            myVertex* neighbor = he->next->source;
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                toVisit.push(neighbor);
            }
            he = he->twin->next;
        } while (he != start);
    }

    if (visited.size() != vertices.size()) {
        cout << "Disconnected components detected! Visited " << visited.size()
             << " out of " << vertices.size() << " vertices." << endl;
    } else {
        cout << "Mesh is fully connected." << endl;
    }
}

void myMesh::computeNormals() {
    // Compute normals for all faces
    for (myFace* face : faces) {
        if (face && face->adjacent_halfedge) {
            face->computeNormal();
        }
    }

    // Initialize vertex normals
    for (myVertex* vertex : vertices) {
        if (vertex && vertex->normal) {
            vertex->normal->dX = 0.0;
            vertex->normal->dY = 0.0;
            vertex->normal->dZ = 0.0;
        }
    }

    // Accumulate face normals for each vertex
    for (myHalfedge* he : halfedges) {
        if (he && he->source && he->adjacent_face && he->adjacent_face->normal) {
            he->source->normal->dX += he->adjacent_face->normal->dX;
            he->source->normal->dY += he->adjacent_face->normal->dY;
            he->source->normal->dZ += he->adjacent_face->normal->dZ;
        }
    }

    // Normalize the normals
    for (myVertex* vertex : vertices) {
        if (vertex && vertex->normal) {
            vertex->normal->normalize();
        }
    }
}

void myMesh::normalize()
{
	if (vertices.size() < 1) return;

	int tmpxmin = 0, tmpymin = 0, tmpzmin = 0, tmpxmax = 0, tmpymax = 0, tmpzmax = 0;

	for (unsigned int i = 0; i < vertices.size(); i++) {
		if (vertices[i]->point->X < vertices[tmpxmin]->point->X) tmpxmin = i;
		if (vertices[i]->point->X > vertices[tmpxmax]->point->X) tmpxmax = i;

		if (vertices[i]->point->Y < vertices[tmpymin]->point->Y) tmpymin = i;
		if (vertices[i]->point->Y > vertices[tmpymax]->point->Y) tmpymax = i;

		if (vertices[i]->point->Z < vertices[tmpzmin]->point->Z) tmpzmin = i;
		if (vertices[i]->point->Z > vertices[tmpzmax]->point->Z) tmpzmax = i;
	}

	double xmin = vertices[tmpxmin]->point->X, xmax = vertices[tmpxmax]->point->X,
		ymin = vertices[tmpymin]->point->Y, ymax = vertices[tmpymax]->point->Y,
		zmin = vertices[tmpzmin]->point->Z, zmax = vertices[tmpzmax]->point->Z;

	double scale = (xmax - xmin) > (ymax - ymin) ? (xmax - xmin) : (ymax - ymin);
	scale = scale > (zmax - zmin) ? scale : (zmax - zmin);

	for (unsigned int i = 0; i < vertices.size(); i++) {
		vertices[i]->point->X -= (xmax + xmin) / 2;
		vertices[i]->point->Y -= (ymax + ymin) / 2;
		vertices[i]->point->Z -= (zmax + zmin) / 2;

		vertices[i]->point->X /= scale;
		vertices[i]->point->Y /= scale;
		vertices[i]->point->Z /= scale;
	}
}


void myMesh::splitFaceTRIS(myFace *f, myPoint3D *p)
{
	/**** TODO ****/
}

void myMesh::splitEdge(myHalfedge *e1, myPoint3D *p)
{

	/**** TODO ****/
}

void myMesh::splitFaceQUADS(myFace *f, myPoint3D *p)
{
	/**** TODO ****/
}


void myMesh::subdivisionCatmullClark()
{
	/**** TODO ****/
}

void myMesh::simplify()
{	
	if (faces.size() < 3) {
		cout << "Mesh has faces with less than 3 halfedges." << endl;
		cout << "Cannot simplify further." << endl;
		return;
	}

	// 1) Calculate lengths of all halfedges
	for (myHalfedge* he : halfedges) {
		if (!he || !he->source || !he->source->point ||
			!he->twin || !he->twin->source || !he->twin->source->point) {
			cout << "Invalid halfedge detected." << endl;
			return;
		}
		he->calculateLength();
	}

	// 2) Find the shortest halfedge
	myHalfedge* shortest = nullptr;
	for (myHalfedge* he : halfedges) {
		if(shortest == nullptr || he->length < shortest->length)
		{
			shortest = he;
		}
	}
	if (!shortest || !shortest->twin) {
        cout << "Invalid shortest halfedge." << endl;
        return;
    }

	myHalfedge* he1 = shortest;
	myHalfedge* he2 = shortest->twin;
	myVertex* vA = he1->source;
	myVertex* vB = he2->source;

	// 3) Create a new vertex at the midpoint of vA and vB
	myVertex* middleVertex = new myVertex();
	middleVertex->point = new myPoint3D(
		(vA->point->X + vB->point->X) / 2,
		(vA->point->Y + vB->point->Y) / 2,
		(vA->point->Z + vB->point->Z) / 2
	);
	middleVertex->originof = he1->prev->twin;
	
	// 4) Update sources of halfedges around vertex A&B
	auto updateSources = [middleVertex](myVertex* v) {
        if (!v || !v->originof) return;
        myHalfedge* start = v->originof;
        myHalfedge* he = start;
        do {
            if (he) he->source = middleVertex;
            he = he->twin ? he->twin->next : nullptr;
        } while (he && he != start);
    };
    updateSources(vA);
    updateSources(vB);

	// Function to remove a halfedge
	auto removeHalfedge = [this](myHalfedge* he) {
		auto it = std::find(halfedges.begin(), halfedges.end(), he);
		halfedges.erase(it);
		delete he;
	};

	// Function to remove a vertex
    auto removeVertex = [this](myVertex* v) {
        auto it = std::find(vertices.begin(), vertices.end(), v);
        vertices.erase(it);
		delete v;
    };

	// Function to remove a triangle face
    auto removeTriangleFace = [this, removeHalfedge](myHalfedge* he) {
		myHalfedge* he1 = he;
		myHalfedge* he2 = he1->next;
		myHalfedge* he3 = he2->next;

		myVertex* v3 = he3->source;

		if (v3->originof == he3) {
			v3->originof = he3->twin->next;
		}

		myHalfedge* twin1 = he1->twin;
		myHalfedge* twin2 = he2->twin;
		myHalfedge* twin3 = he3->twin;

		twin2->twin = twin3;
		twin3->twin = twin2;

        auto it = std::find(faces.begin(), faces.end(), he->adjacent_face);
        faces.erase(it);
        delete he->adjacent_face;

		removeHalfedge(he1);
		removeHalfedge(he2);
		removeHalfedge(he3);
    };

	// Function to simplify a non-triangle face
    auto updateNonTriangleFace = [removeHalfedge](myHalfedge* he) {
        if (he->adjacent_face->adjacent_halfedge == he) {
            he->adjacent_face->adjacent_halfedge = he->next;
        }
        he->next->prev = he->prev;
        he->prev->next = he->next;
		removeHalfedge(he);
    };
	

	// Check if shortest->adjacent_face is a triangle
	if (he1->adjacent_face && he1->adjacent_face->countEdges() == 3) {
		removeTriangleFace(he1);
	} else {
		updateNonTriangleFace(he1);
	}

	// Check if shortest->twin->adjacent_face is a triangle
	if (he2->adjacent_face && he2->adjacent_face->countEdges() == 3) {
		removeTriangleFace(he2);
	} else {
		updateNonTriangleFace(he2);
	}

    vertices.push_back(middleVertex);    
	removeVertex(vA);
    removeVertex(vB);

	checkMesh();
	logMeshStatistics();
}

void myMesh::triangulate()
{
    int length = faces.size();
	for (int i = 0; i < length; i++) {
		myFace* f = faces[i];
		triangulate(f);
	}

    logMeshStatistics();
}

bool myMesh::triangulate(myFace* f) {
    // No triangulation needed
    if (f->countEdges() == 3) {
        return false;
    }

    myHalfedge* start = f->adjacent_halfedge;
    myVertex* commonSource = start->source;
    myHalfedge* current = start;

    // Keep triangulating until the remaining face has only 3 vertices
    while (current->next->next->next != current) {
        // Diagonal and twin halfedges
        myHalfedge* he1 = new myHalfedge();
        myHalfedge* he2 = new myHalfedge();

        he1->twin = he2;
        he2->twin = he1;

        // Define the direction of each halfedge
        he1->source = current->next->next->source;
        he2->source = commonSource;

        // Create a new face (triangle)
        myFace* newFace = new myFace();
        newFace->adjacent_halfedge = he2;

        // Save the next and previous edges for he2 (triangle side)
        myHalfedge* he2Next = current->next->next;
        myHalfedge* he2Prev = current->prev;

        // Link the new diagonal (he1)
        he1->next = current;
        he1->prev = current->next;
        current->next->next = he1;
        current->prev = he1;

        // Link the twin halfedge (he2)
        he2->next = he2Next;
        he2->prev = he2Prev;
        he2Next->prev = he2;
        he2Prev->next = he2;

        // Assign the new face to its 3 halfedges
        myHalfedge* cur = he2;
        do {
            cur->adjacent_face = newFace;
            cur = cur->next;
        } while (cur != he2);

        // Assign the current face to the diagonal
        he1->adjacent_face = f;


        faces.push_back(newFace);
        halfedges.push_back(he1);
        halfedges.push_back(he2);

        // Move on to the remaining polygon
        current = he2;
    }

    return true;
}

void myMesh::logMeshStatistics() {
    int triangleCount = 0;
    int squareCount = 0;
    int otherPolygonCount = 0;

    for (myFace* face : faces) {
        if (!face || !face->adjacent_halfedge) continue;

        int vertexCount = face->countEdges();

        // Classify the face based on the vertex count
        if (vertexCount == 3) {
            triangleCount++;
        } else if (vertexCount == 4) {
            squareCount++;
        } else {
            otherPolygonCount++;
        }
    }

    // Print the counts
    cout << endl;
    cout << "Mesh statistics:" << endl;
    cout << "  Triangles: " << triangleCount << endl;
    cout << "  Squares: " << squareCount << endl;
    cout << "  Other polygons: " << otherPolygonCount << endl;
}