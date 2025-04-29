#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <GL/glew.h>
#include "myVector3D.h"

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

/*
void myMesh::checkMesh()
{	
	vector<myVertex *>::iterator it1;
	
	//Verify every vertices, every halfedges, every faces.
	vector<myHalfedge *>::iterator it2;
	for (it2 = halfedges.begin(); it2 != halfedges.end(); it2++)
	{
		if ((*it2)->twin == NULL)
			break;
		
			//if
	}
	if (it2 != halfedges.end())
		cout << "Error! Not all edges have their twins!\n";
	else cout << "Each edge has a twin!\n";

	vector<myFace *>::iterator it3;
}
*/

/*
vérifier en partant de origineof

pour chaque vertex, on calcul le vecteur moyen des normales des faces si il est bien égal au vecteur normal du sommet
    - on fait la somme de tous les vecteurs normaux divisé par le nombre de face
    - on vérifie si le vecteur normal du sommet est égal à celui calculé précédemment 

Demi-arêtes : 
    
Sommet : 

Face : Faire une boucle sur toutes les faces, et on fait la somme de tous les vecteurs normaux divisé par le nombre de face
*/

void myMesh::checkMesh()
{
    // Check for null pointers and closed circuits
    for (myHalfedge* he : halfedges) {
        if (!he || !he->next || !he->prev || !he->source || !he->adjacent_face) {
            cout << "Error: Null pointer found in half-edge structure!" << endl;
            return;
        }

        // Check for closed circuit (limit to 15 iterations)
        myHalfedge* temp = he;
        int count = 0;
        while (temp && count < 15) {
            temp = temp->next;
            count++;
        }
        if (count >= 15) {
            cout << "Error: Half-edge circuit is not closed!" << endl;
            return;
        }
    }

    // Check if adjacent half-edges are correctly connected to the face
    for (myFace* face : faces) {
        if (!face || !face->adjacent_halfedge) {
            cout << "Error: Null pointer found in face structure!" << endl;
            return;
        }

        myHalfedge* start = face->adjacent_halfedge;
        myHalfedge* current = start;
        do {
            if (current->adjacent_face != face) {
                cout << "Error: Half-edge not correctly connected to its face!" << endl;
                return;
            }
            current = current->next;
        } while (current && current != start);
    }

    // Check if the vertex normal matches the average of adjacent face normals
    for (myVertex* vertex : vertices) {
        if (!vertex || !vertex->point || !vertex->normal) {
            cout << "Error: Null pointer found in vertex structure!" << endl;
            return;
        }

        myVector3D averageNormal(0.0, 0.0, 0.0);
        myHalfedge* start = vertex->originof;
        myHalfedge* current = start;
        int faceCount = 0;

        if (!current) continue; // Skip if the vertex has no outgoing half-edge

        do {
            if (current->adjacent_face && current->adjacent_face->normal) {
                averageNormal.dX += current->adjacent_face->normal->dX;
                averageNormal.dY += current->adjacent_face->normal->dY;
                averageNormal.dZ += current->adjacent_face->normal->dZ;
                faceCount++;
            }
            current = current->twin ? current->twin->next : nullptr;
        } while (current && current != start);

        if (faceCount > 0) {
            averageNormal = averageNormal/faceCount;
            averageNormal.normalize();
            if (*vertex->normal != averageNormal) {
                cout << "Error: Vertex normal does not match the average of adjacent face normals!" << endl;
                return;
            }
        }
    }

    // Check if the sum of face normals is consistent
    myVector3D totalNormal(0.0, 0.0, 0.0);
    for (myFace* face : faces) {
        if (!face || !face->normal) {
            cout << "Error: Null pointer found in face structure!" << endl;
            return;
        }
        totalNormal.dX += face->normal->dX;
        totalNormal.dY += face->normal->dY;
        totalNormal.dZ += face->normal->dZ;
    }
    totalNormal = totalNormal / faces.size();
    totalNormal.normalize();

    cout << "Mesh check completed: No errors found!" << endl;
}


bool myMesh::readFile(std::string filename)
{
	string s, t, u;
	vector<int> faceids;
	myHalfedge **hedges;

	ifstream fin(filename);
	if (!fin.is_open()) {
		cout << "Unable to open file!\n";
		return false;
	}
	name = filename;

	map<pair<int, int>, myHalfedge *> twin_map;
	map<pair<int, int>, myHalfedge *>::iterator it;

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

			myPoint3D* point3D = new myPoint3D(x, y, z);
			myVertex* vertex = new myVertex;
			vertex->point = point3D;
			this->vertices.push_back(vertex);
			
		}
		else if (t == "mtllib") {}
		else if (t == "usemtl") {}
		else if (t == "s") {}
		else if (t == "f")
		{
			faceids.clear();
			while (myline >> u) // Read indices of vertices from a face
				faceids.push_back(atoi((u.substr(0, u.find("/"))).c_str()) - 1);
			if (faceids.size() < 3) // Ignore degenerate faces
				continue;

			hedges = new myHalfedge *[faceids.size()]; // Allocation
			for (unsigned int i = 0; i < faceids.size(); i++) 
				hedges[i] = new myHalfedge(); // Allocation

			myFace *f = new myFace(); // Allocation
			f->adjacent_halfedge = hedges[0]; // Connect the face with its incident edge

			for (unsigned int i = 0; i < faceids.size(); i++)
			{
				int iplusone = (i + 1) % faceids.size();
				int iminusone = (i - 1 + faceids.size()) % faceids.size();

				 // Connect prev and next half-edges
				hedges[i]->prev = hedges[iminusone];
				hedges[i]->next = hedges[iplusone];

				hedges[i]->source = vertices[faceids[i]];

				// Map the edge to create its twin
				pair<int, int> edge_key(faceids[i], faceids[iplusone]);
				pair<int, int> twin_key(faceids[iplusone], faceids[i]);

				auto it = twin_map.find(twin_key);
				if (it != twin_map.end()) {
					// Twin exists, assign it to the current half-edge
					hedges[i]->twin = it->second; //pair
					it->second->twin = hedges[i]; //pair
					twin_map.erase(it);
				} else {
					// Twin does not exist, add this edge to the map
					twin_map[edge_key] = hedges[i];
				}
				hedges[i]->adjacent_face = f;

				// Add the half-edge to the mesh
				halfedges.push_back(hedges[i]);
			}

			// Add the face to the mesh
			faces.push_back(f);
			delete[] hedges; // Free the allocated memory for half-edges
		}
	}

	checkMesh();
	normalize();

	return true;
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
	/**** TODO ****/
}

void myMesh::simplify(myVertex *)
{
	/**** TODO ****/
}

void myMesh::triangulate()
{   
    vector<myFace *> originalFaces = faces;

    for (myFace *f : originalFaces) {
        triangulate(f);
    }
}

bool myMesh::triangulate(myFace *f)
{
    myHalfedge *start = f->adjacent_halfedge;
    myHalfedge *current = start->next->next;

    // Return false if the face is already a triangle
    if (current->next == start) {
        return false;
    }

    // Split faces
    while (current != start) {
        myFace *newFace = new myFace();
        myHalfedge *he1 = new myHalfedge();
        myHalfedge *he2 = new myHalfedge();
        myHalfedge *he3 = new myHalfedge();

        he1->source = start->source;
        he2->source = current->prev->source;
        he3->source = current->source;

        he1->next = he2; he2->next = he3; he3->next = he1;
        he1->prev = he3; he2->prev = he1; he3->prev = he2;

        he1->adjacent_face = newFace;
        he2->adjacent_face = newFace;
        he3->adjacent_face = newFace;

        newFace->adjacent_halfedge = he1;

        // Add new objects to the mesh
        faces.push_back(newFace);
        halfedges.push_back(he1);
        halfedges.push_back(he2);
        halfedges.push_back(he3);

        current = current->next;
    }

    return true;
}

