#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <GL/glew.h>
#include "myVector3D.h"
#include <algorithm>

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
    /**** TODO ****/
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
			while (myline >> u)
				faceids.push_back(atoi((u.substr(0, u.find("/"))).c_str()) - 1);

			if (faceids.size() < 3)
				continue;

			hedges = new myHalfedge * [faceids.size()];
			for (unsigned int i = 0; i < faceids.size(); i++)
				hedges[i] = new myHalfedge();

			myFace* f = new myFace();
			f->adjacent_halfedge = hedges[0];

			for (unsigned int i = 0; i < faceids.size(); i++)
			{
				int nextindex = (i + 1) % faceids.size();
				int previousindex = (i - 1 + faceids.size()) % faceids.size();

				//connexion of halfedges
				hedges[i]->next = hedges[nextindex];
				hedges[i]->prev = hedges[previousindex];

				//attribution of faces
				hedges[i]->adjacent_face = f;

				//definitions origin of the halfedges
				if (faceids[i] >= vertices.size()) continue;
				hedges[i]->source = vertices[faceids[i]];
				vertices[faceids[i]]->originof = hedges[i];

				pair<int, int> edge_key(faceids[i], faceids[nextindex]);
				pair<int, int> twin_key(faceids[nextindex], faceids[i]);
				twin_map[edge_key] = hedges[i];

				if (twin_map.find(twin_key) != twin_map.end())
				{
					hedges[i]->twin = twin_map[twin_key];
					twin_map[twin_key]->twin = hedges[i];
				}

				//add halfedges to mesh
				halfedges.push_back(hedges[i]);
			}

			//add faces to mesh
			faces.push_back(f);
		}

	}
	checkMesh();
	normalize();
    logMeshStatistics();

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
    std::vector<myFace*> original_faces = faces;

    for (myFace* f : original_faces)
    {
        triangulate(f);
    }

    logMeshStatistics();
}



bool myMesh::triangulate(myFace* f)
{
    //recupere chaque halfedges
    std::vector<myHalfedge*> hedge_list;
    myHalfedge* start = f->adjacent_halfedge;
    myHalfedge* current = start;

    do {
        hedge_list.push_back(current);
        current = current->next;
    } while (current != start);

    int n = hedge_list.size();
    if (n == 3) return false; //si un triangle passe

    //premier sommet point fixe
    myVertex* v0 = hedge_list[0]->source;

    auto it = std::find(faces.begin(), faces.end(), f);
    if (it != faces.end()) faces.erase(it);

    delete f;

    //cration triangle
    for (int i = 1; i < n - 1; ++i)
    {
        myFace* new_face = new myFace();

        //creation halfedges
        myHalfedge* he1 = new myHalfedge(); // v0 -> vi
        myHalfedge* he2 = new myHalfedge(); // vi -> vi+1
        myHalfedge* he3 = new myHalfedge(); // vi+1 -> v0

        he1->source = v0;
        he2->source = hedge_list[i]->source;
        he3->source = hedge_list[i + 1]->source;

        //connexion halfedges
        he1->next = he2; he2->next = he3; he3->next = he1;
        he1->prev = he3; he2->prev = he1; he3->prev = he2;

        //association halfedges a leur face
        he1->adjacent_face = new_face;

        if (!he1->source->originof) he1->source->originof = he1;

        new_face->adjacent_halfedge = he1;

        //ajoutdans mesh
        faces.push_back(new_face);
        halfedges.push_back(he3);
    }

    return true;
}

void myMesh::logMeshStatistics() {
    int triangleCount = 0;
    int squareCount = 0;
    int otherPolygonCount = 0;

    for (myFace* face : faces) {
        if (!face || !face->adjacent_halfedge) continue;

        // Count the number of vertices in the face
        int vertexCount = 0;
        myHalfedge* start = face->adjacent_halfedge;
        myHalfedge* current = start;
        do {
            vertexCount++;
            current = current->next;
        } while (current != start);

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