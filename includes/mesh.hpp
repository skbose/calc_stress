#ifndef __mesh_hpp__
#define __mesh_hpp__

#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Geometry>

struct mesh_vertex
{
	Eigen::Vector3d pos;
	double area;
	double ratio;
	double neighborAreaDeformation;
	double neighborAreaOriginal;
	double ux, uy, uz;
	double lx, ly, lz;
	int neighbours;

	// constructor for the mesh vertex
	mesh_vertex(Eigen::Vector3d pos_) : pos(pos_) { area = 0.0; ux = uy = uz = 0.0; neighbours = 0; neighborAreaDeformation = 0.0; neighborAreaOriginal = 0.0;};
};

struct mesh_triangle
{
	// store vertex indices here
	int v_ind[3];
	double area;
	double deformedArea;

	mesh_triangle(int i, int j, int k, double area_)
	{
		v_ind[0] = i;
		v_ind[1] = j;
		v_ind[2] = k;

		area = area_;
		deformedArea = 0.0;
	}
};

class triangle_mesh_t
{
	public:
		std::vector<mesh_triangle *> triangles;
		std::vector<mesh_vertex *> vertices;

		// mesh destructor
		~triangle_mesh_t() { clear(); }

		// mesh constructor
		triangle_mesh_t(std::string const & path) { if (!load(path)) throw "Mesh: Could not load file"; calculatePtArea();}

		// load mesh data
		bool load(std::string const & path);

		// clear the mesh data
		void clear();

		std::string get_name(void) const { return std::string("mesh"); }

		void calculatePtArea();

		void calculateAreaDeformationOfMeshTriangles();

		void calculateLocalAreaDeformationAboutPoints();

		bool savePtArea(std::string const &path);

		bool saveMassRatio(std::string const & in_file_path, std::string const & out_file_path);

		bool saveMassRatioAsFeature(std::string const & path);

		bool loadDeformation(std::string const & path);

		bool loadDeformation(double * u);

		void flushVertexData();

		bool saveAbsLaplacianValuesAsFeature(std::string const & path);

		bool saveLaplacianVectorNormAsFeature(std::string const & path, int num_points=-1);

		bool saveLocalAreaDeformationAboutPointsAsFeature(std::string const &path, int max_points=-1);

		bool laplacianOnDeformations();

		bool savePointMovementMagnitudeAsFeature(std::string const &path, int max_points=-1);

		bool savePointAndDeformation(std::string const &path, int max_points=-1);

		// operator overloading to display the mesh deformation.
		void operator-(triangle_mesh_t const &m);
};

#endif