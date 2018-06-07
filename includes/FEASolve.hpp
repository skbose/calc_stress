#ifndef __FEASolve__hpp__
#define __FEASolve__hpp__

#ifndef PARDISO_SOLVER_IS_AVAILABLE
#endif

#define PARDISO_SOLVER_IS_AVAILABLE
#define NUM_COMPUTE_THREADS 1

#include "VegaHeaders.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <ctime>
#include <Eigen/Dense>

#include "SimulatorApp.hpp"
#include "mesh.hpp"


using namespace std;

class FEASolve
{
	private:
		/* 3d model pointers (Volumetric data) */
		VolumetricMesh *volumetricMesh;
		SparseMatrix *massMatrix;
		TetMesh *tetMesh;

		double mass;					// mass of the 3d model. set to 0 on init.

		/* solver parameters: */
		int numFixedVertices;
		int *constrainedVertexList;		// list of constrained vertex IDs
		double totalStrainEnergy, maxStrainEnergy;
		double *u;						// list of deformation values. u_x, u_y, u_z for each v_i.
		double *pvEnergy;				// per vertex energy calculated from the deformation values (u).

		// support vertices management
		bool loaded_s_verts;
		vector<int> num_s_verts;
		vector<int *> supportVertices;
		double * w;
		
		/* volumetric solver (slow) */
		CorotationalLinearFEM *deformableModel;

		/* The data (.cub) required for fast simulation is precalculated
		and stored in the following variables */

		// linear modes
		bool linearModesAvailable, frequenciesAvailable;
		int rLin;
		ModalMatrix * linearModalMatrix;
		double *frequencies;

		// this is an external parameter for the fast simulation. set to 20 by default.
		int numRigidModes;

		// modal derivatives
		bool modalDerivativesAvailable;
  		int numDeriv;
		ModalMatrix * modalDerivativesMatrix;

		// non linear modes
		bool nonLinearModesAvailable;
		ModalMatrix *nonLinearModalMatrix;
		int rNonLin, numComputedNonLinearModes;

		// cubic polynomials
		bool cubicPolynomialAvailable;
		StVKReducedInternalForces *cubicPolynomials;

		// rendering modal matrix and deformation mesh pointers.
		// this is required to interpolate "u" from the solver output.
		ModalMatrix * renderingModalMatrix;
		SceneObjectReducedCPU * deformableObjectRenderingMeshCPU;
		SceneObjectReduced * deformableObjectRenderingMeshReduced;

		// STVK variables required for the fast solver
		StVKReducedInternalForces *stVKReducedInternalForces;
	    StVKReducedStiffnessMatrix *stVKReducedStiffnessMatrix;
	   	
	   	ReducedStVKForceModel *reducedStVKForceModel;
	  	ReducedLinearStVKForceModel *reducedLinearStVKForceModel;
	  	ReducedForceModel *reducedForceModel;
	  	
	  	double * massMatrix_f, *u_prev;
	  	ImplicitNewmarkDense *implicitNewmarkDense;
	
		int nRendering;					// number of vertices in the rendering mesh.
		int r;							// set to 20 by default. this stores the dimension of the reduced space.
		float * URenderingFloat;		// rendering modal matrix. [currently obtained from largemodaldeformationfactory]

		// misc items required for the solver.
		int computationRunning;			// parallelize the pre-computation. (not sure why required - ref Vega)
		ObjMesh *visualMesh;			// original surface mesh. The deformation is later applied on this to visualize the output.
		int oneIndexed;					// veg file is 0 indexed or 1 indexed.
		triangle_mesh_t * m;
		
		// stores the proportion of mass that should be used for a particular vertex
		// in order to calculate the gravitational force on it.
		// the approach used for calculating the ratio is not fixed and is frequently experimented.
		// EXPERIMENTAL
		std::vector<double> ratios;

		// force vector. f_x, f_y, f_z for every vertex in the 3d model.
		double *f;

		// reduced force vector and spring force vector (reduced model)
		double *fq, *springForce;

		// misc bool variables for error checking.
		bool isInitImplicitBackwardEulerSparse;
		bool isInitImplicitNewmarkDense;

		// Applet pointer - free later.
		SimulatorApp *appPtr;


		// used to init common data - not required as of now.
		// void initCommon();
		bool readFixedVertices();
		bool calculateMass();
		// calculate the displacements for the implicit backward euler solver.
		bool calculateDisplacements(ImplicitBackwardEulerSparse *implicitBackwardEulerSparse);
		// calculate the displacements for the implicit newmark dense solver.
		bool calculateDisplacements(ImplicitNewmarkDense *implicitNewmarkDense);
		
		void applyDeformation(std::string const &path, double *u, bool is_surf) const;

		// misc helpers
		void resetVector(double *vec, int size);


		// pre-processing pipeline related functions.
		bool preProcessMesh();
		void LinearModesWorker(int &numDesiredModes, int * r, double ** frequencies_, double ** modes_);
		void ComputeModalDerivatives(int *code, double **modalDerivatives);
		void * ComputeCubicPolynomials(int * code, StVKReducedInternalForces ** newCubicPolynomial);
		void * ComputeNonLinearModes(int * code, int dataOrigin, int numNonLinearModes, double ** modes_);

	protected:
		void RemoveSixRigidModes(int numVectors, double * x);

	public:
		FEASolve(SimulatorApp * ptr);	// param constructor - ptr read to access opts.
		~FEASolve();					// destructor

		bool initImplicitBackwardEulerSparse();
		bool destroyImplicitBackwardEulerSparse();
		bool runImplicitBackwardEulerSparse();

		bool initImplicitNewmarkDense();
		bool destroyImplicitNewmarkDense();
		bool runImplicitNewmarkDense();
		void flushImplicitNewmarkDenseData();

		// getters
		int getNumVertices() const;
		// num elements is only relevant to volumetric meshes.
		int getNumElements() const;
		// returns the total absolute strain energy of the 3d model under strain.
		double getAbsoluteStrainEnergy() const;
		// returns the magnitude of maximum strain energy of all vertices.
		double getMaximumStrainEnergy() const;
		// returns the original mass (in case of volumetric mesh)
		// returns the mass estimate (in case of surface mesh)
		double getMass() const;

		double getStress() const;
		// returns the importance scores of the regions
		double * getWeightVector(int &n_) const;

		double getNormOfWeightVector() const;
		
		int getNumSamplesForOptimizer() const;

		void updateWeightVector(double * dw);
		void setWeightVector(double * w_);
		void setMass(double const &mass_);

		// calculate the per vertex energy for implicit backward euler solver
		bool calculatePerVertexEnergy();
		bool readSupportVertices();

		// save methods
		// provide filename as "<filename>.pts"
		bool saveFixedVerticesPointCloud(std::string const &path) const;
		bool saveSurfaceMesh(std::string const &path) const;
		bool saveVolumetricMesh(std::string const &path) const {;}
		void applyDeformationAndSaveVolumetricMesh(std::string const &path) const;
		void applyDeformationAndSaveSurfaceMesh(std::string const &path) const;
		bool saveModalDerivatives(std::string const &path) const;
		bool saveLinearModes(std::string const &path) const;
		bool savePerVertexEnergy(std::string const &path) const;
		bool saveMeanDeformationsPerVertex(std::string const &path) const;
		// EXPERIMENTAL: can be modified to save deformations according to need.
		bool saveCustomDeformationsPerVertex(std::string const &path) const {;}
		bool saveDeformationsPerVertex(std::string const &path) const;
};

// misc helper classes
class Timer
{
	public:
	    Timer() { clock_gettime(CLOCK_REALTIME, &beg_); }

	    double elapsed() {
	        clock_gettime(CLOCK_REALTIME, &end_);
	        return end_.tv_sec - beg_.tv_sec +
	            (end_.tv_nsec - beg_.tv_nsec) / 1000000000.;
	    }

	    void reset() { clock_gettime(CLOCK_REALTIME, &beg_); }

	private:
	    timespec beg_, end_;
};


#endif 	// FEASolve.hpp header ends here.