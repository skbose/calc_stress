#ifndef __FeaSolve__hpp__
#define __FeaSolve__hpp__

#ifndef PARDISO_SOLVER_IS_AVAILABLE
#endif
#define PARDISO_SOLVER_IS_AVAILABLE

#define NUM_COMPUTE_THREADS 1

#include "volumetricMeshLoader.h"
#include "corotationalLinearFEM.h"
#include "renderVolumetricMesh.h"
#include "StVKElementABCDLoader.h"
#include "StVKElementABCD.h"
#include "StVKInternalForces.h"
#include "corotationalLinearFEMForceModel.h"
#include "generateMassMatrix.h"
#include "implicitBackwardEulerSparse.h"
#include "loadList.h"
#include "implicitNewmarkDense.h"
#include "StVKReducedInternalForces.h"
#include "StVKReducedStiffnessMatrix.h"
#include "reducedForceModel.h"
#include "reducedLinearStVKForceModel.h"
#include "reducedStVKForceModel.h"
#include "modalMatrix.h"
#include "StVKStiffnessMatrix.h"
// #include "/home/sourav/Experimentation/VegaFEM-v3.1/utilities/largeModalDeformationFactory/largeModalDeformationFactory.h"
#include "ARPACKSolver.h"
#include "insertRows.h"
#include "matrixPCA.h"
#include "computeStiffnessMatrixNullspace.h"
#include "StVKHessianTensor.h"
#include "StVKReducedInternalForcesWX.h"
#include "matrixIO.h"
#include "sparseSolvers.h"
#include "sceneObjectReducedCPU.h"
#include "sceneObjectReducedGPU.h"
#include "sceneObjectReduced.h"
#include "objMesh.h"

// basic c++ headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
using namespace std;

class FEASolve
{
	private:
		string i_file;
		string o_file;
		string fixedVerticesFile;
		string modesFilename;
		string deformableObjectFilename;

		// 3d model params
		VolumetricMesh *volumetricMesh;
		SparseMatrix *massMatrix;
		TetMesh *tetMesh;
		double mass;
		CorotationalLinearFEM *deformableModel;
		int numFixedVertices;
		int *constrainedVertexList;
		int numTimeSteps;
		double totalStrainEnergy, maxStrainEnergy;
		double *u;
		double *pvEnergy;

		// linear modes info
		bool linearModesAvailable, frequenciesAvailable;
		int rLin;
		ModalMatrix * linearModalMatrix;
		double *frequencies;

		int numRigidModes;

		// modal derivatives info
		bool modalDerivativesAvailable;
  		int numDeriv;
		ModalMatrix * modalDerivativesMatrix;

		// non linear modes info
		bool nonLinearModesAvailable;
		ModalMatrix *nonLinearModalMatrix;
		int rNonLin, numComputedNonLinearModes;

		// cubic polynomials info
		bool cubicPolynomialAvailable;
		StVKReducedInternalForces *cubicPolynomials;

		// misc
		int computationRunning;
		ObjMesh *visualMesh;
		std::vector<double> ratios;

		// fast reduced deformation variables
		ModalMatrix * renderingModalMatrix;
		SceneObjectReducedCPU * deformableObjectRenderingMeshCPU;
		SceneObjectReduced * deformableObjectRenderingMeshReduced;
		int nRendering, r;
		float * URenderingFloat;

	protected:
		void RemoveSixRigidModes(int numVectors, double * x);

	public:
		// constructors
		FEASolve();
		FEASolve(string i_file_, string o_file_, string fixedVerticesFile_) : i_file(i_file_), o_file(o_file_), fixedVerticesFile(fixedVerticesFile_) {};

		// destructor
		~FEASolve();
		
		// member functions
		void init(std::string modalMatrixFilePath, std::string triangleMeshFilePath, std::string massRatioFilePath);
		void initSolver(int steps);
		void initSolverNewmarkSparse(int steps);
		void initSolverNewmarkDense(string cubicPolyFilename, int steps);
		void simulateImplicitNewmarkDense(int steps, ImplicitNewmarkDense *implicitNewmarkDense);

		// getters
		int getNumVertices();
		int getNumElements();
		double getStrainEnergy();
		double getAbsStrainEnergy();
		double getMaxStrainEnergy() {return maxStrainEnergy;}
		double getMass() { return mass; }

		// read vertices
		int * readFixedVertices(string filename, int &numConstrainedDOFS, int oneIndexed);
		bool saveFixedVerticesPointCloud(std::string const &path);

		// setters
		void setGravity();

		// misc functions
		void calculateMass();
		void calculateDisplacements(ImplicitBackwardEulerSparse *implicitBackwardEulerSparse);
		void calculateDisplacements(ImplicitNewmarkDense *implicitNewmarkDense);
		void calculatePerVertexEnergy();
		void saveDeformedMesh(string filename);
		void resetVector(double *vec, int num)
		{
			memset(vec, 0, num);
		}
		void preprocessMesh();
		void LinearModesWorker(int &numDesiredModes, int * r, double ** frequencies_, double ** modes_);
		void ComputeModalDerivatives(int *code, double **modalDerivatives);
		void * ComputeCubicPolynomials(int * code, StVKReducedInternalForces ** newCubicPolynomial);
		void * ComputeNonLinearModes(int * code, int dataOrigin, int numNonLinearModes, double ** modes_ );


		// save methods
		void SaveModalDerivatives(string path);
		void SaveLinearModes(string path);
		void savePerVertexEnergy(string path)
		{
			ofstream output(path.c_str());

			for (int i = 0; i < getNumVertices(); i++)
				output << pvEnergy[i] << endl;
			output.close();

			cout << "Saved per-vertex energy data to file!" << endl;
		}

		void saveMeanDeformations(string path)
		{
			ofstream output(path.c_str());

			for (int i = 0; i < nRendering; i++)
			{
				double mean = fabs(u[i]) + fabs(u[i+1]) + fabs(u[i+2]);
				mean = mean / 3.0;

				output << mean << endl;
			}

			output.close();
		}

		void saveMaxDeformations(string path)
		{

			ofstream output(path.c_str());

			for (int i = 0; i < nRendering; i++)
			{
				double max = fmax(fabs(u[i]), fabs(u[i+1]));
				max = fmax(max, fabs(u[i+2]));

				output << max << endl;
			}

			output.close();	
		}

		void saveDeformations(string path)
		{
			ofstream output(path.c_str());

			for (int i = 0; i < nRendering; i++)
				output << u[i] << " " << u[i + 1] << " " << u[i + 2] << endl;
			output.close();
		}

		// display
		void showPerVertexEnergy()
		{
			for (int i = 0; i < getNumVertices(); i++)
				cout << pvEnergy[i] << endl;
		}

		void verifyOrder()
		{
			for (int i = 0; i < nRendering; i++)
			{
				Vec3d pos_ = visualMesh->getPosition(i);
				int index = deformableObjectRenderingMeshCPU->GetClosestVertex(pos_);

				if (i != index)
					cout << "May Day!" << endl;
			}
		}
};

#endif