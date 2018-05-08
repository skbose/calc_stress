// personal headers
#include "FeaSolve.hpp"

int FEASolve::getNumVertices()
{
	return volumetricMesh->getNumVertices();
}

int FEASolve::getNumElements()
{
	return volumetricMesh->getNumElements();
}

void FEASolve::init(std::string modalMatrixFilePath, std::string triangleMeshFilePath, std::string massRatioFilePath)
{
	volumetricMesh = VolumetricMeshLoader::load(i_file.c_str());
	
	pvEnergy = NULL;
	u = NULL;
	massMatrix = NULL;
	tetMesh = NULL;
	constrainedVertexList = NULL;
	
	// linear modes info
	linearModesAvailable = frequenciesAvailable = false;
	linearModalMatrix = NULL;
	rLin = -1;
	frequencies = NULL;

	// modal derivates info
	modalDerivativesAvailable = false;
	modalDerivativesMatrix = NULL;
	numDeriv = -1;

	// non-linear modes info
	nonLinearModalMatrix = NULL;
	nonLinearModesAvailable = false;
	rNonLin = -1, numComputedNonLinearModes = -1;

	// cubic polynomial info
	cubicPolynomialAvailable = false;
	cubicPolynomials = NULL;

	// fast reduced deformation varaibles init to NULL
	renderingModalMatrix = NULL;
	deformableObjectRenderingMeshCPU = NULL;
	deformableObjectRenderingMeshReduced = NULL;


	if (volumetricMesh == NULL)
	{
		printf("Error: failed to load mesh.\n");
		return;
	}
	else
		printf("Success. Number of vertices: %d . Number of elements: %d .\n", volumetricMesh->getNumVertices(), volumetricMesh->getNumElements());

	calculateMass();
	
	if (volumetricMesh->getElementType() == VolumetricMesh::TET)
		tetMesh = (TetMesh*) volumetricMesh; // such down-casting is safe in Vega
	else
	{
		printf("Error: not a tet mesh.\n");
		exit(1);
	}

	deformableModel = new CorotationalLinearFEM(tetMesh);

	// modesFilename.assign("/home/sourav/Experimentation/VegaFEM-v3.1/models/sims/male_sitting_stool/2/male_sitting_stool.U");
	// deformableObjectFilename.assign("/home/sourav/Experimentation/VegaFEM-v3.1/models/sims/male_sitting_stool/mesh.obj");

	modesFilename = modalMatrixFilePath;
	deformableObjectFilename = triangleMeshFilePath;

	ReadMatrixFromDisk_(modesFilename.c_str(),&nRendering,&r,&URenderingFloat);
	nRendering /= 3;
	cout << "Number of vertices (read from Modal Matrix): " << nRendering << ", r: " << r << endl;

	double * URendering = (double*) malloc (sizeof(double) * 3 * nRendering * r);
	
	for(int i=0; i < 3 * nRendering * r; i++)
		URendering[i] = URenderingFloat[i];
	free(URenderingFloat);
	
	renderingModalMatrix = new ModalMatrix(nRendering,r,URendering);
	free(URendering); // ModalMatrix made an internal copy
	
	deformableObjectRenderingMeshCPU = new SceneObjectReducedCPU(deformableObjectFilename.c_str(), renderingModalMatrix);
	deformableObjectRenderingMeshReduced = deformableObjectRenderingMeshCPU;

	// load triangle mesh
	visualMesh = new ObjMesh(deformableObjectFilename);

	// load mass ratio
	ifstream f(massRatioFilePath.c_str());
	std:string line;

	double ratio = 0.0;
	// Mass Ratio file - uncomment later
	// while (std::getline(f, line))
	// {
	// 	std::istringstream linestream(line);
	// 	linestream >> ratio;
	// 	ratios.push_back(ratio);
	// }
	double *vertexVolumes;
	vertexVolumes = (double *) malloc(sizeof(double) * getNumVertices());
	volumetricMesh->getVertexVolumes(vertexVolumes);
	double totalVolume = 0.0;
	
	for (int i = 0; i < getNumVertices(); i++)
	{
		totalVolume += vertexVolumes[i];
	}

	for (int i = 0; i < getNumVertices(); i++)
	{
		ratios.push_back(vertexVolumes[i] / totalVolume);
	}
	free(vertexVolumes); vertexVolumes = NULL;

	// verifyOrder();
}

void FEASolve::initSolver(int steps)
{
	ForceModel * forceModel = new CorotationalLinearFEMForceModel(deformableModel);
	int r = 3 * getNumVertices(); // total number of DOFs
	double timestep = 0.0333;

	GenerateMassMatrix::computeMassMatrix(tetMesh, &massMatrix, true);

	double dampingMassCoef = 0.0; // "underwater"-like damping (here turned off)
	double dampingStiffnessCoef = 0.0; // (primarily) high-frequency damping
	
	int numConstrainedDOFs = 0;
	int *constrainedDOFs = readFixedVertices(fixedVerticesFile, numConstrainedDOFs, 0);

	// initialize the integrator
	ImplicitBackwardEulerSparse * implicitBackwardEulerSparse = new
	ImplicitBackwardEulerSparse(r, timestep, massMatrix, forceModel, numConstrainedDOFs, constrainedDOFs, dampingMassCoef, dampingStiffnessCoef);

	// alocate buffer for external forces
	double * f = (double*) malloc (sizeof(double) * r);
	implicitBackwardEulerSparse->SetExternalForcesToZero();
	volumetricMesh->computeGravity(f, 9.81, true);

	// mapping forces to the reduced domain


	for (int i = 0; i < steps; i++)
	{
		cout << "Simulating step " << i << endl;
		if (i == 0) // set some force at the first timestep
		{
			implicitBackwardEulerSparse->SetExternalForces(f);
		}
		implicitBackwardEulerSparse->DoTimestep();
	}

	calculateDisplacements(implicitBackwardEulerSparse);
	// calculatePerVertexEnergy();

	// free up force vector
	free(f);
}

void FEASolve::calculateDisplacements(ImplicitBackwardEulerSparse *implicitBackwardEulerSparse)
{
	int r = getNumVertices() * 3;
	u = (double*) malloc (sizeof(double) * r);
	implicitBackwardEulerSparse->GetqState(u);
}

void FEASolve::calculateDisplacements(ImplicitNewmarkDense *implicitNewmarkDense)
{
	// workplace
	// int r = 20;
	// u = (double*) malloc (sizeof(double) * r);
	// implicitNewmarkDense->GetqState(u);


	double * q = (double *) malloc (sizeof(double) * r);

	memcpy(q, implicitNewmarkDense->Getq(), sizeof(double) * r);

	deformableObjectRenderingMeshReduced->Setq(q);
  	deformableObjectRenderingMeshReduced->Compute_uUq();

  	if (u == NULL)
  		u = (double *) malloc(sizeof(double) * nRendering * 3);
  	deformableObjectRenderingMeshReduced->Getu(u);

  	// for (int i = 0; i < nRendering * 3; i++)
  	// 	cout << u[i] << endl;
	free(q);
}

double FEASolve::getAbsStrainEnergy()
{
	return totalStrainEnergy;
}

void FEASolve::calculatePerVertexEnergy()
{
	StVKElementABCD * stvk_element = StVKElementABCDLoader::load(volumetricMesh);
	StVKInternalForces * stvk = new StVKInternalForces(volumetricMesh, stvk_element, true);

	pvEnergy = (double*)malloc(sizeof(double) * tetMesh->getNumVertices());

	int elementLow = 0;
	int elementHigh = getNumElements();
	int numVertices = getNumVertices();

	double *buffer = (double *)malloc(sizeof(double) * (getNumVertices() * 3));
	double *energy_ = (double *)malloc(sizeof(double) * (getNumVertices() * 3));

	stvk->AddLinearTermsContribution(u, buffer, elementLow, elementHigh);
	for (int i = 0; i < 3 * numVertices; i++)
    	energy_[i] += 0.5 * buffer[i] * u[i];

    resetVector(buffer, numVertices * 3);
	stvk->AddQuadraticTermsContribution(u, buffer, elementLow, elementHigh);
	double oneThird = 1.0 / 3;
	for (int i = 0; i < 3 * numVertices; i++)
    	energy_[i] += oneThird * buffer[i] * u[i];

    resetVector(buffer, numVertices * 3);
	stvk->AddCubicTermsContribution(u, buffer, elementLow, elementHigh);
	double oneQuarter = 1.0 / 4;
	for (int i = 0; i < 3 * numVertices; i++)
    	energy_[i] += oneQuarter * buffer[i] * u[i];

    totalStrainEnergy = 0;
    double maxEnergy = -1;
    for (int i = 0; i < numVertices; i += 1)
    {
    	double abs_energy = fabs(energy_[i*3]) + fabs(energy_[i*3 + 1]) + fabs(energy_[i*3 + 2]);
    	pvEnergy[i] = abs_energy;

    	if (abs_energy > maxEnergy)
    		maxEnergy = abs_energy;
    	totalStrainEnergy += abs_energy;
    }
    maxStrainEnergy = maxEnergy;

    // free up buffers
    free(buffer); buffer = NULL;
    free(energy_); energy_ = NULL;
}

void FEASolve::calculateMass()
{
	mass = 0;
	int numElements = volumetricMesh->getNumElements();

	for(int el = 0; el < numElements; el++)
  	{
	    double volume = volumetricMesh->getElementVolume(el);
	    double density = volumetricMesh->getElementDensity(el);
	    double e_mass = density * volume;
	    mass += e_mass;
	}
}

int * FEASolve::readFixedVertices(string filename, int &numConstrainedDOFs, int oneIndexed)
{
	int numConstrainedVertices;
	constrainedVertexList = NULL;

	LoadList::load(filename.c_str(), &numConstrainedVertices, &constrainedVertexList);

	numConstrainedDOFs = numConstrainedVertices * 3;
	int *constrainedDOFs = (int *)malloc(sizeof(int) * numConstrainedDOFs);

	std::cout << "Number of constrained vertices read: " << numConstrainedVertices << std::endl;

	for (int i = 0, pos = 0; i < numConstrainedVertices; i++, pos++)
	{
		// cout << i << endl;
		// to make the vertex list start from 0 index
		constrainedDOFs[pos*3 + 0] = (constrainedVertexList[i] - 1 + oneIndexed) * 3 + oneIndexed;
		constrainedDOFs[pos*3 + 1] = (constrainedVertexList[i] - 1 + oneIndexed) * 3 + 1 + oneIndexed;
		constrainedDOFs[pos*3 + 2] = (constrainedVertexList[i] - 1 + oneIndexed) * 3 + 2 + oneIndexed;
	}

	// cout << "Exits loop" << endl;
	// update the global variable
	numFixedVertices = numConstrainedVertices;

	return constrainedDOFs;
}

void FEASolve::saveDeformedMesh(string filename)
{
	volumetricMesh->applyDeformation(u);
	volumetricMesh->save(filename.c_str());
}

FEASolve::~FEASolve()
{
	if (pvEnergy)
	{
		free(pvEnergy);
		pvEnergy = NULL;
	}
	if (u)
	{
		free(u);
		u = NULL;
	}
	if (constrainedVertexList)
	{	
		free(constrainedVertexList);
		constrainedVertexList = NULL;
	}	
}

bool FEASolve::saveFixedVerticesPointCloud(std::string const &path)
{
	std::ofstream f(path.c_str());
  
	if (!f)
	{
		std::cerr << "Couldn't open file " << path << std::endl;
		return false;
	}

	for (int i = 0; i < numFixedVertices; i++)
	{		
		int indexOfConstrainedVertex = constrainedVertexList[i] - 1;
		if (indexOfConstrainedVertex < 4128)
			continue;
		Vec3d pos(visualMesh->getPosition(indexOfConstrainedVertex));

		f << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	}

	return true;

}
void FEASolve::initSolverNewmarkSparse(int steps)
{
	ForceModel * forceModel = new CorotationalLinearFEMForceModel(deformableModel);
	int r = 3 * getNumVertices(); // total number of DOFs
	double timestep = 0.0333;

	GenerateMassMatrix::computeMassMatrix(tetMesh, &massMatrix, true);

	double dampingMassCoef = 0.0; // "underwater"-like damping (here turned off)
	double dampingStiffnessCoef = 0.0; // (primarily) high-frequency damping
	
	int numConstrainedDOFs = 0;
	int *constrainedDOFs = readFixedVertices(fixedVerticesFile, numConstrainedDOFs, 0);

	// initialize the integrator
	ImplicitNewmarkSparse * implicitNewmarkSparse = new ImplicitNewmarkSparse(r, timestep, massMatrix, forceModel, numConstrainedDOFs, constrainedDOFs, dampingMassCoef, dampingStiffnessCoef);

	// alocate buffer for external forces
	double * f = (double*) malloc (sizeof(double) * r);
	implicitNewmarkSparse->SetExternalForcesToZero();
	volumetricMesh->computeGravity(f, 9.81, true);

	for (int i = 0; i < steps; i++)
	{
		cout << "Simulating step " << i << endl;
		if (i == 0) // set some force at the first timestep
		{
			implicitNewmarkSparse->SetExternalForces(f);
		}
		implicitNewmarkSparse->DoTimestep();
	}

	// calculateDisplacements(implicitNewmarkSparse);
	// calculatePerVertexEnergy();

	// free up force vector
	free(f); f = NULL;
}

void FEASolve::initSolverNewmarkDense(string cubicPolyFilename, int steps)
{
	// workplace
	// r = 20;
	double * massMatrix = (double*) calloc(r*r,sizeof(double));
	double timestep = 0.0003;

	int numConstrainedDOFs = 0;
	int *constrainedDOFs = readFixedVertices(fixedVerticesFile, numConstrainedDOFs, 0);

  	for(int i = 0; i < r; i++)
  	{	
    	massMatrix[ELT(r,i,i)] = 1.0;
  	}

    StVKReducedInternalForces *stVKReducedInternalForces = new StVKReducedInternalForces(cubicPolyFilename.c_str());
    StVKReducedStiffnessMatrix *stVKReducedStiffnessMatrix = new StVKReducedStiffnessMatrix(stVKReducedInternalForces);

   	ReducedStVKForceModel *reducedStVKForceModel = new ReducedStVKForceModel(stVKReducedInternalForces, stVKReducedStiffnessMatrix);
  	ReducedLinearStVKForceModel *reducedLinearStVKForceModel = new ReducedLinearStVKForceModel(stVKReducedStiffnessMatrix);
  	ReducedForceModel *reducedForceModel = reducedStVKForceModel;

  	ImplicitNewmarkDense *implicitNewmarkDense = new ImplicitNewmarkDense(r, timestep, massMatrix, reducedForceModel);

 	
	double * f = (double*) malloc (sizeof(double) * 3 * nRendering);
	// reduced force vector
	double * fq = (double *) malloc (sizeof(double) * r);
	double * springForce = (double *) malloc(sizeof(double) * 3 * nRendering);

	implicitNewmarkDense->SetExternalForcesToZero();
	
	// flush the force vector to zero
	for (int i = 0; i < 3 * nRendering; i++)
	{
		f[i] = 0; springForce[i] = 0;
	}

	// volumetricMesh->computeGravity(f, 9.81 * 10, true);
	
	// TODO - recheck gravity calculations.
	double gravity_force_per_vertex = (mass * 9.81 * 10) / (double)nRendering;
	for (int i = 0; i < nRendering; i++)
		f[i + 1] = -mass * ratios[i] * 9.81 * 5;

	// init fq to zero
	for (int i = 0; i < r; i++)
		fq[i] = 0;

	// TODO - remove this later, for debugging purposes.
	// std::string fqFilePath = "test/fq.data";
	// std::string line;

	// ifstream fqFileObj(fqFilePath.c_str());
	// getline(fqFileObj, line);

	// stringstream ss(line);

	// double force;
	// int ind = 0;
	// while (ss >> force)
	// {
	// 	fq[ind++] = force;
	// }

	renderingModalMatrix->ProjectSingleVertex(0, f[0], f[1], f[2], fq);
	
	for (int pulledVertex = 1; pulledVertex < nRendering; pulledVertex++)
		renderingModalMatrix->AddProjectSingleVertex(pulledVertex,
           	f[pulledVertex], f[pulledVertex + 1], f[pulledVertex + 2], fq);

	float newmarkBeta = 0.25;
	float newmarkGamma = 0.5;
	float timeStep = 0.0003;

	implicitNewmarkDense->SetTimestep(timeStep);
	implicitNewmarkDense->SetNewmarkBeta(newmarkBeta);
	implicitNewmarkDense->SetNewmarkGamma(newmarkGamma);

	for (int i = 0; i < steps; i++)
	{
		// cout << "Simulating step " << i << endl;
		if (i == 0) // set some force at the first timestep
		{
			implicitNewmarkDense->SetExternalForces(fq);
		}
		implicitNewmarkDense->DoTimestep();

		// flush existing forces to zero
		// for (int k = 0; k < r; k++)
		// 	fq[k] = 0;

		// // update the forces acting on the fixed vertices everytime.
		// calculateDisplacements(implicitNewmarkDense);

		// for (int vertex_i = 0; vertex_i < numFixedVertices; vertex_i++)
		// {
		// 	int v = constrainedVertexList[vertex_i];
		// 	springForce[v + 0] = -2 * u[v + 0];
		// 	springForce[v + 1] = -2 * u[v + 1];
		// 	springForce[v + 2] = -2 * u[v + 2];
 	//  	}

 	//  	// TODO: clean this up. Adding up spring force and gravity force
 	//  	for (int vertex_i = 0; vertex_i < getNumVertices(); vertex_i++)
 	//  	{	
 	//  		springForce[vertex_i + 0] += f[vertex_i + 0];
 	//  		springForce[vertex_i + 1] += f[vertex_i + 1];
 	//  		springForce[vertex_i + 2] += f[vertex_i + 2];
 	//  	}

 	//  	renderingModalMatrix->ProjectSingleVertex(0, springForce[0], springForce[1], springForce[2], fq);
	
		// for (int pulledVertex = 1; pulledVertex < nRendering; pulledVertex++)
		// 	renderingModalMatrix->AddProjectSingleVertex(pulledVertex,
  //          		springForce[pulledVertex], springForce[pulledVertex + 1], springForce[pulledVertex + 2], fq);
 	 	
 	//  	implicitNewmarkDense->SetExternalForcesToZero();
 	//  	implicitNewmarkDense->SetExternalForces(fq);
	}

	free(f); f = NULL;
	free(fq); fq = NULL;

	calculateDisplacements(implicitNewmarkDense);
	// calculatePerVertexEnergy();

	// TODO - add a method later.
  	visualMesh->deform(u);
  	visualMesh->save("deformedMesh.obj");
}

void FEASolve::simulateImplicitNewmarkDense(int steps, ImplicitNewmarkDense *implicitNewmarkDense)
{
	;
}

void FEASolve::preprocessMesh()
{
	int newr;
	int numLinearModes = 16;
	double *newFrequencies = NULL;
	double *newLinearModes = NULL;

	LinearModesWorker(numLinearModes, &newr, &newFrequencies, &newLinearModes);
	
	if (newr != numLinearModes)
	{
		cout << "Linear Modes calculation failed!" << endl;
	}
	else
	{
		cout << "Linear Modes successfully calculated!" << endl;
		linearModesAvailable = true;
		delete(linearModalMatrix);
		rLin = newr;

		linearModalMatrix = new ModalMatrix(volumetricMesh->getNumVertices(), rLin, newLinearModes);
      	free(newLinearModes);
      	free(frequencies);
      	frequencies = newFrequencies;
      	frequenciesAvailable = true;


      	// set the rigid modes based on the number of fixed vertices
		if (numFixedVertices >= 3)
		{
			numRigidModes = 0;
		}
		else if (numFixedVertices == 2)
		{
			numRigidModes = 1;
		}
		else if (numFixedVertices == 1)
		{
			numRigidModes = 3;
		}
		else
		{
			numRigidModes = 6;
		}
	}

	if (!linearModesAvailable)
  	{
    	cout << "Cannot compute modal derivatives. Linear modes not available!" << endl;
    	return;
  	}

  	// linear modes are available here. Start modal derivatives calculation.
	double * modalDerivatives = NULL;
	int code;
	ComputeModalDerivatives(&code, &modalDerivatives);
	
	if (code != 0)
	{
		cout << "Modal derivatives calculation failed!" << endl;
		free(modalDerivatives);
		return;
	}
	else
	{
		cout << "Modal derivative computation succeeded." << endl;
		modalDerivativesAvailable = true;

		delete(modalDerivativesMatrix);
		modalDerivativesMatrix = new ModalMatrix(
	    volumetricMesh->getNumVertices(), 
	    (rLin - numRigidModes) * (rLin - numRigidModes + 1) / 2, 
	    modalDerivatives);
		free(modalDerivatives);
	}

	// linear modes and modal derivatives calculation done. Starting non-linear modes calculates.
	if (!(linearModesAvailable && modalDerivativesAvailable && frequenciesAvailable))
	{
		cout << "Non-Linear modes calculation failed!" << endl;
		return;
	}

	double * newNonLinearModes = NULL;
	int dataOrigin = 0;

	// set the number of non linear modes here
	numComputedNonLinearModes = 20;
	ComputeNonLinearModes(&code, dataOrigin, numComputedNonLinearModes, &newNonLinearModes);

	if (code != 0)
    {
      char s[96];
      if (code < 0)
        sprintf(s, "Nonlinear mode computation failed. Memory allocation error.");
      else
        sprintf(s, "Nonlinear mode computation failed. dgesvd exit code: %d\n", code);
      cout << s << endl;
      free(newNonLinearModes);
      return;
    }
    else
    {
      nonLinearModesAvailable = true;
 
      delete(nonLinearModalMatrix);
      rNonLin = numComputedNonLinearModes;
      nonLinearModalMatrix = new ModalMatrix(
        volumetricMesh->getNumVertices(), 
        rNonLin, newNonLinearModes);
      free(newNonLinearModes);

      cout << "\nNon-Linear modes successfully calculated!" << endl;
    }

	// linear modes, modal derivatives and non-linear modes calculation done. Starting cubic polynomial calculation.
    StVKReducedInternalForces * newCubicPolynomials = NULL;
    ComputeCubicPolynomials(&code, &newCubicPolynomials);

    if (code != 0)
    {
    	cout << "\nCubic Polynomial calculation failed!" << endl;
		free(newCubicPolynomials);
		return;
    }
    else
    {
		delete(cubicPolynomials);
		cubicPolynomials = newCubicPolynomials;
		cubicPolynomialAvailable = true;

		cout << "\nCubic Polynomials successfully calculated!" << endl;
    }
}

void FEASolve::LinearModesWorker(int &numDesiredModes, int * r, double ** frequencies_, double ** modes_)
{
	*r = -1;

	// create mass matrix
	SparseMatrix * massMatrix;
	GenerateMassMatrix::computeMassMatrix(volumetricMesh, &massMatrix, true);

	// create stiffness matrix
	StVKElementABCD * precomputedIntegrals = StVKElementABCDLoader::load(volumetricMesh);
	StVKInternalForces * internalForces = 
		new StVKInternalForces(volumetricMesh, precomputedIntegrals);

	SparseMatrix * stiffnessMatrix;
	StVKStiffnessMatrix * stiffnessMatrixClass = new StVKStiffnessMatrix(internalForces);
	stiffnessMatrixClass->GetStiffnessMatrixTopology(&stiffnessMatrix);
	double * zero = (double*) calloc(3 * volumetricMesh->getNumVertices(), sizeof(double));
	stiffnessMatrixClass->ComputeStiffnessMatrix(zero, stiffnessMatrix);

	free(zero);
	delete(precomputedIntegrals);
	delete(stiffnessMatrixClass);
	delete(internalForces);

	// constrain the degrees of freedom
	// int * constrainedDOFs = (int*) malloc (sizeof(int) * 3 * numConstrainedVertices);
	// set<int> :: iterator iter;
	// int i = 0;
	// for(iter = precomputationState.fixedVertices.begin(); iter != precomputationState.fixedVertices.end(); iter++)
	// {
	// 	constrainedDOFs[3*i+0] = 3 * (*iter) + 1;
	// 	constrainedDOFs[3*i+1] = 3 * (*iter) + 2;
	// 	constrainedDOFs[3*i+2] = 3 * (*iter) + 3;
	// 	i++;
	// }
	
	// read and update constrained vertices
	int numConstrainedDOFs;
	int oneIndexed = 1;
	int *constrainedDOFs = readFixedVertices(fixedVerticesFile, numConstrainedDOFs, oneIndexed);
	int numConstrainedVertices = numFixedVertices;

	massMatrix->RemoveRowsColumns(
	3 * numConstrainedVertices, constrainedDOFs, oneIndexed);

	stiffnessMatrix->RemoveRowsColumns(
	3 * numConstrainedVertices, constrainedDOFs, oneIndexed);

	// call ARPACK

	double * frequenciesTemp = (double*) malloc (sizeof(double) * numDesiredModes);
	int numRetainedDOFs = stiffnessMatrix->Getn();
	double * modesTemp = (double*) malloc (sizeof(double) * numDesiredModes * numRetainedDOFs);

	printf("Computing linear modes using ARPACK: ...\n");
	// PerformanceCounter ARPACKCounter;
	double sigma = -1.0;

	int numLinearSolverThreads = 4;
	if (numLinearSolverThreads > 3)
		numLinearSolverThreads = 3; // diminished returns in solver beyond 3 threads

	//massMatrix->Save("MFactory");
	//stiffnessMatrix->Save("KFactory");

	ARPACKSolver generalizedEigenvalueProblem;
	int nconv = generalizedEigenvalueProblem.SolveGenEigShInv
	(stiffnessMatrix, massMatrix, 
	 numDesiredModes, frequenciesTemp, 
	 modesTemp, sigma, numLinearSolverThreads);

	// ARPACKCounter.StopCounter();
	// double ARPACKTime = ARPACKCounter.GetElapsedTime();
	// printf("ARPACK time: %G s.\n", ARPACKTime); fflush(NULL);

	if (nconv < numDesiredModes)
	{
		free(modesTemp);
		free(frequenciesTemp);
		*r = -3;
		free(constrainedDOFs);
		delete(massMatrix);
		delete(stiffnessMatrix);
		return;
	}

	int n3 = 3 * volumetricMesh->getNumVertices();
	*frequencies_ = (double*) calloc (numDesiredModes, sizeof(double));
	*modes_ = (double*) calloc (numDesiredModes * n3, sizeof(double));

	for(int i=0; i<numDesiredModes; i++)
	{
		// insert zero rows into the computed modes
		int oneIndexed = 1;
		InsertRows(n3, &modesTemp[numRetainedDOFs*i], &((*modes_)[n3*i]), 
		  3 * numConstrainedVertices, constrainedDOFs, oneIndexed);
	}

	for(int i=0; i<numDesiredModes; i++)
	{
		if (frequenciesTemp[i] <= 0)
		  (*frequencies_)[i] = 0.0;
		else
		  (*frequencies_)[i] = sqrt((frequenciesTemp)[i]) / (2 * M_PI);
	}

	free(modesTemp);
	free(frequenciesTemp);
	free(constrainedDOFs);

	delete(massMatrix);
	delete(stiffnessMatrix);

	*r = numDesiredModes;
	linearModesAvailable = true;

	return;
}

void FEASolve::ComputeModalDerivatives(int *code, double **modalDerivatives)
{
	*code = 0;

	// create stiffness matrix
	StVKElementABCD * precomputedIntegrals = StVKElementABCDLoader::load(volumetricMesh);
	StVKInternalForces * internalForces = 
	new StVKInternalForces(volumetricMesh, precomputedIntegrals); 

	// create stiffness matrix
	int n3 = 3 * volumetricMesh->getNumVertices();
	SparseMatrix * stiffnessMatrix;
	StVKStiffnessMatrix * stiffnessMatrixClass = new StVKStiffnessMatrix(internalForces);
	stiffnessMatrixClass->GetStiffnessMatrixTopology(&stiffnessMatrix);
	double * zero = (double*) calloc(n3, sizeof(double));
	stiffnessMatrixClass->ComputeStiffnessMatrix(zero, stiffnessMatrix);
	free(zero);

	// now, the stiffness matrix is computed
	// constrain the degrees of freedom
	int numConstrainedVertices = (int) (numFixedVertices);
	int * constrainedDOFs = (int*) malloc (sizeof(int) * 3 * numConstrainedVertices);
	for (int i = 0, pos = 0; i < numFixedVertices; i++, pos++)
	{
		int * iter = constrainedVertexList + i;
		constrainedDOFs[3*pos+0] = 3 * (*iter) + 1;
		constrainedDOFs[3*pos+1] = 3 * (*iter) + 2;
		constrainedDOFs[3*pos+2] = 3 * (*iter) + 3;
	}

	int oneIndexed = 1;
	stiffnessMatrix->RemoveRowsColumns(
	3 * numConstrainedVertices, constrainedDOFs, oneIndexed);

	int numRetainedDOFs = stiffnessMatrix->Getn();

	// generate rhs side
	bool computeHessianAtZero = (volumetricMesh->getNumElements() < 5000);

	if (computeHessianAtZero)
		printf("Hessian at zero will be computed explicitly.\n");
	else
		printf("Hessian at zero will not be computed explicitly.\n");

	StVKHessianTensor * stVKStiffnessHessian = new 
	StVKHessianTensor(stiffnessMatrixClass);

	int numUsedLinearModes = rLin - numRigidModes;
	numDeriv = numUsedLinearModes * (numUsedLinearModes + 1) / 2;
	double * rhs = (double*) malloc (sizeof(double) * n3 * numDeriv);
	if (!rhs)
	{
		printf("Error: could not allocate space for all modal derivatives.\n");
		*code = 1;
		delete(precomputedIntegrals);
		delete(stiffnessMatrixClass);
		delete(internalForces);
		return;
	}

	if (computeHessianAtZero)
	{
		// compute hessian at zero
		if (stVKStiffnessHessian->ComputeHessianAtZero() != 0)
		{
			printf("Error: failed to evaluate the Hessian at the origin.\n");
			*code = 1;
			delete(precomputedIntegrals);
			delete(stiffnessMatrixClass);
			delete(internalForces);
			return;
		}
	}

	printf("Preparing to compute %d modal derivatives...\n", numDeriv);

	double * Ulin = linearModalMatrix->GetMatrix();
	if (computeHessianAtZero)
	{
		printf("Using the high-memory version.\n");
		int pos = 0;
		for(int i=0; i<numUsedLinearModes; i++)
		{
			printf("%d: ",i);fflush(NULL);
			for(int j=i; j<numUsedLinearModes; j++)
			{
				printf("%d ",j);fflush(NULL);
		    	stVKStiffnessHessian->EvaluateHessianQuadraticForm(
		      	&Ulin[n3*(numRigidModes + i)], &Ulin[n3*(numRigidModes + j)], &rhs[ELT(n3,0,pos)]);
		    
		    	for(int k=0; k<n3; k++) //multiply by -1
		      		rhs[ELT(n3,k,pos)] *= -1.0;

		    	pos++;
		  	}
		  	printf("\n");
		}
	}
	else
	{
		printf("Using the low-memory version.\n");
		stVKStiffnessHessian->EvaluateHessianQuadraticFormDirectAll(
		  Ulin,rLin,rhs,numRigidModes);

		if (n3 * numDeriv < 0)
		{
			printf("Error: data too large to be indexed with the word size of your machine.\n");
			*code = 2;
			delete(stVKStiffnessHessian);
			delete(stiffnessMatrix);
			free(rhs);
			delete(precomputedIntegrals);
			delete(stiffnessMatrixClass);
			delete(internalForces);
			return;
		}

	/*
	if ((n3 > 200000) || (precomputationState.numDeriv > 1000))
	{
	  printf("Warning: size of data %d might be too large to be indexed with the word size of your machine.\n",n3*precomputationState.numDeriv);
	}
	*/

	//multiply by -1
	for(int i=0; i<n3*numDeriv; i++)
	  rhs[i] *= -1.0;
	}

	printf("Right-hand sides for modal derivatives computed.\n");fflush(NULL);

	delete(stVKStiffnessHessian);
	delete(precomputedIntegrals);
	delete(stiffnessMatrixClass);
	delete(internalForces);

	// create mass matrix
	SparseMatrix * massMatrix;
	GenerateMassMatrix::computeMassMatrix(volumetricMesh, &massMatrix, true);

	if (numRigidModes < 6)
	{
		double * buffer0 = (double*) malloc (sizeof(double) * n3);
		for(int i=0; i<numDeriv; i++)
		{
			massMatrix->MultiplyVector(&rhs[ELT(n3,0,i)], buffer0);
			for(int j=0; j<numRigidModes; j++)
			{
				// rhs -= <rhs, rigid mode j> * rigid mode j
				// rigid modes are mass-orthonormal
				double dotp = 0.0;
				for(int k=0; k<n3; k++)
				  dotp += buffer0[k] * Ulin[ELT(n3,k,j)]; 
				for(int k=0; k<n3; k++)
				  rhs[ELT(n3,k,i)] -= dotp * Ulin[ELT(n3,k,j)];
			}
		}
		free(buffer0);
	}
	else
	{
		RemoveSixRigidModes(numDeriv, rhs);
	}

	// constrain rhs
	double * rhsConstrained = (double*) malloc (sizeof(double) * numRetainedDOFs * numDeriv);
	for(int i=0; i<numDeriv; i++)
	RemoveRows(n3, &rhsConstrained[numRetainedDOFs*i], 
	  &rhs[n3*i], 3 * numConstrainedVertices, constrainedDOFs, oneIndexed);

	free(rhs);

	// make room for (uninflated) derivatives
	double * modalDerivativesConstrained = (double*) malloc (sizeof(double) * numDeriv * numRetainedDOFs);

	// solve K * modesTemp = rhs
	printf("Factoring the %d x %d stiffness matrix...\n", numRetainedDOFs, numRetainedDOFs);
	//SPOOLESSolver * solver = new SPOOLESSolver(stiffnessMatrix);

	LinearSolver * solver;

	#ifdef PARDISO_SOLVER_IS_AVAILABLE
	//int directIterative = 0;
	int numThreads = 4;
	PardisoSolver * pardisoSolver = new PardisoSolver(stiffnessMatrix, numThreads, PardisoSolver::REAL_SYM_INDEFINITE);
	pardisoSolver->FactorMatrix(stiffnessMatrix);
	solver = pardisoSolver;
	#elif defined(SPOOLES_SOLVER_IS_AVAILABLE)
	int numThreads = 4;
	if (numThreads > 1)
	  solver = new SPOOLESSolverMT(stiffnessMatrix, numThreads);
	else
	  solver = new SPOOLESSolver(stiffnessMatrix);
	#else
	solver = new CGSolver(stiffnessMatrix);
	#endif

	#ifdef PARDISO_SOLVER_IS_AVAILABLE
	// with Pardiso, we cannot parallelize the multiple solves using OpenMP, as Pardiso is not thread-safe in this way
	printf("Solving for modal derivatives using PARDISO.\n"); fflush(NULL);
	pardisoSolver->SolveLinearSystemMultipleRHS(modalDerivativesConstrained, rhsConstrained, numDeriv);
	#else
	#ifdef USE_OPENMP
	  #ifndef SPOOLES_SOLVER_IS_AVAILABLE
	    // we can parallelize the CG solver using OpenMP:
	    #pragma omp parallel for
	  #endif
	#endif
	for(int i=0; i< numDeriv; i++)
	{
	  printf("Solving for derivative #%d out of %d.\n", i + 1, numDeriv); fflush(NULL);
	  solver->SolveLinearSystem(&modalDerivativesConstrained[ELT(numRetainedDOFs, 0, i)], &rhsConstrained[ELT(numRetainedDOFs, 0, i)]);
	}
	#endif

	free(rhsConstrained);
	delete(solver);
	delete(stiffnessMatrix);

	*modalDerivatives = (double*) malloc (sizeof(double) * numDeriv * n3);

	// insert zero rows into the computed derivatives
	for(int i=0; i<numDeriv; i++)
	{
	InsertRows(n3, &modalDerivativesConstrained[numRetainedDOFs*i], 
	  &((*modalDerivatives)[n3*i]), 
	  3 * numConstrainedVertices, constrainedDOFs, oneIndexed);
	}

	free(modalDerivativesConstrained);

	// remove rigid modes from modal derivatives
	if (numRigidModes > 0)
	printf("Removing rigid modes from modal derivatives...\n");fflush(NULL);

	if (numRigidModes < 6)
	{
		double * buffer1 = (double*) malloc (sizeof(double) * n3);
		for(int i=0; i<numDeriv; i++)
		{
		  	massMatrix->MultiplyVector(&((*modalDerivatives)[n3 * i]), buffer1);
		  	for(int j=0; j<numRigidModes; j++)
		  	{
			    // rhs -= <rhs, rigid mode j> * rigid mode j
			    // rigid modes are mass-orthonormal
			    double dotp = 0.0;
			    for(int k=0; k<n3; k++)
			      dotp += buffer1[k] * Ulin[ELT(n3,k,j)]; 
			    for(int k=0; k<n3; k++)
			      (*modalDerivatives)[ELT(n3,k,i)] -= dotp * Ulin[ELT(n3,k,j)];
			}
		}
		free(buffer1);
	}
	else
	{
		RemoveSixRigidModes(numDeriv, *modalDerivatives);
	}

	//printf("Mass-normalizing modal derivatives...\n");fflush(NULL);

	// mass-normalize modal derivatives
	//for(int i=0; i < precomputationState.numDeriv; i++)
	//massMatrix->NormalizeVector(&((*modalDerivatives)[n3 * i]));

	delete(massMatrix);

	free(constrainedDOFs);
}

void FEASolve::RemoveSixRigidModes(int numVectors, double * x)
{
  int n3 = 3 * volumetricMesh->getNumVertices();

  // remove six rigid modes from rhs
  double * nullspace6 = (double*) malloc (sizeof(double) * n3 * 6);
  double * defoPos6 = (double*) malloc (sizeof(double) * n3);
  for(int i=0; i<n3/3; i++)
  {
    Vec3d restPos = volumetricMesh->getVertex(i);
    for(int j=0; j<3; j++)
      defoPos6[3*i+j] = restPos[j];
  }

  ComputeStiffnessMatrixNullspace::ComputeNullspace(n3 / 3, defoPos6, nullspace6, 1, 1);
  free(defoPos6);

  for(int i=0; i<numVectors; i++)
  {
    ComputeStiffnessMatrixNullspace::RemoveNullspaceComponent(n3 / 3, 6, nullspace6, &x[ELT(n3,0,i)]);
  }

  free(nullspace6);
}

void * FEASolve::ComputeNonLinearModes(int * code, int dataOrigin, int numNonLinearModes, double ** modes_)
{
	int n3 = 3 * volumetricMesh->getNumVertices();
	int numDataVectors = 0;

	// compute the lumped mass matrix
	SparseMatrix * massMatrix; // will be sparse n3 x n3
	GenerateMassMatrix::computeMassMatrix(volumetricMesh, &massMatrix, true); // exitCode will always be 0

	// mass-normalize modal derivatives
	double * modalDerivatives = modalDerivativesMatrix->GetMatrix();
	double * normalizedModalDerivatives = (double*) malloc (sizeof(double) * n3 * numDeriv);
	memcpy(normalizedModalDerivatives, modalDerivatives, sizeof(double) * n3 * numDeriv);

	for(int i=0; i < numDeriv; i++)
		massMatrix->NormalizeVector(&(normalizedModalDerivatives[n3 * i]));

	double * dataMatrix;
	if (dataOrigin == 0)
	{
		// use linear modes and derivatives
		// construct PCA data matrix:
		int numUsedLinearModes = rLin - numRigidModes;
		numDataVectors = numUsedLinearModes + numUsedLinearModes * (numUsedLinearModes + 1) / 2;
		printf("Number of PCA datamatrix columns: %d.\n", numDataVectors);
		dataMatrix = (double*) malloc (sizeof(double) * n3 * numDataVectors);

		printf("Generating datamatrix for SVD...\n");

		double lambda0 = frequencies[numRigidModes] * frequencies[numRigidModes];

		// scale linear modes 
		double * Ulin = linearModalMatrix->GetMatrix();
		for(int i=0; i<numUsedLinearModes; i++)
		{
		  double lambda = frequencies[numRigidModes + i] * frequencies[numRigidModes + i];
		  double factor = lambda0 / lambda;
		  for(int vertex=0; vertex< n3; vertex++)
		    dataMatrix[ELT(n3,vertex,i)] = factor * Ulin[ELT(n3,vertex,numRigidModes + i)];
		}

		// scale modal derivatives
		int pos = 0;
		for(int i=0; i<numUsedLinearModes; i++)
		{
		  for(int j=i; j<numUsedLinearModes; j++)
		  {
		    double lambdai = frequencies[numRigidModes + i] * frequencies[numRigidModes + i];
		    double lambdaj = frequencies[numRigidModes + j] * frequencies[numRigidModes + j];
		    double factor = lambda0 * lambda0 / (lambdai * lambdaj);

		    for(int vertex=0; vertex < n3; vertex++)
		      dataMatrix[ELT(n3, vertex, numUsedLinearModes + pos)] = 
		        factor * normalizedModalDerivatives[ELT(n3, vertex, pos)];
		    pos++;
		  }
		}
	}

	free(normalizedModalDerivatives);

	// do lumped-mass-PCA on dataMatrix ( n3 x numDataVectors )

	double * ones = (double*) malloc (sizeof(double) * n3);
	for(int i=0; i<n3; i++)
	ones[i] = 1.0;

	double * LTDiagonal = (double*) malloc (sizeof(double) * n3);
	massMatrix->MultiplyVector(ones, LTDiagonal);
	free(ones);
	delete(massMatrix);

	// sqrt
	for(int i=0; i<n3; i++)
	LTDiagonal[i] = sqrt(LTDiagonal[i]);

	// number of retained dimensions can't be more than num linear modes + num derivatives
	if (numComputedNonLinearModes > numDataVectors)
		numComputedNonLinearModes = numDataVectors;

	// premultiply by LT
	for(int i=0; i<n3; i++)
	for(int j=0; j < numDataVectors; j++)
	  dataMatrix[ELT(n3, i, j)] *= LTDiagonal[i];

	// do SVD on dataMatrix ( n3 x numDataVectors ), retain uiState.numComputedNonLinearModes modes

	ThresholdingSpecification thresholdingSpecification;
	thresholdingSpecification.tresholdingType = ThresholdingSpecification::numberOfModesBased;
	thresholdingSpecification.rDesired = numComputedNonLinearModes;

	int outputr;
	int matrixPCACode = 0;
	if ( ((matrixPCACode = MatrixPCA(
	&thresholdingSpecification, n3, numDataVectors, dataMatrix, &outputr)) != 0) 
	|| (outputr != numComputedNonLinearModes))
	{
		printf("Error performing SVD. Code: %d\n", matrixPCACode);
		*code = matrixPCACode;
		free(dataMatrix);
		free(LTDiagonal);
		return NULL;
	}

	// solve L^T U = V
	for(int i=0; i<n3; i++)
	for(int j=0; j < numComputedNonLinearModes; j++)
	  dataMatrix[ELT(n3, i, j)] /= LTDiagonal[i];

	free(LTDiagonal);

	// export data
	*modes_ = (double*) realloc (dataMatrix, 
	sizeof(double) * n3 * numComputedNonLinearModes);

	computationRunning = -1;

	*code = 0;
	return NULL;
}

void * FEASolve::ComputeCubicPolynomials(int * code, StVKReducedInternalForces ** newCubicPolynomial)
{
  // compute cubic polynomials
  StVKElementABCD * precomputedABCDIntegrals = StVKElementABCDLoader::load(volumetricMesh);

  int numComputationThreads = NUM_COMPUTE_THREADS;

  if (numComputationThreads == 1)
  {
    *newCubicPolynomial = new StVKReducedInternalForces(
      rNonLin,
      nonLinearModalMatrix->GetMatrix(),
      volumetricMesh, precomputedABCDIntegrals); 
  }
  else
  {
    *newCubicPolynomial = new StVKReducedInternalForcesWX(
      rNonLin,
      nonLinearModalMatrix->GetMatrix(),
      volumetricMesh, precomputedABCDIntegrals, false, numComputationThreads); 
  }

  delete(precomputedABCDIntegrals);

  computationRunning = -1;

  *code = 0;
  return NULL;
}


// save precomputation to disk
void FEASolve::SaveModalDerivatives(string path)
{
	const char * filename = path.c_str();
	int code = WriteMatrixToDisk((char*)filename, 
		3 * modalDerivativesMatrix->Getn(), 
		modalDerivativesMatrix->Getr(), 
		modalDerivativesMatrix->GetMatrix());

	if (code != 0)
	{
		cout << "Unable to save modal derivatives to file!" << endl;
		return;
	}
}

void FEASolve::SaveLinearModes(string path)
{
	const char * filename = path.c_str();
	int code = WriteMatrixToDisk((char*)filename,
		3 * linearModalMatrix->Getn(), 
		linearModalMatrix->Getr(), 
		linearModalMatrix->GetMatrix());

	if (code != 0)
	{
		cout << "Could not save linear modes to disk!" << endl;
		return;
	}
}