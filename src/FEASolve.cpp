#include "../includes/FEASolve.hpp"

FEASolve::FEASolve(SimulatorApp * ptr)
{
	// pointer to the simulator app object
	appPtr = ptr;

	volumetricMesh = NULL;
	pvEnergy = NULL;
	u = NULL;
	massMatrix = NULL;
	tetMesh = NULL;
	constrainedVertexList = NULL;
	deformableModel = NULL;

	// linear modes
	linearModesAvailable = frequenciesAvailable = false;
	linearModalMatrix = NULL;
	rLin = -1;
	frequencies = NULL;

	// modal derivates
	modalDerivativesAvailable = false;
	modalDerivativesMatrix = NULL;
	numDeriv = -1;

	// non-linear modes
	nonLinearModalMatrix = NULL;
	nonLinearModesAvailable = false;
	rNonLin = -1, numComputedNonLinearModes = -1;

	// cubic polynomial
	cubicPolynomialAvailable = false;
	cubicPolynomials = NULL;

	// fast reduced deformation varaibles init to NULL
	renderingModalMatrix = NULL;
	deformableObjectRenderingMeshCPU = NULL;
	deformableObjectRenderingMeshReduced = NULL;

	// other misc variables
	mass = 0.0;
	numFixedVertices = -1;
	totalStrainEnergy = maxStrainEnergy = -1;
	visualMesh = NULL;
	f = NULL;		// needs a "free" call.
	fq = springForce = NULL;	// needs a "free" call

	loaded_s_verts = false;
	oneIndexed = appPtr->opts.oneIndexed ? 1 : 0;

	nRendering = -1;
	URenderingFloat = NULL;
	r = 20;

	// TODO: set a default value for this variable by consulting the 
	// computationRunning = <some-default-value>

	cout << "\nSolver object initialized.\n";
}

bool FEASolve::initImplicitBackwardEulerSparse()
{
	isInitImplicitBackwardEulerSparse = true;

	volumetricMesh = VolumetricMeshLoader::load((appPtr->opts).in_veg_file_path.c_str());

	if (volumetricMesh == NULL)
	{
		printf("Error: failed to load mesh.\n");
		// free solver variables.
		destroyImplicitBackwardEulerSparse();
		return false;
	}
	else
		printf("Success. Number of vertices: %d . Number of elements: %d .\n", volumetricMesh->getNumVertices(), volumetricMesh->getNumElements());

	if (volumetricMesh->getElementType() == VolumetricMesh::TET)
		tetMesh = (TetMesh*) volumetricMesh; // such down-casting is safe in Vega
	else
	{
		printf("Error: not a tet mesh.\n");
		// free solver variables.
		destroyImplicitBackwardEulerSparse();
		return false;
	}

	deformableModel = new CorotationalLinearFEM(tetMesh);
	// read fixed vertices "0" indexed.
	bool success = readFixedVertices();
	
	if (!success)
	{
		printf("Could not read fixed vertices file. Please check filename.\n");
		// free solver variables.
		destroyImplicitBackwardEulerSparse();
		return false;
	}

	GenerateMassMatrix::computeMassMatrix(tetMesh, &massMatrix, true);
	r = 3 * getNumVertices();

	// allocate memory for force vector
	f = (double*) malloc (sizeof(double) * r);

	// compute the gravitational pull on every vertex.	
	volumetricMesh->computeGravity(f, 9.81, true);

	return true;
}

bool FEASolve::destroyImplicitBackwardEulerSparse()
{
	if (volumetricMesh)
	{	
		free(volumetricMesh);
		volumetricMesh = tetMesh = NULL;
	}

	if (deformableModel)
	{
		free(deformableModel);
		deformableModel = NULL;
	}

	if (constrainedVertexList)
	{
		free(constrainedVertexList);
		constrainedVertexList = NULL;
	}

	if (f)
	{
		free(f);
		f = NULL;
	}

	// reset the boolean flag for successive calls.
	isInitImplicitBackwardEulerSparse = false;

	return true;
}

bool FEASolve::runImplicitBackwardEulerSparse()
{
	if (!isInitImplicitBackwardEulerSparse)
		if	(!initImplicitBackwardEulerSparse())
		{	
			printf("Error: init method failed for implicit backward euler sparse solver!\n");
			return false;
		}

	ForceModel * forceModel = new CorotationalLinearFEMForceModel(deformableModel);

	double dampingMassCoef = 0.0; // "underwater"-like damping (here turned off)
	double dampingStiffnessCoef = 0.0; // (primarily) high-frequency damping

	// expand the constrained DOFs from the constrained vertices list. for now, all the x, y, z dirs are fixed.
	int numConstrainedDOFs = numFixedVertices * 3;
	int *constrainedDOFs = (int *)malloc(sizeof(int) * numConstrainedDOFs);

	// int oneIndexed = 0;
	for (int i = 0, pos = 0; i < numFixedVertices; i++, pos++)
	{
		constrainedDOFs[pos*3 + 0] = (constrainedVertexList[i] - 1 + oneIndexed) * 3 + oneIndexed;
		constrainedDOFs[pos*3 + 1] = (constrainedVertexList[i] - 1 + oneIndexed) * 3 + 1 + oneIndexed;
		constrainedDOFs[pos*3 + 2] = (constrainedVertexList[i] - 1 + oneIndexed) * 3 + 2 + oneIndexed;
	}

	// initialize the integrator
	ImplicitBackwardEulerSparse * implicitBackwardEulerSparse = new
	ImplicitBackwardEulerSparse(r, (appPtr->opts).in_timestep, massMatrix, forceModel, numConstrainedDOFs, constrainedDOFs, dampingMassCoef, dampingStiffnessCoef);

	implicitBackwardEulerSparse->SetExternalForcesToZero();

	// simulation main loop
	for (int i = 0; i < (appPtr->opts).in_steps; i++)
	{
		// cout << "Simulating step " << i << endl;
		if (i == 0) // set some force at the first timestep
		{
			implicitBackwardEulerSparse->SetExternalForces(f);
		}
		implicitBackwardEulerSparse->DoTimestep();
	}

	// simulation ended. calculate the deformations.
	calculateDisplacements(implicitBackwardEulerSparse);

	// no more in use - solver ended execution.
	free(constrainedDOFs);
	// should have destructor specified - check and uncomment later.
	// free(implicitBackwardEulerSparse);
	// free(forceModel);
}

bool FEASolve::readFixedVertices()
{
	if (constrainedVertexList)
	{
		printf("Fixed vertices already read.\n");
		return false;
	}

	LoadList::load((appPtr->opts).in_fixed_vertices_file_path.c_str(), &numFixedVertices, &constrainedVertexList);
	std::cout << "Number of constrained vertices read: " << numFixedVertices << std::endl;

	return true;
}

bool FEASolve::calculateDisplacements(ImplicitBackwardEulerSparse *implicitBackwardEulerSparse)
{
	if (u)
	{
		free(u); u = NULL;
	}

	u = (double*) malloc (sizeof(double) * r);
	implicitBackwardEulerSparse->GetqState(u);

	return true;
}

bool FEASolve::calculatePerVertexEnergy()
{
	if (pvEnergy)
	{
		free(pvEnergy);
		pvEnergy = NULL;
	}

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

// getters defination
int FEASolve::getNumVertices() const
{
	if (volumetricMesh)
	{
		return volumetricMesh->getNumVertices();
	}
	else if (nRendering > 0)
	{
		return nRendering;
	}
	else
	{
		cout << "Error: No model loaded yet!\n";
	}
}

int FEASolve::getNumElements() const
{
	if (volumetricMesh)
	{
		return volumetricMesh->getNumElements();
	}
	else
	{
		cout << "\nError: Surface deformation active.\n";
		return -1;
	}
}

void FEASolve::resetVector(double *vec, int num)
{
	memset(vec, 0, num);
}

FEASolve::~FEASolve()
{
	if ((appPtr->opts).in_fast)
	{
		destroyImplicitNewmarkDense();
	}
	else
	{
		destroyImplicitBackwardEulerSparse();
	}

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

double FEASolve::getAbsoluteStrainEnergy() const
{
	return totalStrainEnergy;
}

double FEASolve::getMaximumStrainEnergy() const
{
	return maxStrainEnergy;
}

void FEASolve::applyDeformation(std::string const &path, double *u, bool is_surf) const
{
	if (!is_surf)
	{
		volumetricMesh->applyDeformation(u);
		volumetricMesh->save(path.c_str());
	}
	else
	{
		visualMesh->deform(u);
  		visualMesh->save(path.c_str());
	}
}

void FEASolve::applyDeformationAndSaveVolumetricMesh(std::string const &path) const
{
	applyDeformation(path, u, false);
}

void FEASolve::applyDeformationAndSaveSurfaceMesh(std::string const &path) const
{
	applyDeformation(path, u, true);
}

double FEASolve::getMass() const
{
	// TOOD
	// mass estimate for surface mesh has to be estimated.
	return mass;
}


bool FEASolve::initImplicitNewmarkDense()
{
	isInitImplicitNewmarkDense = true;

	ReadMatrixFromDisk_((appPtr->opts).in_modal_matrix_file_path.c_str(),&nRendering,&r,&URenderingFloat);
	nRendering /= 3;
	cout << "Number of vertices (read from Modal Matrix): " << nRendering << ", r: " << r << endl;

	double * URendering = (double*) malloc (sizeof(double) * 3 * nRendering * r);
	
	for(int i=0; i < 3 * nRendering * r; i++)
		URendering[i] = URenderingFloat[i];
	free(URenderingFloat);
	
	renderingModalMatrix = new ModalMatrix(nRendering,r,URendering);
	free(URendering); // ModalMatrix made an internal copy
	
	deformableObjectRenderingMeshCPU = new SceneObjectReducedCPU((appPtr->opts).in_rendering_mesh_file_path.c_str(), renderingModalMatrix);
	deformableObjectRenderingMeshReduced = deformableObjectRenderingMeshCPU;

	// store a separate copy for viewing the simulation output.
	visualMesh = new ObjMesh((appPtr->opts).in_rendering_mesh_file_path);

	// use the fixed vertices reading module to read the support vertices
	// remember to set them to NULL after copying data.
	// vertices "0" indexed.
	bool success;
	if (appPtr->opts.num_support_regions > 0)
		success = readSupportVertices();
	
	if (!success)
	{
		printf("Could not read fixed vertices file. Please check filename.\n");
		// free solver variables.
		destroyImplicitNewmarkDense();
		return false;
	}

	f = (double*) malloc (sizeof(double) * 3 * nRendering);
	fq = (double *) malloc (sizeof(double) * r);
	springForce = (double *) malloc(sizeof(double) * 3 * nRendering);
	
	// flush the force vector to zero
	for (int i = 0; i < 3 * nRendering; i++)
	{
		f[i] = 0;
		springForce[i] = 0;
	}

	// init fq to zero
	for (int i = 0; i < r; i++)
		fq[i] = 0;
	
	// TODO: SDF based gravity force.
	double const gp_per_vertex = -(mass * 9.81 * 5) / (double)nRendering;
	for (int i = 0; i < 3 * nRendering; i += 3)
		f[i + 1] = gp_per_vertex;

	renderingModalMatrix->ProjectSingleVertex(0, f[0], f[1], f[2], fq);
	
	for (int pulledVertex = 3; pulledVertex < 3 * nRendering; pulledVertex += 3)
		renderingModalMatrix->AddProjectSingleVertex(pulledVertex/3,
           	f[pulledVertex], f[pulledVertex + 1], f[pulledVertex + 2], fq);

	return true;
}

bool FEASolve::destroyImplicitNewmarkDense()
{
	if (f)
	{
		free(f); f = NULL;
	}

	if (fq)
	{
		free(fq); fq = NULL;
	}

	if (springForce)
	{
		free(springForce); springForce = NULL;
	}

	if (loaded_s_verts)
	{
		int n = appPtr->opts.num_support_regions;
		for (int i = 0; i < n; i++)
		{
			free(supportVertices[i]);
			supportVertices[i] = NULL;
		}

		free(w);
		w = NULL;
	}

	isInitImplicitNewmarkDense = false;
}

bool FEASolve::runImplicitNewmarkDense()
{

	if (!isInitImplicitNewmarkDense)
		if	(!initImplicitNewmarkDense())
		{	
			printf("Error: init method failed for implicit newmark dense solver!\n");
			return false;
		}

	double * massMatrix = (double*) calloc(r*r,sizeof(double));

  	for(int i = 0; i < r; i++)
  	{	
    	massMatrix[ELT(r,i,i)] = 1.0;
  	}

    StVKReducedInternalForces *stVKReducedInternalForces = new StVKReducedInternalForces((appPtr->opts).in_cubic_polynomial_file_path.c_str());
    StVKReducedStiffnessMatrix *stVKReducedStiffnessMatrix = new StVKReducedStiffnessMatrix(stVKReducedInternalForces);

   	ReducedStVKForceModel *reducedStVKForceModel = new ReducedStVKForceModel(stVKReducedInternalForces, stVKReducedStiffnessMatrix);
  	ReducedLinearStVKForceModel *reducedLinearStVKForceModel = new ReducedLinearStVKForceModel(stVKReducedStiffnessMatrix);
  	ReducedForceModel *reducedForceModel = reducedStVKForceModel;

	float const newmarkBeta = 0.25;
	float const newmarkGamma = 0.5;
	float const timestep = appPtr->opts.in_timestep;
  	
  	ImplicitNewmarkDense *implicitNewmarkDense = new ImplicitNewmarkDense(r, timestep, massMatrix, reducedForceModel);

  	double * u_prev = (double *) malloc(sizeof(double) * 3 * nRendering);

  	implicitNewmarkDense->SetExternalForcesToZero();
	implicitNewmarkDense->SetTimestep(timestep);
	implicitNewmarkDense->SetNewmarkBeta(newmarkBeta);
	implicitNewmarkDense->SetNewmarkGamma(newmarkGamma);

	Timer tmr;

	double springconst = 10000;	// spring constant
	double velocitydampingconst = 5;	// velocity damping constant

	for (int i = 0; i < (appPtr->opts).in_steps; i++)
	{
		// cout << "Simulating step " << i << endl;
		if (i == 0) // set some force at the first timestep
		{
			implicitNewmarkDense->SetExternalForces(fq);
		}
		
		implicitNewmarkDense->DoTimestep();

		if (appPtr->opts.num_support_regions > 0)
		{	// flush existing forces to zero
			for (int k = 0; k < r; k++)
				fq[k] = 0;

			// // flush spring forces
			for (int vertex_i = 0; vertex_i < nRendering; vertex_i++)
			{
				springForce[3*vertex_i + 0] = 0;
				springForce[3*vertex_i + 1] = 0;
				springForce[3*vertex_i + 2] = 0;
			}

			// update the forces acting on the fixed vertices everytime.
			calculateDisplacements(implicitNewmarkDense);

			int num_regions = appPtr->opts.num_support_regions;

			for (int iter_n = 0; iter_n < num_regions; iter_n++)
			{
				int * supportVerticesList = supportVertices[iter_n];
				int numSupportVertices = num_s_verts[iter_n];

				for (int vertex_i = 0; vertex_i < numSupportVertices; vertex_i++)
				{
					int v = supportVerticesList[vertex_i] - 1;
					springForce[3*v + 0] = -springconst * u[3*v + 0] - velocitydampingconst * (u[3*v + 0] - u_prev[3*v + 0]) / timestep;
					springForce[3*v + 1] = -springconst * u[3*v + 1] - velocitydampingconst * (u[3*v + 1] - u_prev[3*v + 1]) / timestep;
					springForce[3*v + 2] = -springconst * u[3*v + 2] - velocitydampingconst * (u[3*v + 2] - u_prev[3*v + 2]) / timestep;

					u_prev[3*v + 0] = u[3*v + 0];
					u_prev[3*v + 1] = u[3*v + 1];
					u_prev[3*v + 2] = u[3*v + 2];
		 	 	}
		 	}

	 	 	// Adding up spring force and gravity force
	 	 	for (int vertex_i = 0; vertex_i < 3 * nRendering; vertex_i += 3)
	 	 	{	
	 	 		springForce[vertex_i + 0] += f[vertex_i + 0];
	 	 		springForce[vertex_i + 1] += f[vertex_i + 1];
	 	 		springForce[vertex_i + 2] += f[vertex_i + 2];
	 	 	}

	 	 	renderingModalMatrix->ProjectSingleVertex(0, springForce[0], springForce[1], springForce[2], fq);
		
			for (int pulledVertex = 3; pulledVertex < 3 * nRendering; pulledVertex += 3)
				renderingModalMatrix->AddProjectSingleVertex(pulledVertex/3,
	           		springForce[pulledVertex], springForce[pulledVertex + 1], springForce[pulledVertex + 2], fq);

	 	 	implicitNewmarkDense->SetExternalForcesToZero();
	 	 	implicitNewmarkDense->SetExternalForces(fq);
	 	}
	}

	double t = tmr.elapsed();
	cout << "Time to run the solver: " << t << endl;

	calculateDisplacements(implicitNewmarkDense);
}

bool FEASolve::calculateDisplacements(ImplicitNewmarkDense *implicitNewmarkDense)
{
	double * q = (double *) malloc (sizeof(double) * r);

	memcpy(q, implicitNewmarkDense->Getq(), sizeof(double) * r);

	deformableObjectRenderingMeshReduced->Setq(q);
  	deformableObjectRenderingMeshReduced->Compute_uUq();

  	if (u == NULL)
  		u = (double *) malloc(sizeof(double) * nRendering * 3);
  	deformableObjectRenderingMeshReduced->Getu(u);

	free(q);

	return true;
}

bool FEASolve::saveDeformationsPerVertex(std::string const &path) const
{
	int n = getNumVertices();

	ofstream out(path.c_str());
	for (int i = 0; i < 3 * n; i+=3)
	{
		out << u[i] << " " << u[i+1] << " " << u[i+2] << endl;
	}
}

void FEASolve::setMass(double const &mass_)
{
	mass = mass_;
}

bool FEASolve::readSupportVertices()
{
	if (loaded_s_verts)
	{
		std::cout << "Support vertices have already been loaded or not been freed. Please check.\n";
		return false;
	}

	int n = appPtr->opts.num_support_regions;
	
	// weight vector - importance score for the regions
	w = (double *)malloc(sizeof(double) * n);

	for (int i = 0; i < n; i++)
	{
		num_s_verts.push_back(0);
		supportVertices.push_back(NULL);

		w[i] = 0.0;	// initially none are important
		LoadList::load((appPtr->opts).support_vertices_files[i].c_str(), &num_s_verts[i], &supportVertices[i]);
	}

	return true;
}

double * FEASolve::getWeightVector(int & n_) const
{
	n_ = appPtr->opts.num_support_regions;
	return w;
}

void FEASolve::setWeightVector(double * w_)
{
	int n = appPtr->opts.num_support_regions;
	for (int i = 0; i < n; i++)
	{
		w[i] = w_[i];
	}
}

void FEASolve::updateWeightVector(double * dw)
{
	int n = appPtr->opts.num_support_regions;
	for (int i = 0; i < n; i++)
	{
		w[i] = dw[i];
	}
}