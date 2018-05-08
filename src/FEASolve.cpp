FEASolve::FEASolve()
{
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

	nRendering = -1;
	URenderingFloat = NULL;
	r = 20;

	// TODO: set a default value for this variable by consulting the 
	// computationRunning = <some-default-value>
}

bool FEASolve::initImplicitBackwardEulerSparse()
{
	isInitImplicitBackwardEulerSparse = true;

	volumetricMesh = VolumetricMeshLoader::load(opts.in_veg_file_path.c_str());

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
	bool success = readFixedVertices(0);
	
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
	int numConstrainedDOFs = numConstrainedVertices * 3;
	int *constrainedDOFs = (int *)malloc(sizeof(int) * numConstrainedDOFs);

	int oneIndexed = 0;
	for (int i = 0, pos = 0; i < numFixedVertices; i++, pos++)
	{
		constrainedDOFs[pos*3 + 0] = (constrainedVertexList[i] - 1 + oneIndexed) * 3 + oneIndexed;
		constrainedDOFs[pos*3 + 1] = (constrainedVertexList[i] - 1 + oneIndexed) * 3 + 1 + oneIndexed;
		constrainedDOFs[pos*3 + 2] = (constrainedVertexList[i] - 1 + oneIndexed) * 3 + 2 + oneIndexed;
	}

	// initialize the integrator
	ImplicitBackwardEulerSparse * implicitBackwardEulerSparse = new
	ImplicitBackwardEulerSparse(r, opts.timestep, massMatrix, forceModel, numConstrainedDOFs, constrainedDOFs, dampingMassCoef, dampingStiffnessCoef);

	implicitBackwardEulerSparse->SetExternalForcesToZero();

	// simulation main loop
	for (int i = 0; i < opts.steps; i++)
	{
		cout << "Simulating step " << i << endl;
		if (i == 0) // set some force at the first timestep
		{
			implicitBackwardEulerSparse->SetExternalForces(f);
		}
		implicitBackwardEulerSparse->DoTimestep();
	}

	// simulation ended. calculate the deformations.
	// calculateDisplacements(implicitBackwardEulerSparse);

	// no more in use - solver ended execution.
	free(constrainedDOFs);
	// should have destructor specified - check and uncomment later.
	// free(implicitBackwardEulerSparse);
	// free(forceModel);
}

bool FEASolve::readFixedVertices(int oneIndexed)
{
	if (constrainedVertexList)
	{
		printf("Fixed vertices already read.\n");
		return false;
	}

	LoadList::load(opts.in_fixed_vertices_file_path.c_str(), &numFixedVertices, &constrainedVertexList);
	std::cout << "Number of constrained vertices read: " << numFixedVertices << std::endl;

	return true;
}

bool FEASolve::calculateDisplacements(ImplicitBackwardEulerSparse *implicitBackwardEulerSparse)
{
	u = (double*) malloc (sizeof(double) * r);
	implicitBackwardEulerSparse->GetqState(u);

	return true;
}

bool FEASolve::calculatePerVertexEnergy(ImplicitBackwardEulerSparse *implicitBackwardEulerSparse)
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