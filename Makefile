ifndef DISPLAYOBJ
DISPLAYOBJ=DISPLAYOBJ

ifndef CLEANFOLDER
CLEANFOLDER=DISPLAYOBJ
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for these utilities
DISPLAYOBJ_OBJECTS=

# the libraries these utilities depend on
# DISPLAYOBJ_LIBS=loadList integratorDense minivector polarDecomposition getopts sparseMatrix corotationalLinearFEM objMesh openGLHelper volumetricMesh imageIO glslPhong camera matrixIO lighting configFile loadList

DISPLAYOBJ_LIBS=modalMatrix reducedForceModel sceneObject lighting sceneObjectReduced reducedElasticForceModel reducedStvk integratorSparse integratorDense integrator stvk renderVolumetricMesh openGLHelper elasticForceModel isotropicHyperelasticFEM stvk forceModel sparseMatrix loadList  matrix insertRows volumetricMesh getopts graph corotationalLinearFEM clothBW polarDecomposition matrixIO    massSpringSystem objMesh configFile minivector imageIO sparseSolver

# the headers in this library
DISPLAYOBJ_HEADERS=

DISPLAYOBJ_LINK=$(addprefix -l, $(DISPLAYOBJ_LIBS)) $(ARPACK_LIB) $(SPOOLES_LIB) $(BLASLAPACK_LIB) $(PARDISO_LIB) $(FORTRAN_LIB) $(GLEW_LIB) $(STANDARD_LIBS) $(GLUI_LIB)

DISPLAYOBJ_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/computeDeformationStress/, $(DISPLAYOBJ_OBJECTS))
DISPLAYOBJ_HEADER_FILENAMES=$(addprefix $(R)/utilities/computeDeformationStress/, $(DISPLAYOBJ_HEADERS))
DISPLAYOBJ_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(DISPLAYOBJ_LIBS))
DISPLAYOBJ_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(DISPLAYOBJ_LIBS))

include $(DISPLAYOBJ_LIB_MAKEFILES)

all: $(R)/utilities/computeDeformationStress/solve #$(R)/utilities/3dfea/objMergeFiles

$(R)/utilities/computeDeformationStress/solve: $(R)/utilities/computeDeformationStress/main.cpp  /home/sourav/Experimentation/VegaFEM-v3.1/utilities/largeModalDeformationFactory/StVKReducedInternalForcesWX.cpp ./src/SimulatorApp.cpp ./src/FEASolve.cpp ./src/Optimizer.cpp $(DISPLAYOBJ_LIB_FILENAMES) $(DISPLAYOBJ_HEADER_FILENAMES)
	$(CXXLD) $(LDFLAGS) $(INCLUDE) $(GLUI_INCLUDE) $(BLASLAPACK_INCLUDE) $(DISPLAYOBJ_OBJECTS) $^ $(DISPLAYOBJ_LINK) -lThea -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options -lm -Wl,-rpath,$(GLUI_DIR)/lib $(IMAGE_LIBS) `$(WX_CONFIG) --cxxflags --libs core,base,gl` -o $@; cp $@ $(R)/utilities/bin/

# $(R)/utilities/3dfea/objMergeFiles: $(R)/utilities/3dfea/objMergeFiles.cpp $(DISPLAYOBJ_LIB_FILENAMES) $(DISPLAYOBJ_HEADER_FILENAMES)
#	$(CXXLD) $(LDFLAGS) $(INCLUDE) $(GLUI_INCLUDE) $(DISPLAYOBJ_OBJECTS) $^ $(DISPLAYOBJ_LINK) -Wl,-rpath,$(GLUI_DIR)/lib $(IMAGE_LIBS) -o $@; cp $@ $(R)/utilities/bin/

$(DISPLAYOBJ_OBJECTS_FILENAMES): %.o: %.cpp $(DISPLAYOBJ_LIB_FILENAMES) $(DISPLAYOBJ_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) `$(WX_CONFIG) --cxxflags` -c $(PARDISO_INCLUDE) $(INCLUDE) $(IMAGE_INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), DISPLAYOBJ)
clean: cleandisplayObj
endif

deepclean: cleandisplayObj

cleandisplayObj:
	$(RM) $(R)/utilities/computeDeformationStress/solve

endif

