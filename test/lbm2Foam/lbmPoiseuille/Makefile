
# test case employed directories
LOG_DIR = log
MESH_DIR = constant/polyMesh
CODESTRAM_DIR = dynamicCode
OUT_DIR = VTK

# log files
MESHLOG = $(LOG_DIR)/mesh.log
SIMLOG = $(LOG_DIR)/lbm.log
OUTLOG = $(LOG_DIR)/vtk.log


# default objective
test: $(LOG_DIR) $(MESHLOG) $(SIMLOG) $(OUTLOG)


# main objectives
$(LOG_DIR):
	mkdir -p $(LOG_DIR)

$(MESHLOG): $(LOG_DIR)
	blockMesh > $(MESHLOG)

$(SIMLOG): $(MESHLOG)
	lbm2Foam > $(SIMLOG)

$(OUTLOG): $(SIMLOG)
	foamToVTK > $(OUTLOG)
	foamListTimes -rm


# aliases
mesh: $(MESHLOG)

simulation: $(SIMLOG)

vtk: $(OUTLOG)


# clean functions
cleanMesh:
	rm -r -f $(MESH_DIR)
	rm -f $(MESHLOG)

cleanVtk:
	rm -r -f $(OUT_DIR)
	rm -f $(OUTLOG)

cleanSim:
	rm -r -f $(CODESTRAM_DIR)
	rm -f $(SIMLOG)

clean: cleanMesh cleanSim cleanVtk
	rm -r -f $(LOG_DIR)
