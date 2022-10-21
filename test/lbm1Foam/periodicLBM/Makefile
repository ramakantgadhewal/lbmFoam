
LOG_DIR = log
MESH_DIR = constant/polyMesh
CODESTRAM_DIR = dynamicCode
OUT_DIR = VTK

MESHLOG = $(LOG_DIR)/mesh.log
SIMLOG = $(LOG_DIR)/lbm.log
OUTLOG = $(LOG_DIR)/vtk.log

test: $(LOG_DIR) $(MESHLOG) $(SIMLOG) $(OUTLOG)


$(LOG_DIR):
	mkdir -p $(LOG_DIR)


$(MESHLOG): $(LOG_DIR)
	blockMesh > $(MESHLOG)

$(SIMLOG): $(MESHLOG)
	lbm1Foam > $(SIMLOG)

$(OUTLOG): $(SIMLOG)
	foamToVTK > $(OUTLOG)
	foamListTimes -rm

mesh: $(MESHLOG)

simulation: $(SIMLOG)

vtk: $(OUTLOG)


meshClean:
	rm -r -f $(MESH_DIR)
	rm -f $(MESHLOG)

vtkClean:
	rm -r -f $(OUT_DIR)
	rm -f $(OUTLOG)

simClean:
	rm -r -f $(CODESTRAM_DIR)
	rm -f $(SIMLOG)

clean: meshClean simClean vtkClean
	rm -r -f $(LOG_DIR)
	rm -r -f $(CODESTRAM_DIR)
