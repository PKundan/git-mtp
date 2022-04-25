#
CC=g++
CFLAGS= -c
Objs = main.o readMesh.o dataStructure.o datastructureMedianDual.o\
		 cell_center_solver.o medianDualSolver.o writeVTK.o edge.o cell.o\
		  Roe_solver.o eulerFunctions.o node.o vectorFunctions.o 

all : output
 
#output: main.o readMesh.o dataStructure.o datastructureMedianDual.o cell_center_solver.o writeVTK.o edge.o cell.o Roe_solver.o eulerFunctions.o node.o vectorFunctions.o 
#	$(CC) main.o readMesh.o dataStructure.o datastructureMedianDual.o cell_center_solver.o writeVTK.o edge.o cell.o Roe_solver.o eulerFunctions.o node.o vectorFunctions.o -o output

output : $(Objs)
	$(CC) $(Objs) -o output

main.o : main.cpp
	$(CC) $(CFLAGS) main.cpp
	
readMesh.o : src\readMesh.cpp include\readMesh.h
	$(CC) $(CFLAGS) src\readMesh.cpp

dataStructure.o : src\dataStructure.cpp include\dataStructure.h
	$(CC) $(CFLAGS) src\dataStructure.cpp

datastructureMedianDual.o : src\dataStructureMedianDual.cpp include\dataStructureMedianDual.h
	$(CC) $(CFLAGS) src\dataStructureMedianDual.cpp
	
cell_center_solver.o : src\cell_center_solver.cpp include\cell_center_solver.h
	$(CC) $(CFLAGS) src\cell_center_solver.cpp

medianDualSolver.o : src\medianDualSolver.cpp include\medianDualSolver.h
	$(CC) $(CFLAGS) src\medianDualSolver.cpp
	
writeVTK.o : src\writeVTK.cpp include\writeVTK.h
	$(CC) $(CFLAGS) src\writeVTK.cpp
	
edge.o : src\edge.cpp include\edge.h
	$(CC) $(CFLAGS) src\edge.cpp
	
cell.o : src\cell.cpp include\cell.h
	$(CC) $(CFLAGS) src\cell.cpp

Roe_solver.o : src\Roe_solver.cpp include\Roe_solver.h
	$(CC) $(CFLAGS) src\Roe_solver.cpp

eulerFunctions.o : src\eulerFunctions.cpp include\eulerFunctions.h
	$(CC) $(CFLAGS) src\eulerFunctions.cpp
	
node.o : src\node.cpp include\node.h
	$(CC) $(CFLAGS) src\node.cpp

vectorFunctions.o : src\vectorFunctions.cpp include\vectorFunctions.h
	$(CC) $(CFLAGS) src\vectorFunctions.cpp
	
clean:
	rm -rf *.o output

result:
	python .\post-processing\converge_plot.py -f post-processing\convergePlot.csv

visualize:
	paraview .\post-processing\solution.vtk