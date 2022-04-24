output: main.o
	g++ main.o -o output

main.o: main.cpp
	g++ -c main.cpp

result:
	python .\post-processing\converge_plot.py -f post-processing\convergePlot.csv

visualize:
	paraview .\post-processing\solution.vtk