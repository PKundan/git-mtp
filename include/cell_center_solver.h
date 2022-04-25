#ifndef _CELL_CENTER_SOLVER_H
#define _CELL_CENTER_SOLVER_H

#include <vector>
#include <algorithm>

#include "..\include\constants.h"
#include "..\include\node.h"
#include "..\include\edge.h"
#include "..\include\cell.h"
#include "..\include\eulerFunctions.h"
#include "..\include\Roe_solver.h"
#include "..\include\vectorFunctions.h"

void residualsCompute(const std::vector<Node> &nodes,
                      std::vector<Edge> &edges,
                      const std::vector<double> &UInlet,
                      const std::vector<std::vector<double>> &Un,
                      std::vector<std::vector<double>> &Residuals);

void residualsCompute(const std::vector<Edge> &edges,
                      const std::vector<double> &edgeLens,
                      const std::vector<Node> &edgeNrmls,
                      const std::vector<double> &UInlet,
                      const std::vector<std::vector<double>> &Un,
                      std::vector<std::vector<double>> &Residuals);

double timeStep(const double &CFL,
                const std::vector<std::vector<double>> &U,
                const std::vector<Node> &nodes,
                std::vector<Edge> &edges,
                std::vector<Cell> &cells,
                const std::vector<double>& cellVolumes);
#endif