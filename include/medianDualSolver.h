#ifndef _MEDIAN_DUAL_SOLVER_H
#define _MEDIAN_DUAL_SOLVER_H

#include "..\include\node.h"
#include "..\include\edge.h"
#include "..\include\cell.h"
#include "..\include\constants.h"
#include "..\include\eulerFunctions.h"

void medianDualVolumes(const std::vector<Node> &nodes,
                       std::vector<Edge> &edges,
                       std::vector<double> &volumes);
                       
void absSurfaceVectors(const std::vector<Node> &nodes,
                                 std::vector<Edge> &edges,
                                 std::vector<Node> &surfVecs);

double medianDualTimeStep(const double &CFL,
                          const std::vector<std::vector<double>> &U,
                          const std::vector<double> &cellVolumes,
                          const std::vector<Node> &surfVecs);

#endif