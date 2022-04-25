
#ifndef _dataStructure_H
#define _dataStructure_H

#include <vector>

#include "..\include\node.h"
#include "..\include\edge.h"
#include "..\include\cell.h"
#include "..\include\vectorFunctions.h"

void edgeStructure2(std::vector<Cell> &cells,
                   std::vector<Edge> &boundaryEdges);

void edgeStructure(std::vector<Cell> &cells,
                   std::vector<Edge> &boundaryEdges);

#endif