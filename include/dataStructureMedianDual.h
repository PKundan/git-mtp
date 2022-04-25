#ifndef _dataStructureMD_H
#define _dataStructureMD_H

#include <vector>

#include "..\include\node.h"
#include "..\include\edge.h"
#include "..\include\cell.h"
#include "..\include\vectorFunctions.h"

/************************************************************************************************/
// EdgeDataStructure for median dual control volumes                                             /
/************************************************************************************************/
void createBoundaryEdges(std::vector<Node> &nodes,
                         const std::vector<Edge> &boundaryEdgesPrimal,
                         std::vector<Edge> &edges);

void createInteriorEdges(std::vector<Node> &nodes,
                         const std::vector<Cell> &cellsPrimal,
                         std::vector<Edge> &edgesPrimal,
                         std::vector<Edge> &edges,
                         const int &nNodesPrimal);

void makeEdgeDataStructVertexBased(
    std::vector<Node> &nodes,
    std::vector<Cell> &cellsPrimal,
    std::vector<Edge> &edgesPrimal,
    std::vector<Edge> &edges);
#endif