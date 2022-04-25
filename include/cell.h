#ifndef _CELL_H
#define _CELL_H


#include <iostream>
#include <vector>

#include "..\include\node.h"

/*****************************************************************/
// ----------------------Cell class-------------------------------/
/*****************************************************************/
/* Attributes:
        VTK_type     : cell's VTK_type ;
                        ex. triangular : 5, quadrilateral: 7 /
        index        : cell index in vector of cells
        nodeIndices  : indices of nodes forming the cell
        boundaryFlag : 1 if cell is at boundary, otherwise 0
/* Methods :
        centroid(nodes) : compute cell centroid
        volume(nodes)   : compute cell volume
/* Operators :
        <<   :  to print cell attributes on console
/******************************************************************/

class Cell
{
public:
    // int VTK_type = 0, // for triangular elements
    int index;
    std::vector<int> nodeIndices;
    bool boundaryFlag = false;

    Cell(){};

    Cell(const int &index_, //const int &VTK_type_,
         const std::vector<int> &nodeIndices_) : // VTK_type(VTK_type_),
                                                 index(index_),
                                                 nodeIndices(nodeIndices_) {}
    Node centroid(const std::vector<Node> &nodes);
    static Node centroid(const Node &v1, const Node &v2, const Node &v3);
    double volume(const std::vector<Node> &nodes);
    static double volume(const Node &v1, const Node &v2, const Node &v3);
    static double volume_quad(const Node &v1, const Node &v2, const Node &v3, const Node &v4);
    friend std::ostream &operator<<(std::ostream &os, const Cell &cell);
};

#endif