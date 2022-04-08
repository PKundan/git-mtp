#ifndef _CELL_H
#define _CELL_H

#include "node.hpp"
#include <iostream>
#include <vector>

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

/** Implementation **/
Node Cell::centroid(const std::vector<Node> &nodes)
{
    Node c;
    Node v1 = nodes[nodeIndices[0]],
         v2 = nodes[nodeIndices[1]],
         v3 = nodes[nodeIndices[2]];
    c.x = (v1.x + v2.x + v3.x) / 3.;
    c.y = (v1.y + v2.y + v3.y) / 3.;
    return c;
}
Node Cell::centroid(const Node &v1, const Node &v2, const Node &v3)
{
    Node c;
    c.x = (v1.x + v2.x + v3.x) / 3.;
    c.y = (v1.y + v2.y + v3.y) / 3.;
    return c;
}

double Cell::volume(const std::vector<Node> &nodes)
{
    double vol_;
    Node v1 = nodes[nodeIndices[0]],
         v2 = nodes[nodeIndices[1]],
         v3 = nodes[nodeIndices[2]];
    vol_ = (v1.x - v2.x) * (v1.y + v2.y);
    vol_ += (v2.x - v3.x) * (v2.y + v3.y);
    vol_ += (v3.x - v1.x) * (v3.y + v1.y);
    return 0.5 * vol_;
}
double Cell::volume(const Node &v1, const Node &v2, const Node &v3)
{
    double vol_;
    vol_ = (v1.x - v2.x) * (v1.y + v2.y);
    vol_ += (v2.x - v3.x) * (v2.y + v3.y);
    vol_ += (v3.x - v1.x) * (v3.y + v1.y);
    return 0.5 * vol_;
}
double Cell::volume_quad(const Node &v1, const Node &v2, const Node &v3, const Node &v4)
{
    return volume(v1, v2, v3) + volume(v1, v3, v4);
}

std::ostream &operator<<(std::ostream &os, const Cell &cell)
{
    os << cell.index << '\t'; // cell.VTK_type << "\t";
    for (const auto &nIndx : cell.nodeIndices)
        os << nIndx << "\t";
    os << cell.boundaryFlag;
    return os;
}

#endif