#ifndef _EDGE_H
#define _EDGE_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

#include "..\include\node.h"
#include "..\include\vectorFunctions.h"

/*****************************************************************/
// ----------------------Edge class-------------------------------/
/*****************************************************************/
/* Attributes:
        index            : edge index in vector of edges
        nodeIndices      : indices of nodes forming the edge
        left_cell_index  : index of cell left to the edge
        right_cell_index : index of right left to the edge
        boundaryFlag     : 1 if edge is at boundary, otherwise 0
/* Methods :
        nrml   : compute outward normal to the edge
        SxSy   : Compute x and y componenets of face vector S
        center : compute the edge center
        len    : compute the lenghth of edge

/* Operators :
        =    : Assignment operator
        ==   : comaprison operator
        <<   :  to print cell attributes on console
/******************************************************************/

class Edge
{
public:
    int index;
    std::vector<int> nodeIndices;
    int left_cell_index = -1, right_cell_index = -1;
    int boundaryFlag = 0;
    // std::string tag = "";

    Edge() {}

    Edge(const int &iNodeOrg, const int &iNodeDest)
    {
        nodeIndices.push_back(iNodeOrg);
        nodeIndices.push_back(iNodeDest);
    }

    Edge(const std::vector<int> &nodeIndices_, const int &index_) :
                                index(index_), nodeIndices(nodeIndices_) {}

    Edge(const Edge &edg);

    Edge operator=(Edge &edg);

    bool operator==(Edge const &edg) const;

    Node nrml(const std::vector<Node> &nodes);

    static Node nrml(const Node &v1, const Node &v2);

    Node SxSy(const std::vector<Node> &nodes);

    Node center(const std::vector<Node> &nodes);
    static Node center(const Node &v1, const Node &v2);

    double len(const std::vector<Node> &nodes);

    static double len(const Node &v1, const Node &v2);

    // inline int getIndex() { return index; }

    friend std::ostream &operator<<(std::ostream &os, const Edge &edg);
};

bool isEdgeEqual(const Edge &edg1, const Edge &edg2);
/* Check whether edges are 'opposite' */
// bool isEdgeOpposite(const Edge &edg1, const Edge &edg2);

/* Check whether vector of edges contains particular edge */
// bool includesEdge(const Edge &edg, const std::vector<Edge> &edges);

/* Check whether the edge is 'equal' to edge present in vector */
bool isIncludedEdgeEqual(const Edge &edg, const std::vector<Edge> &edges);

/* Check whether the edge is 'opposite' to edge present in vector */
// bool isIncludedEdgeOpposite(const Edge &edg, const std::vector<Edge> &edges);

/* Find the index of the edge in a vector */
int findEdgeIndex(const Edge &edg, std::vector<Edge> &edges);

int findEdgeIndex2(const Edge &edg, std::vector<Edge> &edges);

Edge findEdge(const Edge &edg, std::vector<Edge> &edges);

/* Reverse edge: by swapping node indices and interchanging\
 left_cell_index and right_cell_index */
Edge reverseEdge(Edge edg);
// Edge reverseEdge(const Edge &edg);

#endif
