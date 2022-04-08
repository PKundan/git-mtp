#ifndef _EDGE_H
#define _EDGE_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

#include "node.hpp"
#include "vectorFunctions.hpp"

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
    int left_cell_index, right_cell_index;
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

Edge::Edge(const Edge &edg)
{
    index = edg.index;
    left_cell_index = edg.left_cell_index;
    right_cell_index = edg.right_cell_index;
    nodeIndices = edg.nodeIndices;
    boundaryFlag = edg.boundaryFlag;
}

Edge Edge::operator=(Edge &edg)
{
    this->index = edg.index;
    this->left_cell_index = edg.left_cell_index;
    this->right_cell_index = edg.right_cell_index;
    this->nodeIndices = edg.nodeIndices;
    this->boundaryFlag = edg.boundaryFlag;
    return *this;
}

bool Edge::operator==(Edge const &edg) const
{
    return (edg.nodeIndices == this->nodeIndices);
}

Node Edge::nrml(const std::vector<Node> &nodes)
{
    Node n, S;
    S = SxSy(nodes);
    double S_abs = sqrt(S.x * S.x + S.y * S.y);
    n.x = S.x / S_abs;
    n.y = S.y / S_abs;
    return n;
}

Node Edge::nrml(const Node &v1, const Node &v2)
{
    Node n;
    double Sx = v2.y - v1.y;
    double Sy = v1.x - v2.x;
    double S_abs = sqrt(Sx * Sx + Sy * Sy);
    n.x = Sx / S_abs;
    n.y = Sy / S_abs;
    return n;
}

Node Edge::SxSy(const std::vector<Node> &nodes)
{
    Node S, v1 = nodes[nodeIndices[0]],
            v2 = nodes[nodeIndices[1]];
    S.x = v2.y - v1.y;
    S.y = v1.x - v2.x;
    return S;
}

Node Edge::center(const std::vector<Node> &nodes)
{
    Node c, v1 = nodes[nodeIndices[0]],
            v2 = nodes[nodeIndices[1]];
    c.x = 0.5 * (v1.x + v2.x);
    c.y = 0.5 * (v1.y + v2.y);
    return c;
}
Node Edge::center(const Node &v1, const Node &v2)
{
    Node c;
    c.x = 0.5 * (v1.x + v2.x);
    c.y = 0.5 * (v1.y + v2.y);
    return c;
}

double Edge::len(const std::vector<Node> &nodes)
{
    Node S = SxSy(nodes);
    return sqrt(S.x * S.x + S.y * S.y);
}

double Edge::len(const Node &v1, const Node &v2)
{
    double Sx = v2.y - v1.y;
    double Sy = v1.x - v2.x;
    return sqrt(Sx * Sx + Sy * Sy);
}

std::ostream &operator<<(std::ostream &os, const Edge &edg)
{
    os << edg.index << '\t' << edg.nodeIndices[0]
       << '\t' << edg.nodeIndices[1]
       << '\t' << edg.left_cell_index
       << '\t' << edg.right_cell_index
       << '\t' << edg.boundaryFlag;
    return os;
}

/* Check whether edges are 'equal' */
bool isEdgeEqual(const Edge &edg1, const Edge &edg2)
{
    return ((edg1.nodeIndices == edg2.nodeIndices));
}
/* Check whether edges are 'opposite' */
bool isEdgeOpposite(const Edge &edg1, const Edge &edg2)
{
    return (edg1.nodeIndices == vectorReverse(edg2.nodeIndices));
}

/* Check whether vector of edges contains particular edge */
bool includesEdge(const Edge &edg, const std::vector<Edge> &edges)
{
    return (std::any_of(edges.begin(),
                        edges.end(), [&edg](auto edg1)
                        { return (isEdgeEqual(edg, edg1) || isEdgeOpposite(edg, edg1)); }));
}

/* Check whether the edge is 'equal' to edge present in vector */
bool isIncludedEdgeEqual(const Edge &edg, const std::vector<Edge> &edges)
{
    return (std::any_of(edges.begin(),
                        edges.end(), [&edg](auto edg1)
                        { return (isEdgeEqual(edg, edg1)); }));
}

/* Check whether the edge is 'opposite' to edge present in vector */
bool isIncludedEdgeOpposite(const Edge &edg, const std::vector<Edge> &edges)
{
    return (std::any_of(edges.begin(),
                        edges.end(), [&edg](auto edg1)
                        { return (isEdgeOpposite(edg, edg1)); }));
}

/* Find the index of the edge in a vector */
int findEdgeIndex(const Edge &edg, std::vector<Edge> &edges)
{
    auto var = std::find(edges.begin(), edges.end(), edg);
    return (*var).index;
}

Edge findEdge(const Edge &edg, std::vector<Edge> &edges)
{
    auto var = std::find(edges.begin(), edges.end(), edg);
    return (*var);
}


/* Reverse edge: by swapping node indices and interchanging\
 left_cell_index and right_cell_index */
Edge reverseEdge(const Edge &edg)
{
    Edge edg2(edg);
    edg2.nodeIndices = vectorReverse(edg.nodeIndices);
    edg2.left_cell_index = edg.right_cell_index;
    edg2.right_cell_index = edg.left_cell_index;
    return edg2;
}

#endif
