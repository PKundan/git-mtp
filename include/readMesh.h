#ifndef _READMESH_H
#define _READMESH_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>

#include "..\include\node.h"
#include "..\include\edge.h"
#include "..\include\cell.h"

/********************************************************************************
/ readCell Function : reads the line containing cell info(su2 format)           /
/ Inputs :- line : string , line containing cell info                           /
/           cells : vector of cells                                             /
/ Result :- cell read from the line will be appended in cells vector            /
/*******************************************************************************/

void readCell(const std::string &line,
              std::vector<Cell> &cells);

/********************************************************************************
/ readNode Function : reads the line containing node info(su2 format)           /
/ Inputs :- line : string , line containing node info                           /
/           nodes : vector of nodes                                             /
/ Result :- node read from the line will be appended in nodes vector            /
/*******************************************************************************/
void readNode(const std::string &line,
              std::vector<Node> &nodes);

/********************************************************************************
/ readEdge Function : reads the line containing boundary edge info(su2 format)  /
/ Inputs :- line : string , line containing boundary edge info                  /
/           maekerTag : tag of the boundary, ex. inlet, outlet, wall            /
/           boundaryEdges : vector of edges at the boundary                     /
/ Result :- edges read from the line will be appended in boundaryEdges vector   /
/*******************************************************************************/
void readEdge(std::string &line,
              const std::string &markerTag,
              std::vector<Edge> &boundaryEdges);

/*******************************************************************************
/ readMesh Function : Modifies input vectors nodes, cells and boundaryedges     /
/                       by reading the meshfile( *.su2 format)                  /
/ Inputs :- fileName : string , mesh file name (*.su2 format)                   /
/           nodes : vector of <Node> objects, size = 0; (initially empty)       /
/           cells : vector of <Cell> objects, size = 0; (initially empty)       /
/           boundaryEdges : vector of <Edge> objects, size=0, (initailly empty) /
/ Result :-  modification in nodes, cells and boundaryEdges vectors             /
/*******************************************************************************/
void readMesh(std::string &fileName,
              std::vector<Node> &nodes,
              std::vector<Cell> &cells,
              std::vector<Edge> &boundaryEdges);
              
#endif