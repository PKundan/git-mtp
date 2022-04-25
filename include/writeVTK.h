#ifndef _WRITEVTK_H
#define _WRITEVTK_H

#include<iostream>
#include<string>
#include<vector>
#include<fstream>

#include "..\include\node.h"
#include "..\include\edge.h"
#include "..\include\cell.h"
#include "..\include\eulerFunctions.h"


void writeVTK(const std::string& fileName,\
              const std::vector<Node>& nodes,\
              const std::vector<Cell> cells,\
              const std::vector<std::vector<double>>& Unp1);


void writeVTK_vertexCentered(const std::string& fileName,\
              const std::vector<Node>& nodes,\
              const std::vector<Cell>& cells,
              const std::vector<std::vector<double>>& Unp1);

#endif