#ifndef _WRITEVTK_H
#define _WRITEVTK_H

#include<iostream>
#include<string>
#include<vector>
#include<fstream>

#include "node.hpp"
#include "edge.hpp"
#include "cell.hpp"
#include "eulerFunctions.hpp"


void writeVTK(const std::string& fileName,\
              const std::vector<Node>& nodes,\
              const std::vector<Cell> cells,\
              const std::vector<std::vector<double>>& Unp1)
{
    std::fstream outFile;
    size_t nCells = cells.size();
    size_t nNodes = nodes.size();
    outFile.open(fileName, std::ios::out);
    outFile << "# vtk DataFile Version 1.0" << std::endl;
    outFile <<"Supersonic flow over a wedge"<< std::endl;
    outFile << "ASCII" << std::endl;
    outFile << "\nDATASET UNSTRUCTURED_GRID" << std::endl;
    outFile<< "POINTS "<<nNodes << " float" << std::endl;
    for(const auto& node : nodes){
        outFile << node.x << " " << node.y << " " << 1.0 << std::endl;
    }

    outFile<< "\nCELLS "<< nCells << " " << nCells*4 << std::endl;
    for(const auto& cell : cells){
        outFile << 3 << " " << cell.nodeIndices[0] <<\
                " "  << cell.nodeIndices[1] << \
               " "  << cell.nodeIndices[2] <<std::endl;
    }
    outFile<< "\nCELL_TYPES "<< nCells << std::endl;
    for(size_t i=0; i!=nCells; ++i){
        outFile << 5 << std::endl; // Triangular elements
    }

    outFile << "CELL_DATA "<< nCells << std::endl;
    outFile << "SCALARS density float" << std::endl;
    outFile << "LOOKUP_TABLE default" << std::endl;
    for(size_t i=0; i!=nCells; ++i){
        outFile << Unp1[i][0] << std::endl; // Triangular elements
    }
    outFile << "SCALARS pressure float" << std::endl;
    outFile << "LOOKUP_TABLE default" << std::endl;
    for(size_t i=0; i!=nCells; ++i){
        outFile << pressure(Unp1[i]) << std::endl; // Triangular elements
    }
    outFile << "\nVECTORS velocity float" << std::endl;
    //outFile << "LOOKUP_TABLE default" << std::endl;
    for(size_t i=0; i!=nCells; ++i){
        outFile << Unp1[i][1]/Unp1[i][0] << \
        " " << Unp1[i][2]/Unp1[i][0] << " 0.0" << std::endl; // Triangular elements
    }
}


void writeVTK_vertexCentered(const std::string& fileName,\
              const std::vector<Node>& nodes,\
              const std::vector<Cell>& cells,
              const std::vector<std::vector<double>>& Unp1)
{
    std::fstream outFile;
    size_t nNodes = nodes.size();
    size_t nCells = cells.size();
    outFile.open(fileName, std::ios::out);
    outFile << "# vtk DataFile Version 1.0" << std::endl;
    outFile <<"Supersonic flow over a wedge"<< std::endl;
    outFile << "ASCII" << std::endl;
    outFile << "\nDATASET UNSTRUCTURED_GRID" << std::endl;
    outFile<< "POINTS "<<nNodes << " float" << std::endl;
    for(const auto& node : nodes){
        outFile << node.x << " " << node.y << " " << 1.0 << std::endl;
    }

    outFile<< "\nCELLS "<< nCells << " " << nCells*4 << std::endl;
    for(const auto& cell : cells){
        outFile << 3 << " " << cell.nodeIndices[0] <<\
                " "  << cell.nodeIndices[1] << \
               " "  << cell.nodeIndices[2] <<std::endl;
    }
    outFile<< "\nCELL_TYPES "<< nCells << std::endl;
    for(size_t i=0; i!=nCells; ++i){
        outFile << 5 << std::endl; // Triangular elements
    }

    outFile << "POINT_DATA "<< nNodes << std::endl;
    outFile << "SCALARS density float" << std::endl;
    outFile << "LOOKUP_TABLE default" << std::endl;
    for(size_t i=0; i!=nNodes; ++i){
        outFile << Unp1[i][0] << std::endl; // Triangular elements
    }
    outFile << "SCALARS pressure float" << std::endl;
    outFile << "LOOKUP_TABLE default" << std::endl;
    for(size_t i=0; i!=nNodes; ++i){
        outFile << pressure(Unp1[i]) << std::endl; // Triangular elements
    }
    outFile << "\nVECTORS velocity float" << std::endl;
    //outFile << "LOOKUP_TABLE default" << std::endl;
    for(size_t i=0; i!=nNodes; ++i){
        outFile << Unp1[i][1]/Unp1[i][0] << \
        " " << Unp1[i][2]/Unp1[i][0] << " 0.0" << std::endl; // Triangular elements
    }
}
#endif