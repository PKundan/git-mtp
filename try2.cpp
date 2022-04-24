#include <iostream>
#include <string>

#include "node.hpp"
#include "edge.hpp"
#include "cell.hpp"
#include "readMesh.hpp"
#include "dataStructure.hpp"

int main()
{
    std::string meshFileName = "meshfiles/rect_coarse.su2";

    std::vector<Node> nodes;
    std::vector<Edge> edges;
    std::vector<Cell> cells;

    /** Read Mesh file **/
    std::cout << "Reading MeshFile" << std::endl;
    readMesh(meshFileName, nodes, cells, edges);
    std::cout << ".\n..\n...\nDone!" << std::endl;

    std::cout << "\n\nCreating edge-based datastructure" << std::endl;
    // edgeStructure(cells, edges);
    std::cout << ".\n..\n...\nDone!" << std::endl;

    for (const auto &elem : edges)
    {
        std::cout << elem << std::endl;
    }

    // for (const auto &elem : cells)
    // {
    //     std::cout << elem << std::endl;
    // }

    return 0;
}
