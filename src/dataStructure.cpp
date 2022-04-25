
#include "..\include\dataStructure.h"

void edgeStructure2(std::vector<Cell> &cells,
                   std::vector<Edge> &boundaryEdges)
{
    std::vector<Edge> interiorEdges;
    int edgeIndx = 0;
    int nBEdges = boundaryEdges.size();
    for (size_t i = 0; i < cells.size(); ++i)
    {
        std::vector<int> cellNodes = cells[i].nodeIndices;
        cellNodes.push_back(cells[i].nodeIndices[0]);
        for (size_t k = 0; k < cellNodes.size() - 1; k++)
        {
            Edge edg;
            edg.nodeIndices.push_back(cellNodes[k]);
            edg.nodeIndices.push_back(cellNodes[k + 1]);
            edg.left_cell_index = i;
            // std::cout << edg << std::endl;

            int iEdgB = findEdgeIndex2(edg, boundaryEdges);

            if (iEdgB != -1)
            {
                boundaryEdges[iEdgB].left_cell_index = i;
                boundaryEdges[iEdgB].right_cell_index = -1;
                cells[i].boundaryFlag = true;
            }
            else if (edgeIndx > 1)
            {
                int iEdgI = findEdgeIndex2(edg, interiorEdges);
                if (iEdgI != -1)
                {
                    interiorEdges[iEdgI].right_cell_index = i;
                }
                else
                {
                    edg.index = edgeIndx;
                    interiorEdges.push_back(edg);
                    edgeIndx += 1;
                }
            }
            else
            {
                edg.index = edgeIndx;
                interiorEdges.push_back(edg);
                edgeIndx += 1;
            }
        }
    }
    for (auto &edg : interiorEdges)
    {
        edg.index += nBEdges;
        boundaryEdges.push_back(edg);
    }
}

void edgeStructure(std::vector<Cell> &cells,
                   std::vector<Edge> &boundaryEdges)
{
    std::vector<Edge> interiorEdges;
    int edgeIndx = 0;
    int nBEdges = boundaryEdges.size();
    for (size_t i = 0; i < cells.size(); ++i)
    {
        std::vector<int> cellNodes = cells[i].nodeIndices;
        cellNodes.push_back(cells[i].nodeIndices[0]);
        for (size_t k = 0; k < cellNodes.size() - 1; k++)
        {
            Edge edg;
            edg.nodeIndices.push_back(cellNodes[k]);
            edg.nodeIndices.push_back(cellNodes[k + 1]);
            edg.left_cell_index = i;
            // std::cout << edg << std::endl;

            int iEdgB = findEdgeIndex2(edg, boundaryEdges);

            if (iEdgB != -1)
            {
                boundaryEdges[iEdgB].left_cell_index = i;
                boundaryEdges[iEdgB].right_cell_index = -1;
                cells[i].boundaryFlag = true;
            }
            else if (edgeIndx > 1)
            {
                int iEdgI = findEdgeIndex2(edg, interiorEdges);
                if (iEdgI != -1)
                {
                    interiorEdges[iEdgI].right_cell_index = i;
                }
                else
                {
                    edg.index = edgeIndx;
                    interiorEdges.push_back(edg);
                    edgeIndx += 1;
                }
            }
            else
            {
                edg.index = edgeIndx;
                interiorEdges.push_back(edg);
                edgeIndx += 1;
            }
        }
    }
    for (auto &edg : interiorEdges)
    {
        edg.index += nBEdges;
        boundaryEdges.push_back(edg);
    }
}
