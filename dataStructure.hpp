
#ifndef _dataStructure_H
#define _dataStructure_H

#include <vector>

#include "node.hpp"
#include "edge.hpp"
#include "cell.hpp"
#include "vectorFunctions.hpp"

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

/************************************************************************************************/
// EdgeDataStructure for median dual control volumes                                              /
/************************************************************************************************/
void createBoundaryEdges(std::vector<Node> &nodes,
                         const std::vector<Edge> &boundaryEdgesPrimal,
                         std::vector<Edge> &edges)
{
    int k = nodes.size(); // index of appended edge midpoints and centroids
    int edgeIndex = 0;
    for (const auto &edge : boundaryEdgesPrimal)
    {
        int iNO = edge.nodeIndices[0];
        int iND = edge.nodeIndices[1];
        Node midNode = Edge::center(nodes[iNO], nodes[iND]); // edge midpoint
        midNode.flag = 2;

        nodes.push_back(midNode); // append midpoint in  nodes
        int iMid = k;
        ++k; // increment node index

        std::vector<Edge> bEdges(2);
        bEdges[0].nodeIndices = {iNO, iMid}; // create two edges from one primal boundary edge
        bEdges[1].nodeIndices = {iMid, iND};

        bEdges[0].boundaryFlag = bEdges[1].boundaryFlag = edge.boundaryFlag; // boundaryFlag
        bEdges[0].left_cell_index = iNO;
        bEdges[1].left_cell_index = iND; // assign left and righ cell indices
        bEdges[0].right_cell_index = bEdges[1].right_cell_index = -1;

        for (auto &edg : bEdges)
        {
            edg.index = edgeIndex;
            ++edgeIndex;
            edges.push_back(edg); // push newly created edges
        }
    }
}

void createInteriorEdges(std::vector<Node> &nodes,
                         const std::vector<Cell> &cellsPrimal,
                         std::vector<Edge> &edgesPrimal,
                         std::vector<Edge> &edges,
                         const int &nNodesPrimal)
{
    int edgeIndexPrimal = edgesPrimal.size();
    int edgeIndex = edges.size();
    int nodeIndex = nodes.size();
    for (int iCell = 0; iCell < cellsPrimal.size(); ++iCell)
    {
        Cell cell = cellsPrimal[iCell];
        std::vector<int> iCellNodes = cell.nodeIndices;
        iCellNodes.push_back(cell.nodeIndices[0]);
        std::vector<int> iEdgCenters(3, 0);

        for (int i = 0; i < iCellNodes.size() - 1; ++i)
        {
            Edge edg;
            edg.boundaryFlag = 0;
            edg.nodeIndices.push_back(iCellNodes[i]);
            edg.nodeIndices.push_back(iCellNodes[i + 1]);

            int indx = findEdgeIndex2(edg, edgesPrimal);
            if(indx != -1){
                iEdgCenters[i] = indx + nNodesPrimal + cellsPrimal.size();
            }
            // if (includesEdge(edg, edgesPrimal))
            // {
            //     if (isIncludedEdgeOpposite(edg, edgesPrimal))
            //     {
            //         edg.nodeIndices = vectorReverse(edg.nodeIndices);
            //     }
            //     int indx = findEdgeIndex(edg, edgesPrimal);
            //     iEdgCenters[i] = indx + nNodesPrimal + cellsPrimal.size();
            // }
            else
            {
                edg.index = edgeIndexPrimal;
                edgeIndexPrimal += 1;
                edgesPrimal.push_back(edg);
                Node mid = edg.center(nodes);
                mid.flag = 2;
                nodes.push_back(mid);
                iEdgCenters[i] = nodeIndex;
                ++nodeIndex;
            }
        }
        for (int k = 0; k != iEdgCenters.size(); ++k)
        {
            int iCentroid = nNodesPrimal + iCell;
            Edge edge;
            edge.nodeIndices.push_back(iEdgCenters[k]);
            edge.nodeIndices.push_back(iCentroid);
            edge.left_cell_index = iCellNodes[k];
            edge.right_cell_index = iCellNodes[k + 1];
            edge.index = edgeIndex;
            ++edgeIndex;
            edges.push_back(edge);
        }
    }
}

void makeEdgeDataStructVertexBased(
    std::vector<Node> &nodes,
    std::vector<Cell> &cellsPrimal,
    std::vector<Edge> &edgesPrimal,
    std::vector<Edge> &edges)
{
    int nNodesPrimal = nodes.size();
    for (auto &cell : cellsPrimal)
    {
        Node centroid = cell.centroid(nodes);
        centroid.flag = 1;
        nodes.push_back(centroid);
    }
    createBoundaryEdges(nodes, edgesPrimal, edges);
    createInteriorEdges(nodes, cellsPrimal, edgesPrimal, edges, nNodesPrimal);
}
#endif