/***------------------------------------------------------------/
/ Code for solution of 2D Euler equations on unstructured grid  /
//////////////////////////////////////////////////////////////**/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>

/** Define constants **/

// const double pi = 3.14129;
// const double GAMMA = 1.4;
// const double gp1 = 2.4;   // ga
// const double gm1 = 0.4;
// const double R_air = 287.0;
//-----------------------------------------/

/** Import other code files **/

#include "include\constants.h"
#include "include\vectorFunctions.h"
#include "include\node.h"
#include "include\edge.h"
#include "include\cell.h"
#include "include\readMesh.h"
#include "include\dataStructure.h"
#include "include\dataStructureMedianDual.h"
#include "include\eulerFunctions.h"
#include "include\Roe_solver.h"
#include "include\cell_center_solver.h"
#include "include\medianDualSolver.h"
#include "include\writeVTK.h"
//------------------------------------------/

int main()
{
    std::string method =
        "CELL_CENTER";
        // "CELL_VERTEX";
    std::string fileName =
        // "meshfiles/rect_coarse.su2";
    "meshfiles/supersonic_wedge_fine_0.05.su2";
    //"meshfiles/naca4412_exp.su2";

    std::vector<Node> nodes;
    std::vector<Edge> edges;
    std::vector<Cell> cells;
    std::vector<Cell> cellsPrimal;

    int nCells, nEdges, nNodes;

    if (method == "CELL_VERTEX")
    {
        std::vector<Edge> edgesPrimal;
        std::cout << "Reading Meshfile" << std::endl;
        readMesh(fileName, nodes, cellsPrimal, edgesPrimal);
        std::cout << ".\n..\n...Done!" << std::endl;

        nCells = nodes.size();

        std::cout << "Creating datastructure" << std::endl;
        makeEdgeDataStructVertexBased(nodes, cellsPrimal, edgesPrimal, edges);
        nNodes = nodes.size();
        nEdges = edges.size();
        std::cout << "nodes : " << nNodes << "\nCells : "
                  << nCells << "\nEdges : " << nEdges << "\n... Done !" << std::endl;
    }
    else if (method == "CELL_CENTER")
    {
        std::cout << "Reading MeshFile" << std::endl;
        readMesh(fileName, nodes, cells, edges);
        std::cout << ".\n..\n...\nDone!" << std::endl;

        std::cout << "\n\nCreating edge-based datastructure" << std::endl;
        edgeStructure2(cells, edges);

        nNodes = nodes.size();
        nEdges = edges.size();
        nCells = cells.size();
        std::cout << "nodes : " << nNodes << "\nCells : "
                  << nCells << "\nEdges : " << nEdges << "\n... Done !" << std::endl;
    }
    else
    {
        std::cout << "DEFINE METHOD CORRECTLY! (CELL_CENTER/ CELL_VERTEX)" << std::endl;
    }
    // int k = 0;
    // for (const auto &elem : edges)
    // {
    //     std::cout << k << "\t" << elem << std::endl;
    //     ++k;
    // }
    std::vector<Node> surfVecs(nCells);
    std::vector<Node> edgeNrmls(nEdges);
    std::vector<double> edgeLens(nEdges, 0.0);
    std::vector<double> cellVolumes(nCells, 0.0);
    // for (const auto &elem : cellVolumes)
    // {
    //     std::cout << elem << std::endl;
    // }
    absSurfaceVectors(nodes, edges, surfVecs);
    int k = 0;
    for (const auto &elem : surfVecs)
    {
        std::cout << k << "\t" << elem << std::endl;
        ++k;
    }

    // compute nrmls and lengths of all edges and store in edgeNrmls & edgeLens vectors resply.
    for (size_t iE = 0; iE != nEdges; ++iE)
    {
        edgeNrmls[iE] = edges[iE].nrml(nodes);
        edgeLens[iE] = edges[iE].len(nodes);
    }

    if (method == "CELL_VERTEX")
    {
        // Quantities associated with median dual control volumes
        medianDualVolumes(nodes, edges, cellVolumes);
    }
    else if (method == "CELL_CENTER")
    {
        // compute volume of all cells and store in cellVolumes vector
        for (size_t iC = 0; iC != cells.size(); ++iC)
        {
            cellVolumes[iC] = cells[iC].volume(nodes);
        }
    }

    //------------------------------------------------------
    /** ------------Upstream Conditions --------------------------- **/
    double p_upstream = 101353.0, // Pa
        T_upstream = 288.9,       // K
        theta_upstream = 0.0,     // rad
        u_upstream = 852.4,       // m/s
        M_upstream = u_upstream / sqrt(GAMMA * R_air * T_upstream);
    // rho_upstream = p_upstream / R_air / T_upstream; // Kg/m-3
    //		s_upstream = p_upstream / pow(rho_upstream, GAMMA);
    std::vector<double> UInlet(4, 0); /// Inlet velocity vector
    conservatives(M_upstream, p_upstream, T_upstream, UInlet, theta_upstream);
    //------------------------------------------------------------------

    /** Create containers for solutions at nth (Un) and (n+1)th (Unp1) time **/
    std::vector<std::vector<double>> Residuals(nCells, std::vector<double>(4, 0.));
    std::vector<std::vector<double>> Un(nCells, std::vector<double>(4, 0.));
    std::vector<std::vector<double>> Unp1(nCells, std::vector<double>(4, 0.));

    /** Intialization of the solution with upstream conditions**/
    for (auto &U : Un)
    {
        vectorAssign(UInlet, U);
    }
    //------------------------------------------------------------------------------

    /** declaration of parameters for computation **/
    double CFL = 0.8;
    int maxIter = 5000;
    double tolerance = 1e-8; /** for convergence **/

    std::fstream fout; /** Open file to write convergence history **/
    fout.open("post-processing/history.csv", std::ios::out);
    fout << "iter"
         << ",RSS" << std::endl;

    /** ------------------ solve ----------------------- **/
    double RSS = 10.0;
    double delta_t = 1e-9;
    int iter = 0;

    while (RSS > tolerance && iter < maxIter)
    { /** stopping criteria **/
        /** Residual computation by iterating over edges **/
        residualsCompute(edges, edgeLens, edgeNrmls, UInlet, Un, Residuals);
        //-------------------------------------------------
        double rho_sum = 0.0;
        /** optimal time-step computation **/
        delta_t = medianDualTimeStep(CFL, Un, cellVolumes, surfVecs);
        //-------------------------------------------------
        /** Iterate over cells **/
        for (int j = 0; j < nCells; j++)
        {
            // double vol = cells[j].volume(nodes);            /** Compute cell volume **/
            double vol = cellVolumes[j];                        /** Compute cell volume **/
            vectorMultScalar(Residuals[j], -delta_t / vol); /** R = R* -delta_t/volume **/
            Unp1[j] = vectorAdd(Un[j], Residuals[j]);       /** Unp1 = Un + R **/

            rho_sum += pow((Unp1[j][0] - Un[j][0]) / Un[j][0], 2.0); /**for RSS **/
            vectorAssign(Unp1[j], Un[j]);                            /** Un = Unp1 ; updating solution **/
        }
        RSS = sqrt(rho_sum / nCells);
        iter += 1;
        fout << iter << "," << RSS << std::endl;
        if (iter%100 == 0)
            std::cout << "iter : " << iter << "\t delta_t :" << delta_t << "\tRSS : " << RSS << std::endl;
    }
    fout.close();
    /** ------------------ solution complete -------------- **/

    /** Write solution in vtk format **/
    if (method == "CELL_VERTEX")
    {
        std::vector<Node> vertices(nCells);
        for (int k = 0; k < nCells; ++k)
        {
            vertices[k] = nodes[k];
        }
        std::string slnFilename = "post-processing/cell_vertex_solution.vtk";
        writeVTK_vertexCentered(slnFilename, vertices, cellsPrimal, Unp1);
    }
    else if (method == "CELL_CENTER")
    {
        std::string slnFilename = "post-processing/cell_center_solution.vtk";
        writeVTK(slnFilename, nodes, cells, Unp1);
    }
    //-------------------------------------------------------
    //#endif

    return 0;
}
