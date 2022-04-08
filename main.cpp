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
// const double gamma = 1.4;
// const double gp1 = 2.4;   // ga
// const double gm1 = 0.4;
// const double R_air = 287.0;
//-----------------------------------------/

/** Import other code files **/

#include "vectorFunctions.hpp"
#include "node.hpp"
#include "edge.hpp"
#include "cell.hpp"
#include "readMesh.hpp"
#include "dataStructure.hpp"
#include "eulerFunctions.hpp"
#include "Roe_solver.hpp"
#include "cell_center_solver.hpp"
#include "writeVTK.hpp"
//------------------------------------------/

int main()
{
    /** Specify Mesh-file name (.su2 format) **/
    std::string meshFileName = //"supersonicWedgeCoarse.su2";
        "meshfiles/supersonicWedge15_0.4.su2";
    //"supersonicWedge15_0.4_0.01.su2";
    //"gmsh01.su2";
    //"untitled.su2";
    /** Initialize containers for nodes, edges nd cells **/
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    std::vector<Cell> cells;

    /** Read Mesh file **/

    std::cout << "Reading MeshFile" << std::endl;
    readMesh(meshFileName, nodes, cells, edges);
    edgeStructure(cells, edges);
    std::cout << ".\n..\n...\nDone!" << std::endl;

    size_t nEdges = edges.size();
    std::vector<Node> edgeNrmls(nEdges);
    std::vector<double> edgeLens(nEdges);
    std::vector<double> cellVolumes(cells.size());

    // compute nrmls and lengths of all edges and store in edgeNrmls & edgeLens vectors resply.
    for (size_t iE = 0; iE != nEdges; ++iE)
    {
        edgeNrmls[iE] = edges[iE].nrml(nodes);
        edgeLens[iE] = edges[iE].len(nodes);
    }
    // compute volume of all cells and store in cellVolumes vector
    for (size_t iC = 0; iC != cells.size(); ++iC)
    {
        cellVolumes[iC] = cells[iC].volume(nodes);
    }

    /** console output edge and cell data structures **/
    // std::cout << "\nEdges" << std::endl;\
    for(const auto& edg : edges)\
        std::cout << edg << std::endl;

    //std::cout << "\nCells" << std::endl;\
    for(const auto& cell : cells)\
        std::cout << cell << std::endl;
    //#if 0
    //------------------------------------------------------
    /** ------------Upstream Conditions --------------------------- **/
    double p_upstream = 101353.0, // Pa
        T_upstream = 288.9,       // K
        theta_upstream = 0.0,     // rad
        u_upstream = 852.4,       // m/s
        M_upstream = u_upstream / sqrt(gamma * R_air * T_upstream);
    // rho_upstream = p_upstream / R_air / T_upstream; // Kg/m-3
    //		s_upstream = p_upstream / pow(rho_upstream, gamma);
    std::vector<double> UInlet(4, 0); /// Inlet velocity vector
    conservatives(M_upstream, p_upstream, T_upstream, UInlet, theta_upstream);
    //------------------------------------------------------------------

    /** Create containers for solutions at nth (Un) and (n+1)th (Unp1) time **/
    size_t nCells = cells.size();
    std::vector<std::vector<double>> Residuals(nCells, std::vector<double>(4, 0));
    std::vector<std::vector<double>> Un(nCells, std::vector<double>(4, 0));
    std::vector<std::vector<double>> Unp1(nCells, std::vector<double>(4, 0));

    /** Intialization of the solution with upstream conditions**/
    for (auto &U : Un)
    {
        vectorAssign(UInlet, U);
    } // init
    //------------------------------------------------------------------------------

    /** declaration of parameters for computation **/
    double CFL = 0.8;
    int maxIter = 5000;
    double tolerance = 1e-8; /** for convergence **/

    std::fstream fout; /** Open file to write convergence history **/
    fout.open("convergePlot.csv", std::ios::out);
    fout << "iter"
         << ",RSS" << std::endl;

    /** ------------------ solve ----------------------- **/
    double RSS = 10.0;
    double delta_t = 1e-9;
    int iter = 0;
    while (RSS > tolerance && iter < maxIter)
    { /** stopping criteria **/
        /** Residual computation by iterating over edges **/
        // residualsCompute(nodes, edges,
        //                  UInlet, Un, Residuals);
        residualsCompute(edges, edgeLens, edgeNrmls, UInlet, Un, Residuals);
        //-------------------------------------------------
        double rho_sum = 0.0;
        /** optimal time-step computation **/
        delta_t = timeStep(CFL, Un, nodes, edges, cells);
        //-------------------------------------------------
        /** Iterate over cells **/
        for (int j = 0; j < nCells; j++)
        {
            // double vol = cells[j].volume(nodes);            /** Compute cell volume **/
            double vol = cellVolumes[j];           /** Compute cell volume **/
            vectorMultScalar(Residuals[j], -delta_t / vol); /** R = R* -delta_t/volume **/
            Unp1[j] = vectorAdd(Un[j], Residuals[j]);       /** Unp1 = Un + R **/

            rho_sum += pow((Unp1[j][0] - Un[j][0]) / Un[j][0], 2.0); /**for RSS **/
            vectorAssign(Unp1[j], Un[j]);                            /** Un = Unp1 ; updating solution **/
        }
        RSS = sqrt(rho_sum / nCells);
        iter += 1;
        fout << iter << "," << RSS << std::endl;
        // std::cout << "iter : " << iter << "\t delta_t :" << delta_t << "\tRSS : " << RSS << std::endl;
    }
    fout.close();
    /** ------------------ solution complete -------------- **/

    /** Write solution in vtk format **/
    writeVTK(meshFileName + ".vtk", nodes, cells, Unp1);
    //-------------------------------------------------------
    //#endif

    return 0;
}
