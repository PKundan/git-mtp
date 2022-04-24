
#include <iostream>
#include <vector>
#include <string>

#include "vectorFunctions.hpp"
#include "readMesh.hpp"
#include "dataStructure.hpp"
#include "node.hpp"
#include "edge.hpp"
#include "cell.hpp"
#include "eulerFunctions.hpp"
#include "Roe_solver.hpp"
#include "cell_center_solver.hpp"
#include "writeVTK.hpp"

void medianDualVolumes(const std::vector<Node> &nodes,
                       std::vector<Edge> &edges,
                       std::vector<double> &volumes)
{
    for (auto &edge : edges)
    {
        int iO, iD, iL, iR;
        iO = edge.nodeIndices[0];
        iD = edge.nodeIndices[1];
        iL = edge.left_cell_index;
        iR = edge.right_cell_index;
        volumes[iL] += Cell::volume(nodes[iL], nodes[iO], nodes[iD]);
        volumes[iR] += Cell::volume(nodes[iO], nodes[iR], nodes[iD]);
    }
}

void medianDualAbsSurafceVectors(const std::vector<Node> &nodes,
                                 std::vector<Edge> &edges,
                                 std::vector<Node> &surfVecs)
{
    for (auto &edge : edges)
    {
        int iL, iR;
        iL = edge.left_cell_index;
        iR = edge.right_cell_index;
        Node S = edge.SxSy(nodes);
        surfVecs[iL].x += 0.5 * std::fabs(S.x);
        surfVecs[iL].y += 0.5 * std::fabs(S.y);
        surfVecs[iR].x += 0.5 * std::fabs(S.x);
        surfVecs[iR].y += 0.5 * std::fabs(S.y);
    }
}

double medianDualTimeStep(const double &CFL,
                          const std::vector<std::vector<double>> &U,
                          const std::vector<double> &cellVolumes,
                          const std::vector<Node> &surfVecs)
{
    double delta_t_min = 1.;
    for (int i = 0; i < cellVolumes.size(); ++i)
    {
        double vol = cellVolumes[i];
        double u = U[i][1] / U[i][0];
        double v = U[i][2] / U[i][0];
        double c = sqrt(GAMMA * pressure(U[i]) / U[i][0]);
        double S_x = surfVecs[i].x;
        double S_y = surfVecs[i].y;
        double eig_x = (std::fabs(u) + c) * S_x;
        double eig_y = (std::fabs(v) + c) * S_y;
        double delta_t_I = CFL * vol / (eig_x + eig_y);
        delta_t_min = delta_t_min > delta_t_I ? delta_t_I : delta_t_min;
    }
    return delta_t_min;
}

int main()
{
    std::string fileName = //"meshfiles/rect_coarse.su2";
        "meshfiles/supersonic_wedge_fine_0.05.su2";
		//"meshfiles/naca4412_exp.su2";

    std::vector<Node> nodes;
    std::vector<Edge> edgesPrimal;
    std::vector<Cell> cellsPrimal;

    std::vector<Edge> edges;

    std::cout << "Reading Meshfile" << std::endl;
    readMesh(fileName, nodes, cellsPrimal, edgesPrimal);
    std::cout << ".\n..\n...Done!" << std::endl;

    int nCells = nodes.size();
    std::cout << "Creating datastructure" << std::endl;
    makeEdgeDataStructVertexBased(nodes, cellsPrimal, edgesPrimal, edges);
    std::cout << ".\n..\n...Done!" << std::endl;

    // int k = 0;
    // for (auto &elem : edges){
    //     std::cout << elem << std::endl;
    //     ++k;
    // }

    size_t nEdges = edges.size();
    std::vector<Node> edgeNrmls(nEdges);
    std::vector<double> edgeLens(nEdges);

    // compute nrmls and lengths of all edges and store in edgeNrmls & edgeLens vectors resply.
    for (size_t iE = 0; iE != nEdges; ++iE)
    {
        edgeNrmls[iE] = edges[iE].nrml(nodes);
        edgeLens[iE] = edges[iE].len(nodes);
    }

    // Quantities associated with median dual control volumes
    std::vector<double> volumes(nCells, 0.0);
    medianDualVolumes(nodes, edges, volumes);

    std::vector<Node> surfVecs(nCells);
    medianDualAbsSurafceVectors(nodes, edges, surfVecs);

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
    std::vector<std::vector<double>> Residuals(nCells, std::vector<double>(4, 0));
    std::vector<std::vector<double>> Un(nCells, std::vector<double>(4, 0));
    std::vector<std::vector<double>> Unp1(nCells, std::vector<double>(4, 0));

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

    // residualsCompute(edges, edgeLens, edgeNrmls, UInlet, Un, Residuals);

    // for(const auto& elem: Residuals){
    //     std::cout << elem[0] << "\t" <<elem[1] << "\t" << elem[2] <<"\t" << elem[3] <<std::endl;
    // }

    std::fstream fout; /** Open file to write convergence history **/
    fout.open("post-processing/convergePlot.csv", std::ios::out);
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
        delta_t = medianDualTimeStep(CFL, Un, volumes, surfVecs);
        //-------------------------------------------------
        /** Iterate over cells **/
        for (int j = 0; j < nCells; j++)
        {
            // double vol = cells[j].volume(nodes);            /** Compute cell volume **/
            double vol = volumes[j];                        /** Compute cell volume **/
            vectorMultScalar(Residuals[j], -delta_t / vol); /** R = R* -delta_t/volume **/
            Unp1[j] = vectorAdd(Un[j], Residuals[j]);       /** Unp1 = Un + R **/

            rho_sum += pow((Unp1[j][0] - Un[j][0]) / Un[j][0], 2.0); /**for RSS **/
            vectorAssign(Unp1[j], Un[j]);                            /** Un = Unp1 ; updating solution **/
        }
        RSS = sqrt(rho_sum / nCells);
        iter += 1;
        fout << iter << "," << RSS << std::endl;
        std::cout << "iter : " << iter << "\t delta_t :" << delta_t << "\tRSS : " << RSS << std::endl;
    }
    fout.close();
    /** ------------------ solution complete -------------- **/

    std::vector<Node> vertices(nCells);
    for(int k=0; k<nCells; ++k){
        vertices[k] = nodes[k];
    }

    /** Write solution in vtk format **/
    writeVTK_vertexCentered("post-processing/solution.vtk", vertices, cellsPrimal, Unp1);
    // //-------------------------------------------------------
    // //#endif

    return 0;
}