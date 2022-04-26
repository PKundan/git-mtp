#include "..\include\medianDualSolver.h"

#include <vector>
#include <string>

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

void absSurfaceVectors(const std::vector<Node> &nodes,
                       std::vector<Edge> &edges,
                       std::vector<Node> &surfVecs)
{
    if (surfVecs.empty())
    {
        std::cout << "initialize surfvecs to zero" << std::endl;
    }
    for (auto &edge : edges)
    {
        int iL, iR;
        iL = edge.left_cell_index;
        iR = edge.right_cell_index;
        Node S = edge.SxSy(nodes);
        if (edge.boundaryFlag != 0)
        {
            surfVecs[iL].x += 0.5 * std::fabs(S.x);
            surfVecs[iL].y += 0.5 * std::fabs(S.y);
        }
        else
        {
            surfVecs[iL].x += 0.5 * std::fabs(S.x);
            surfVecs[iL].y += 0.5 * std::fabs(S.y);
            surfVecs[iR].x += 0.5 * std::fabs(S.x);
            surfVecs[iR].y += 0.5 * std::fabs(S.y);
        }
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