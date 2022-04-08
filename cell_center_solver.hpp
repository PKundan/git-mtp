#ifndef _CELL_CENTER_SOLVER_H
#define _CELL_CENTER_SOLVER_H

#include <vector>
#include <algorithm>

#include "node.hpp"
#include "edge.hpp"
#include "cell.hpp"
#include "eulerFunctions.hpp"
#include "Roe_solver.hpp"
#include "vectorFunctions.hpp"

void residualsCompute(const std::vector<Node> &nodes,
                      std::vector<Edge> &edges,
                      const std::vector<double> &UInlet,
                      const std::vector<std::vector<double>> &Un,
                      std::vector<std::vector<double>> &Residuals)
{
    for (auto &res : Residuals)
    {
        std::fill(res.begin(), res.end(), 0.);
    }

    for (auto &edg : edges)
    {
        if (edg.boundaryFlag != 0)
        {
            auto iL = edg.left_cell_index;
            std::vector<double> UR(4, 0.);
            std::vector<double> Fc(4, 0.);
            Node nrml = edg.nrml(nodes);
            double edgLength = edg.len(nodes);
            /** @Inlet **/
            if (edg.boundaryFlag == 1)
            {
                vectorAssign(UInlet, UR);
                Roe_flux_2D(Un[iL], UR, Fc, nrml.x, nrml.y);
            }
            /** @Outlet **/
            else if (edg.boundaryFlag == 2)
            {
                vectorAssign(Un[iL], UR);
                Roe_flux_2D(Un[iL], UR, Fc, nrml.x, nrml.y);
            }
            /** @Wall **/
            else if (edg.boundaryFlag == 3)
            {
                conservativesWall(Un[iL], Fc, nrml);
            }
            vectorMultScalar(Fc, edgLength);
            vectorAssignPlus(Fc, Residuals[iL]);
        }
        else
        {
            auto iL = edg.left_cell_index;
            auto iR = edg.right_cell_index;

            std::vector<double> Fc(4, 0.0);
            Node nrml = edg.nrml(nodes);
            double edgLength = edg.len(nodes);

            Roe_flux_2D(Un[iL], Un[iR], Fc, nrml.x, nrml.y); /** Roe: computes Fc ; F.n  **/
            vectorMultScalar(Fc, edgLength);                 /** Fc *= S ; multiply with edge length **/
            vectorAssignPlus(Fc, Residuals[iL]);             /** R[iL] += Fc **/
            vectorAssignNegative(Fc, Residuals[iR]);         /** R[iR] -= Fc **/
        }
    }
}

void residualsCompute(const std::vector<Edge> &edges,
                      const std::vector<double> &edgeLens,
                      const std::vector<Node> &edgeNrmls,
                      const std::vector<double> &UInlet,
                      const std::vector<std::vector<double>> &Un,
                      std::vector<std::vector<double>> &Residuals)
{
    for (auto &res : Residuals)
    {
        std::fill(res.begin(), res.end(), 0.);
    }

    // for (auto &edg : edges)
    for (size_t iE = 0; iE != edges.size(); ++iE)
    {
        if (edges[iE].boundaryFlag != 0)
        {
            auto iL = edges[iE].left_cell_index;
            std::vector<double> UR(4, 0.);
            std::vector<double> Fc(4, 0.);
            Node nrml = edgeNrmls[iE];
            double edgLength = edgeLens[iE];
            /** @Inlet **/
            if (edges[iE].boundaryFlag == 1)
            {
                vectorAssign(UInlet, UR);
                Roe_flux_2D(Un[iL], UR, Fc, nrml.x, nrml.y);
            }
            /** @Outlet **/
            else if (edges[iE].boundaryFlag == 2)
            {
                vectorAssign(Un[iL], UR);
                Roe_flux_2D(Un[iL], UR, Fc, nrml.x, nrml.y);
            }
            /** @Wall **/
            else if (edges[iE].boundaryFlag == 3)
            {
                conservativesWall(Un[iL], Fc, nrml);
            }
            vectorMultScalar(Fc, edgLength);
            vectorAssignPlus(Fc, Residuals[iL]);
        }
        else
        {
            auto iL = edges[iE].left_cell_index;
            auto iR = edges[iE].right_cell_index;

            std::vector<double> Fc(4, 0.0);
            Node nrml = edgeNrmls[iE];
            double edgLength = edgeLens[iE];

            Roe_flux_2D(Un[iL], Un[iR], Fc, nrml.x, nrml.y); /** Roe: computes Fc ; F.n  **/
            vectorMultScalar(Fc, edgLength);                 /** Fc *= S ; multiply with edge length **/
            vectorAssignPlus(Fc, Residuals[iL]);             /** R[iL] += Fc **/
            vectorAssignNegative(Fc, Residuals[iR]);         /** R[iR] -= Fc **/
        }
    }
}

double timeStep(const double &CFL,
                const std::vector<std::vector<double>> &U,
                const std::vector<Node> &nodes,
                std::vector<Edge> &edges,
                std::vector<Cell> &cells)
{
    double delta_t_min = 1.;
    for (auto &cell : cells)
    {
        int i = cell.index;
        double vol = cell.volume(nodes);
        double u = U[i][1] / U[i][0];
        double v = U[i][2] / U[i][0];
        double c = sqrt(gamma * pressure(U[i]) / U[i][0]); /// define soundSpeed function
        double S_x = 0.;
        double S_y = 0.;
        std::vector<int> nIs = cell.nodeIndices;
        nIs.push_back(cell.nodeIndices[0]);
        for (int j = 0; j != nIs.size() - 1; ++j)
        {
            double sx = nodes[nIs[j + 1]].y - nodes[nIs[j]].y;
            double sy = nodes[nIs[j]].x - nodes[nIs[j + 1]].x;
            S_x += 0.5 * std::fabs(sx);
            S_y += 0.5 * std::fabs(sy);
        }

        double eig_x = (std::fabs(u) + c) * S_x;
        double eig_y = (std::fabs(v) + c) * S_y;
        double delta_t_I = CFL * vol / (eig_x + eig_y);
        delta_t_min = delta_t_min > delta_t_I ? delta_t_I : delta_t_min;
    }
    return delta_t_min;
}

#endif