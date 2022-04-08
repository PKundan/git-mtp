#ifndef _ROE_SOLVER_H
#define _ROE_SOLVER_H

#include <vector>
#include "eulerFunctions.hpp"


void Roe_flux_2D(const std::vector<double>& U_L, \
	const std::vector<double>& U_R, std::vector<double>& F_c, \
	const double nx, const double ny)
{
	//-----------------------------------------
	// Left state variables
	//-----------------------------------------
	double rho_L, u_L, v_L, p_L, E_L, H_L,
		u_nL, u_tL;

	rho_L = U_L[0]; E_L = U_L[3];
	u_L = U_L[1] / U_L[0]; v_L = U_L[2] / U_L[0];
	p_L = gm1 * (E_L - 0.5 * rho_L * (u_L * u_L + v_L * v_L));
	H_L = (E_L + p_L) / rho_L;

	u_nL = u_L * nx + v_L * ny;
	u_tL = -u_L * ny + v_L * nx;

	//-----------------------------------------
	// Right state variables
	//-----------------------------------------
	double rho_R, u_R, v_R, p_R, E_R, H_R,
		u_nR, u_tR;

	rho_R = U_R[0]; E_R = U_R[3];
	u_R = U_R[1] / U_R[0]; v_R = U_R[2] / U_R[0];
	p_R = gm1 * (E_R - 0.5 * rho_R * (u_R * u_R + v_R * v_R));
	H_R = (E_R + p_R) / rho_R;

	u_nR = u_R * nx + v_R * ny;
	u_tR = -u_R * ny + v_R * nx;

	//-----------------------------------------
	// Roe Averaged variables
	//-----------------------------------------
	double rho_tilde, u_n_tilde, u_t_tilde, H_tilde, a_tilde;

	rho_tilde = sqrt(rho_L * rho_R);

	double denominator = sqrt(rho_L) + sqrt(rho_R);

	u_n_tilde = ((sqrt(rho_L) * u_nL) + (sqrt(rho_R) * u_nR)) / denominator;
	u_t_tilde = ((sqrt(rho_L) * u_tL) + (sqrt(rho_R) * u_tR)) / denominator;

	H_tilde = ((sqrt(rho_L) * H_L) + (sqrt(rho_R) * H_R)) / denominator;

	a_tilde = sqrt(gm1 * (H_tilde - 0.5 * (u_n_tilde * u_n_tilde + u_t_tilde * u_t_tilde)));

	// Eigen values and eigen vectors
	double eigVal[4], eigVecMat[4][4];

	eigVal[0] = u_n_tilde - a_tilde;
	eigVal[1] = u_n_tilde;
	eigVal[2] = u_n_tilde;
	eigVal[3] = u_n_tilde + a_tilde;

	//double eigVec1[4] = {
	eigVecMat[0][0] = 1.0;
	eigVecMat[1][0] = u_n_tilde - a_tilde;
	eigVecMat[2][0] = u_t_tilde;
	eigVecMat[3][0] = H_tilde - u_n_tilde * a_tilde;

	//double eigVec2[4] = {
	eigVecMat[0][1] = 1.0;
	eigVecMat[1][1] = u_n_tilde;
	eigVecMat[2][1] = u_t_tilde;
	eigVecMat[3][1] = 0.5 * (u_n_tilde * u_n_tilde + u_t_tilde * u_t_tilde);


	//double eigVec3[4] = {
	eigVecMat[0][2] = 0.0;
	eigVecMat[1][2] = 0.0;
	eigVecMat[2][2] = 1.0;
	eigVecMat[3][2] = u_t_tilde;

	//double eigVec4[4] = {
	eigVecMat[0][3] = 1.0;
	eigVecMat[1][3] = u_n_tilde + a_tilde;
	eigVecMat[2][3] = u_t_tilde;
	eigVecMat[3][3] = H_tilde + u_n_tilde * a_tilde;// };


	//Coefficients
	double alpha_tilde[4];

	double delta_p = p_R - p_L,
		delta_u_n = u_nR - u_nL;


	alpha_tilde[0] = 0.5 * (delta_p - rho_tilde *
		a_tilde * delta_u_n) / (a_tilde * a_tilde);
	alpha_tilde[1] = (rho_R - rho_L) - (delta_p / (a_tilde * a_tilde));
	alpha_tilde[2] = rho_tilde * (u_tR - u_tL);
	alpha_tilde[3] = 0.5 * (delta_p + rho_tilde *
		a_tilde * delta_u_n) / (a_tilde * a_tilde);

	double F[4];

	F[0] = 0.5 * (rho_L * u_nL + rho_R * u_nR);
	F[1] = 0.5 * ((rho_L * u_nL * u_nL + p_L) + (rho_R * u_nR * u_nR + p_R));
	F[2] = 0.5 * (rho_L * u_nL * u_tL + rho_R * u_nR * u_tR);
	F[3] = 0.5 * ((E_L + p_L) * u_nL + (E_R + p_R) * u_nR);

	for (int j = 0; j < 4; j++) {
		for (int i = 0; i < 4; i++) {
			F[j] += -0.5 * alpha_tilde[i] * fabs(eigVal[i]) * eigVecMat[j][i];
		}
	}

	F_c[0] = F[0];
	F_c[1] = F[1] * nx - F[2] * ny;
	F_c[2] = F[1] * ny + F[2] * nx;
	F_c[3] = F[3];
}

#endif