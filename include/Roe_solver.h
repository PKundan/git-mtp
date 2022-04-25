#ifndef _ROE_SOLVER_H
#define _ROE_SOLVER_H

#include <vector>
#include "..\include\eulerFunctions.h"


void Roe_flux_2D(const std::vector<double>& U_L, \
	const std::vector<double>& U_R, std::vector<double>& F_c, \
	const double nx, const double ny);

#endif