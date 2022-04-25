#ifndef _EulerFunctions_H
#define _EulerFunctions_H

#include <vector>
#include<cmath>

#include "..\include\node.h"
#include "..\include\constants.h"

void conservatives(const double &M, const double &p, const double &T,
				   std::vector<double> &cons_, double theta);

double pressure(const std::vector<double> &U);

void conservativesWall(const std::vector<double> &U_wall,
					   std::vector<double> &FcWall,
					   const Node &nrml);

#endif
