#include "..\include\eulerFunctions.h"

void conservatives(const double &M, const double &p, const double &T,
				   std::vector<double> &cons_, double theta = 0.0)
{
	double rho = p / T / R_air;
	double V = M * sqrt(GAMMA * R_air * T);
	double u = V * cos(theta);
	double v = V * sin(theta);

	cons_[0] = rho;
	cons_[1] = rho * u;
	cons_[2] = rho * v;
	cons_[3] = rho * ((R_air / gm1) * T + 0.5 * (u * u + v * v));
}

double pressure(const std::vector<double> &U)
{
	return gm1 * (U[3] -
				  0.5 / U[0] * (U[1]*U[1] + U[2]*U[2]));
}

void conservativesWall(const std::vector<double> &U_wall,
					   std::vector<double> &FcWall,
					   const Node &nrml)
{
	double pWall = pressure(U_wall);
	FcWall = {0., 0., 0., 0.};
	FcWall[1] = pWall * nrml.x;
	FcWall[2] = pWall * nrml.y;
}


