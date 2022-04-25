#ifndef _VECTORFUNCTIONS_H
#define _VECTORFUNCTIONS_H

#include <vector>

std::vector<double> vectorMultScalar2(std::vector<double> &v1, double k);

void vectorMultScalar(std::vector<double> &v1, double k);

std::vector<double> vectorAdd(const std::vector<double> &v1, const std::vector<double> &v2);

template <typename T>
std::vector<T> vectorReverse(const std::vector<T> &v1);

template <class T>
T vectorSum(const std::vector<T> &v1);

void vectorAssign(const std::vector<double> &v1, std::vector<double> &v2);

void vectorAssignPlus(const std::vector<double> &v1, std::vector<double> &v2);

void vectorAssignNegative(const std::vector<double> &v1, std::vector<double> &v2);

#endif