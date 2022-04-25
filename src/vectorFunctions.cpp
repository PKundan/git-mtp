#include "..\include\vectorFunctions.h"

std::vector<double> vectorMultScalar2(std::vector<double>& v1, double k) {
	for (unsigned int i = 0; i < v1.size(); i++) {
		v1[i] = k * v1[i];
	}
	return v1;
}

void vectorMultScalar(std::vector<double>& v1, double k) {
	for (unsigned int i = 0; i < v1.size(); i++) {
		v1[i] = k * v1[i];
	}
}

std::vector<double> vectorAdd(const std::vector<double>& v1, const std::vector<double>& v2) {
	//assert(v1.size() == v2.size());
	std::vector<double> v3(v1.size());
	for (unsigned int i = 0; i < v1.size(); i++) {
		v3[i] = v1[i] + v2[i];
	}
	return v3;
}

template<class T>
std::vector<T> vectorReverse(const std::vector<T>& v1){
    std::vector<T> v2(v1.size());
    int i = 0;
    for(auto it = v1.end()-1; it > v1.begin()-1; --it){
            v2[i] = *it;
            i += 1;
    }
    return v2;
}

template<class T>
T vectorSum(const std::vector<T>& v1){
    T sum = 0;
    for(const auto& elem : v1)
        sum += elem;
    return sum;
}

void vectorAssign(const std::vector<double>& v1, std::vector<double>& v2) {
    int i = 0;
	for (const auto& elem : v1){
		v2[i] = elem; i += 1;
    }
}

void vectorAssignPlus(const std::vector<double>& v1, std::vector<double>& v2) {
    int i = 0;
	for (const auto& elem : v1){
		v2[i] += elem; i += 1;
    }
}

void vectorAssignNegative(const std::vector<double>& v1, std::vector<double>& v2) {
    int i = 0;
	for (const auto& elem : v1){
		v2[i] -= elem; i += 1;
    }
}