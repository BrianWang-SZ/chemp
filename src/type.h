#ifndef TYPE_H
#define TYPE_H

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

#include <vector>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

template <class T>;
using Vec = std::vector<T>;
using Vec2d = std::vector<>
#endif /* TYPE_H */