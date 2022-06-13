#ifndef TYPE_H
#define TYPE_H

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

#include <vector>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef std::vector<double> Vector;
typedef std::vector<std::vector<double>> Vec2d;
#endif /* TYPE_H */