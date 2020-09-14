#ifndef _EIGENDEFINITIONS_H_
#define _EIGENDEFINITIONS_H_

// in debug mode: initialize all eigen stuff as nan, to find uninitialized stuff
// (if this file is included first)
#ifndef NDEBUG
#define EIGEN_INITIALIZE_MATRICES_BY_NAN
#endif

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include <list>

typedef float scalar; // Matrix elements NOTE: float for opengl compatibility

// Fixed size matrices
typedef Eigen::Matrix<int, 2, 1> Vector2i;
typedef Eigen::Matrix<int, 3, 1> Vector3i;
typedef Eigen::Matrix<int, 4, 1> Vector4i;
typedef Eigen::Matrix<scalar, 2, 1> Vector2s;
typedef Eigen::Matrix<scalar, 3, 1> Vector3s;
typedef Eigen::Matrix<scalar, 4, 1> Vector4s;
typedef Eigen::Matrix<scalar, 2, 2> Matrix2s;
typedef Eigen::Matrix<scalar, 3, 3> Matrix3s;
typedef Eigen::Matrix<scalar, 4, 4> Matrix4s;

// Dynamic size matrices
typedef Eigen::VectorXi VectorXi;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> VectorXs;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXXs;
typedef Eigen::MatrixXi MatrixXXi;
typedef Eigen::VectorXi VectorXi;

// OpenGL specific
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixGLf;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixGLi;

// Templated matrices
template <int N, int M>
using MatrixNMs = Eigen::Matrix<scalar, N, M>; // e.g. MatrixNMs<6,3>
template <int N>
using VectorNs = Eigen::Matrix<scalar, N, 1>; // e.g. VectorNs<10>

// STL of fixed size matrices // e.g. AlignedVector<Vector<11>>
template <class T>
using AlignedVector = std::vector<T, Eigen::aligned_allocator<T>>;
template <class T>
using AlignedList = std::list<T, Eigen::aligned_allocator<T>>;

// sparse matrices
typedef Eigen::Triplet<scalar> Triplet;
typedef Eigen::SparseVector<scalar> SparseVs;
typedef Eigen::SparseMatrix<scalar> SparseXXs;
typedef Eigen::SparseMatrix<scalar, Eigen::RowMajor> SparseXXsRM;

#endif /* end of include guard: _EIGENDEFINITIONS_H_ */
