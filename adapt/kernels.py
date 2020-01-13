from firedrake import *
try:
    from firedrake.slate.slac.compiler import PETSC_ARCH
except ImportError:
    import os

    PETSC_ARCH = os.environ.get('PETSC_ARCH')
    PETSC_DIR = os.environ.get('PETSC_DIR')
    PETSC_ARCH = os.path.join(PETSC_DIR, PETSC_ARCH)
    if not os.path.exists(os.path.join(PETSC_ARCH, 'include/eigen3')):
        PETSC_ARCH = '/usr/local'

from adapt_utils.options import * 


__all__ = ["eigen_kernel", "get_eigendecomposition", "get_reordered_eigendecomposition",
           "set_eigendecomposition", "intersect", "anisotropic_refinement",
           "metric_from_hessian", "scale_metric", "include_dir",
           "gemv", "matscale", "singular_value_decomposition", "get_maximum_length_edge"]


include_dir = ["%s/include/eigen3" % PETSC_ARCH]


def eigen_kernel(kernel, *args, **kwargs):
    """
    Helper function to easily pass Eigen kernels to Firedrake via PyOP2.
    """
    return op2.Kernel(kernel(*args, **kwargs), kernel.__name__, cpp=True, include_dirs=include_dir)

def get_eigendecomposition(d):
    return """
#include <Eigen/Dense>

using namespace Eigen;

void get_eigendecomposition(double EVecs_[%d], double EVals_[%d], const double * M_) {
  Map<Matrix<double, %d, %d, RowMajor> > EVecs((double *)EVecs_);
  Map<Vector%dd> EVals((double *)EVals_);
  Map<Matrix<double, %d, %d, RowMajor> > M((double *)M_);
  SelfAdjointEigenSolver<Matrix<double, %d, %d, RowMajor>> eigensolver(M);
  Matrix<double, %d, %d, RowMajor> Q = eigensolver.eigenvectors();
  Vector%dd D = eigensolver.eigenvalues();
  EVecs = Q;
  EVals = D;
}
""" % (d*d, d, d, d, d, d, d, d, d, d, d, d)

def get_reordered_eigendecomposition_2d():
    return  """
#include <Eigen/Dense>

using namespace Eigen;

void get_reordered_eigendecomposition(double EVecs_[4], double EVals_[2], const double * M_) {
  Map<Matrix<double, 2, 2, RowMajor> > EVecs((double *)EVecs_);
  Map<Vector2d> EVals((double *)EVals_);
  Map<Matrix<double, 2, 2, RowMajor> > M((double *)M_);
  SelfAdjointEigenSolver<Matrix<double, 2, 2, RowMajor>> eigensolver(M);
  Matrix<double, 2, 2, RowMajor> Q = eigensolver.eigenvectors().transpose();
  Vector2d D = eigensolver.eigenvalues();
  if (fabs(D(0)) > fabs(D(1))) {
    EVecs(0,0) = Q(1,0);
    EVecs(0,1) = Q(1,1);
    EVecs(1,0) = Q(0,0);
    EVecs(1,1) = Q(0,1);
    EVals(0) = D(0);
    EVals(1) = D(1);
  } else {
    EVecs(0,0) = Q(0,0);
    EVecs(0,1) = Q(0,1);
    EVecs(1,0) = Q(1,0);
    EVecs(1,1) = Q(1,1);
    EVals(0) = D(1);
    EVals(1) = D(0);
  }
}
"""

def get_reordered_eigendecomposition(d):
    if d == 2:
        return get_reordered_eigendecomposition_2d()
    else:
        raise NotImplementedError  # TODO: 3d case

def set_eigendecomposition(d):
    return """
#include <Eigen/Dense>

using namespace Eigen;

void set_eigendecomposition(double M_[%d], const double * EVecs_, const double * EVals_) {
  Map<Matrix<double, %d, %d, RowMajor> > M((double *)M_);
  Map<Matrix<double, %d, %d, RowMajor> > EVecs((double *)EVecs_);
  Map<Vector%dd> EVals((double *)EVals_);
  M = EVecs.transpose() * EVals.asDiagonal() * EVecs;
}
""" % (d*d, d, d, d, d, d)

def intersect(d):
    return """
#include <Eigen/Dense>

using namespace Eigen;

void intersect(double M_[%d], const double * A_, const double * B_) {
  Map<Matrix<double, %d, %d, RowMajor> > M((double *)M_);
  Map<Matrix<double, %d, %d, RowMajor> > A((double *)A_);
  Map<Matrix<double, %d, %d, RowMajor> > B((double *)B_);
  SelfAdjointEigenSolver<Matrix<double, %d, %d, RowMajor>> eigensolver(A);
  Matrix<double, %d, %d, RowMajor> Q = eigensolver.eigenvectors();
  Matrix<double, %d, %d, RowMajor> D = eigensolver.eigenvalues().array().sqrt().matrix().asDiagonal();
  Matrix<double, %d, %d, RowMajor> Sq = Q * D * Q.transpose();
  Matrix<double, %d, %d, RowMajor> Sqi = Q * D.inverse() * Q.transpose();
  SelfAdjointEigenSolver<Matrix<double, %d, %d, RowMajor>> eigensolver2(Sqi.transpose() * B * Sqi);
  Q = eigensolver2.eigenvectors();
  D = eigensolver2.eigenvalues().array().max(1).matrix().asDiagonal();
  M = Sq.transpose() * Q * D * Q.transpose() * Sq;
}
""" % (d*d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d)

def anisotropic_refinement(d, direction):
    assert d in (2, 3)
    scale = 4 if d == 2 else 8
    return """
#include <Eigen/Dense>

using namespace Eigen;

void anisotropic_refinement(double A_[%d]) {
  Map<Matrix<double, %d, %d, RowMajor> > A((double *)A_);
  SelfAdjointEigenSolver<Matrix<double, %d, %d, RowMajor>> eigensolver(A);
  Matrix<double, %d, %d, RowMajor> Q = eigensolver.eigenvectors();
  Vector%dd D = eigensolver.eigenvalues();
  Array%dd D_array = D.array();
  D_array(%d) *= %f;
  A = Q * D_array.matrix().asDiagonal() * Q.transpose();
}
""" % (d*d, d, d, d, d, d, d, d, d, direction, scale)

def metric_from_hessian(d, noscale=False, op=Options()):
    scale = 'false' if noscale else 'true'
    if op.norm_order is None:
        return linf_metric_from_hessian(d, scale)
    else:
        return lp_metric_from_hessian(d, scale, op.norm_order)

def linf_metric_from_hessian(d, scale):
    return """
#include <Eigen/Dense>

using namespace Eigen;

void metric_from_hessian(double A_[%d], double * f, const double * B_)
{
  Map<Matrix<double, %d, %d, RowMajor> > A((double *)A_);
  Map<Matrix<double, %d, %d, RowMajor> > B((double *)B_);

  double mean_diag;
  int i,j;
  for (i=0; i<%d-1; i++) {
    for (j=i+1; i<%d; i++) {
      mean_diag = 0.5*(B(i,j) + B(j,i));
      B(i,j) = mean_diag;
      B(j,i) = mean_diag;
    }
  }

  SelfAdjointEigenSolver<Matrix<double, %d, %d, RowMajor>> eigensolver(B);
  Matrix<double, %d, %d, RowMajor> Q = eigensolver.eigenvectors();
  Vector%dd D = eigensolver.eigenvalues();

  /* Normalisation */
  for (i=0; i<%d; i++) D(i) = fmax(1e-10, abs(D(i)));
  if (%s) {
    double det = 1.0;
    for (i=0; i<%d; i++) det *= D(i);
    *f += sqrt(det);
  }
  A += Q * D.asDiagonal() * Q.transpose();
}
""" % (d*d, d, d, d, d, d, d, d, d, d, d, d, d, scale, d)

def lp_metric_from_hessian(d, scale, p):
    return """
#include <Eigen/Dense>

using namespace Eigen;

void metric_from_hessian(double A_[%d], double * f, const double * B_)
{
  Map<Matrix<double, %d, %d, RowMajor> > A((double *)A_);
  Map<Matrix<double, %d, %d, RowMajor> > B((double *)B_);

  double mean_diag;
  int i,j;
  for (i=0; i<%d-1; i++) {
    for (j=i+1; i<%d; i++) {
      mean_diag = 0.5*(B(i,j) + B(j,i));
      B(i,j) = mean_diag;
      B(j,i) = mean_diag;
    }
  }

  SelfAdjointEigenSolver<Matrix<double, %d, %d, RowMajor>> eigensolver(B);
  Matrix<double, %d, %d, RowMajor> Q = eigensolver.eigenvectors();
  Vector%dd D = eigensolver.eigenvalues();

  /* Normalisation */
  for (i=0; i<%d; i++) D(i) = fmax(1e-10, abs(D(i)));
  double scaling = 1.0;
  if (%s) {
    double det = 1.0;
    for (i=0; i<%d; i++) det *= D(i);
    scaling = pow(det, -1 / (2 * %d + 2));
    *f += pow(det, %d / (2 * %d + %d));
  }
  A += scaling * Q * D.asDiagonal() * Q.transpose();
}
""" % (d*d, d, d, d, d, d, d, d, d, d, d, d, d, scale, d, p, p, p, d)

def scale_metric(d, op=Options()):
    return """
#include <Eigen/Dense>

using namespace Eigen;

void scale_metric(double A_[%d])
{
  Map<Matrix<double, %d, %d, RowMajor> > A((double *)A_);
  SelfAdjointEigenSolver<Matrix<double, %d, %d, RowMajor>> eigensolver(A);
  Matrix<double, %d, %d, RowMajor> Q = eigensolver.eigenvectors();
  Vector%dd D = eigensolver.eigenvalues();
  for (int i=0; i<%d; i++) D(i) = fmin(%f, fmax(%f, abs(D(i))));
  double max_eig = fmax(D(0), D(1));
  if (%d == 3) max_eig = fmax(max_eig, D(2));
  for (int i=0; i<%d; i++) D(i) = fmax(D(i), %f * max_eig);
  A = Q * D.asDiagonal() * Q.transpose();
}
""" % (d*d, d, d, d, d, d, d, d, d, pow(op.h_min, -2), pow(op.h_max, -2), d, d, pow(op.max_anisotropy, -2))

def gemv(d, alpha=1.0, beta=0.0, tol=1e-8):
    return """
#include <Eigen/Dense>

using namespace Eigen;

void gemv(double y_[%d], const double * A_, const double * x_) {
  Map<Vector%dd> y((double *)y_);
  Map<Matrix<double, %d, %d, RowMajor> > A((double *)A_);
  Map<Vector%dd> x((double *)x_);
  double alpha = %f;
  double beta = %f;
  double tol = %f;

  if (fabs(beta) < tol) y *= beta;
  if (fabs(alpha-1.0) < tol) {
    y += A * x;
  } else {
    y += alpha * A * x;
  }
}
""" % (d, d, d, d, d, alpha, beta, tol)

def matscale(d):
    return """
#include <Eigen/Dense>

using namespace Eigen;

void matscale(double B_[%d], const double * A_, const double * alpha_) {
  Map<Matrix<double, %d, %d, RowMajor> > A((double *)A_);
  Map<Matrix<double, %d, %d, RowMajor> > B((double *)B_);
  B += *alpha_ * A;
}
""" % (d*d, d, d, d, d)

def singular_value_decomposition(d):
    return """
#include <Eigen/Dense>

using namespace Eigen;

void singular_value_decomposition(double A_[%d], const double * B_) {
  Map<Matrix<double, %d, %d, RowMajor> > A((double *)A_);
  Map<Matrix<double, %d, %d, RowMajor> > B((double *)B_);
  JacobiSVD<Matrix<double, %d, %d, RowMajor> > svd(B, ComputeFullV);

  A += svd.matrixV() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();
}""" % (d*d, d, d, d, d, d, d)

def get_maximum_length_edge(d):
    if d == 2:
        return """
for (int i=0; i<max_vector.dofs; i++) {
  int max_index = 0;
  if (edges[1] > edges[max_index]) max_index = 1;
  if (edges[2] > edges[max_index]) max_index = 2;
  max_vector[0] = vectors[2*max_index];
  max_vector[1] = vectors[2*max_index+1];
}
"""
    elif d == 3:
        raise NotImplementedError  # TODO
    else:
        raise NotImplementedError
