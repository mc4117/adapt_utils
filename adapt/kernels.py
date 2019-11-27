try:
    from firedrake.slate.slac.compiler import PETSC_ARCH
except:
    import os

    PETSC_ARCH = os.environ.get('PETSC_ARCH')
    PETSC_DIR = os.environ.get('PETSC_DIR')
    PETSC_ARCH = os.path.join(PETSC_DIR, PETSC_ARCH)
    if not os.path.exists(os.path.join(PETSC_ARCH, 'include/eigen3')):
        PETSC_ARCH = '/usr/local'

from adapt_utils.options import DefaultOptions


__all__ = ["get_eigendecomposition_kernel", "get_reordered_eigendecomposition_kernel",
           "set_eigendecomposition_kernel", "intersect_kernel", "anisotropic_refinement_kernel",
           "metric_from_hessian_kernel", "scale_metric_kernel", "include_dir"]


include_dir = ["%s/include/eigen3" % PETSC_ARCH]


def get_eigendecomposition_kernel(d):
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

get_reordered_eigendecomposition_kernel_2d = """
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

def get_reordered_eigendecomposition_kernel(d):
    if d == 2:
        return get_reordered_eigendecomposition_kernel_2d
    else:
        raise NotImplementedError  # TODO: 3d case

def set_eigendecomposition_kernel(d):
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

def intersect_kernel(d):
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

def anisotropic_refinement_kernel(d, direction):
    assert d in (2, 3)
    scale = 4 if d == 2 else 8
    return """
#include <Eigen/Dense>

using namespace Eigen;

void anisotropic(double A_[%d]) {
  Map<Matrix<double, %d, %d, RowMajor> > A((double *)A_);
  SelfAdjointEigenSolver<Matrix<double, %d, %d, RowMajor>> eigensolver(A);
  Matrix<double, %d, %d, RowMajor> Q = eigensolver.eigenvectors();
  Vector%dd D = eigensolver.eigenvalues();
  Array%dd D_array = D.array();
  D_array(%d) *= %f;
  A = Q * D_array.matrix().asDiagonal() * Q.transpose();
}
""" % (d*d, d, d, d, d, d, d, d, d, direction, scale)

def metric_from_hessian_kernel(d, noscale=False, op=DefaultOptions()):
    p = op.norm_order
    scale = 'false' if noscale or op.normalisation == 'error' else 'true'
    if p is None:
        return linf_metric_from_hessian_kernel(d, scale)
    else:
        return lp_metric_from_hessian_kernel(d, scale, p)

def linf_metric_from_hessian_kernel(d, scale):
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

def lp_metric_from_hessian_kernel(d, scale, p):
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
    *f += pow(det, %d / (2 * %d + 2));
  }
  A += scaling * Q * D.asDiagonal() * Q.transpose();
}
""" % (d*d, d, d, d, d, d, d, d, d, d, d, d, d, scale, d, p, p, p)

def scale_metric_kernel(d, op=DefaultOptions()):
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
