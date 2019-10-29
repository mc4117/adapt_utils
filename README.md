### Mesh adaptation in Firedrake

In this code, anisotropic and isotropic goal-oriented mesh adaptation is applied to solving the shallow water and tracer transport problems using the coastal, estuarine and ocean modelling solver provided by [Thetis][2]. Thetis is built upon the [Firedrake][1] project, which enables efficient FEM solution in Python by automatic generation of [PETSc][3] code. Anisotropic mesh adaptivity is achieved using [PRAgMaTIc][4]. A continuous adjoint solver is provided for advection-diffusion problems and the discrete adjoint code [Pyadjoint][5] can be used to generate adjoint solutions for more general problems. This is research of the Applied Modelling and Computation Group ([AMCG][6]) at Imperial College London.

### Versions

* `v1.0`: 'Anisotropic Goal-Oriented Mesh Adaptation in Firedrake': [![DOI](https://zenodo.org/badge/169627287.svg)](https://zenodo.org/badge/latestdoi/169627287)

### User instructions

* If using `v1.0`, download Firedrake, PETSc and Pragmatic using [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3250888.svg)](https://doi.org/10.5281/zenodo.3250888).

For a development version:
* Clone this repository and make it accessible to the `PYTHONPATH` environment variable.
* Set the environment variable
  ``export PETSC_CONFIGURE_OPTIONS="--download-pragmatic --with-cxx-dialect=C++11"``
  and install [Firedrake][1] with the flags ``--install thetis`` and ``--install pyadjoint``.
* Fetch and checkout the remote branch
    * ``https://github.com/jwallwork23/firedrake`` for firedrake, fork ``joe/meshadapt``
    and call ``make`` in ``firedrake/src/firedrake`` to enable pragmatic drivers.


#### For feedback, comments and questions, please email j.wallwork16@imperial.ac.uk.

[1]: http://firedrakeproject.org/ "Firedrake"
[2]: http://thetisproject.org/index.html "Thetis"
[3]: https://www.mcs.anl.gov/petsc/ "PETSc"
[4]: https://github.com/meshadaptation/pragmatic "PRAgMaTIc"
[5]: https://bitbucket.org/dolfin-adjoint/pyadjoint/src "Pyadjoint"
[6]: http://www.imperial.ac.uk/earth-science/research/research-groups/amcg/ "AMCG"
