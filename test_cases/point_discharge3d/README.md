## "Point Discharge with Diffusion 3D"

Extension of "Point Discharge with Diffusion" test case (as described in TELEMAC-2D validation
document version 7.0 [1]) to the 3D case. Again, the source is offset by one unit away from the
boundary. 3D extension as described in [2].


### Parameters:
  * Domain: [0, 50]x[0, 10]x[0, 10]
  * Fluid velocity: (1.0, 0.0)
  * Diffusivity coefficient: 0.1
  * Source location: (2.0, 5.0, 5.0)
  * Delta function parametrisation:
    * Centred case: 0.07980
    * Offset case: 0.07972

### Boundary conditions:
  * Dirichlet zero on inflow x=0
  * Outflow condition on x=50
  * Neumann zero elsewhere

### Quantity of interest:
Integral of tracer concentration over region of interest, given by a sphere of radius 0.5, centred at
  * Centred case: (20.0, 5.0, 5.0)
  * Offset case: (20.0, 7.5, 7.5)


[1] Riadh A., Cedric G., Jean M. "TELEMAC modeling system: 2D hydrodynamics TELEMAC-2D software
    release 7.0 user manual." Paris:  R&D, Electricite de France, p. 134, 2014.

[2]  Wallwork, J., Barral, N., Ham, D., Piggott, M.: Anisotropic goal-oriented mesh adaptation in
     Firedrake. In: 28th Intl Meshing Roundtable (2019). (to appear)