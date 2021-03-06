from thetis import *
from thetis.configuration import *

from adapt_utils.swe.options import ShallowWaterOptions


__all__ = ["MorphOptions"]


class MorphOptions(ShallowWaterOptions):
    """
    Parameter class for general morphological problems.
    """

    def __init__(self, **kwargs):
        self.slope_eff = False
        self.angle_correction = False
        self.convective_vel_flag = True
        self.wetting_and_drying = False
        self.conservative = False
        self.depth_integrated = False
        self.suspended = True
        self.bedload = True
        self.implicit_source = False
        self.fixed_tracer = None
        self.solve_tracer = True

        super(MorphOptions, self).__init__(**kwargs)

    """
    def set_tracer_init(self, fs):
        if self.fixed_tracer is not None:
            return self.fixed_tracer
        else:
            if self.conservative:
                tracer_init = project(self.depth*self.ceq/self.coeff, fs)
            else:
                divisor = interpolate(self.ceq/self.coeff, self.ceq.function_space())
                tracer_init = project(divisor, fs)
            return tracer_init
    """

    def set_up_suspended(self, mesh, tracer=None):
        P1 = FunctionSpace(mesh, "CG", 1)
        P1DG = FunctionSpace(mesh, "DG", 1)
        P1_vec = VectorFunctionSpace(mesh, "CG", 1)
        P1DG_vec = VectorFunctionSpace(mesh, "DG", 1)

        self.viscosity_ref = Constant(10**(-6))

        R = Constant(2650/1000 - 1)
        self.dstar = Constant(self.average_size*((self.g*R)/(self.viscosity_ref**2))**(1/3))
        if max(self.dstar.dat.data[:] < 1):
            print('ERROR: dstar value less than 1')
        elif max(self.dstar.dat.data[:] < 4):
            self.thetacr = Constant(0.24*(self.dstar**(-1)))
        elif max(self.dstar.dat.data[:] < 10):
            self.thetacr = Constant(0.14*(self.dstar**(-0.64)))
        elif max(self.dstar.dat.data[:] < 20):
            self.thetacr = Constant(0.04*(self.dstar**(-0.1)))
        elif max(self.dstar.dat.data[:] < 150):
            self.thetacr = Constant(0.013*(self.dstar**(0.29)))
        else:
            self.thetacr = Constant(0.055)

        self.taucr = Constant((2650-1000)*self.gravity*self.average_size*self.thetacr)

        if not hasattr(self, "settling_velocity"):
            if self.average_size <= 100*(10**(-6)):
                self.settling_velocity = Constant(9.81*(self.average_size**2)*((2650/1000)-1)/(18*self.viscosity_ref))
            elif self.average_size <= 1000*(10**(-6)):
                self.settling_velocity = Constant((10*self.viscosity_ref/self.average_size)*(sqrt(1 + 0.01*((((2650/1000) - 1)*9.81*(self.average_size**3))/(self.viscosity_ref**2)))-1))
            else:
                self.settling_velocity = Constant(1.1*sqrt(9.81*self.average_size*((2650/1000) - 1)))

        self.u_cg = interpolate(self.uv_d, P1_vec)
        self.horizontal_velocity = interpolate(self.u_cg[0], P1)
        self.vertical_velocity = interpolate(self.u_cg[1], P1)
        self.elev_cg = interpolate(self.eta_d, P1)
        
        self.unorm = ((self.horizontal_velocity**2) + (self.vertical_velocity**2))

        if self.t_old.dat.data[:] == 0.0:
            self.set_bathymetry(P1)

        if self.wetting_and_drying:
            H = interpolate(self.elev_cg + self.bathymetry, P1)
            self.depth = interpolate(H + (0.5 * (sqrt(H ** 2 + self.wetting_and_drying_alpha ** 2) - H)), P1)
        else:
            self.depth = interpolate(self.elev_cg + self.bathymetry, P1)

        self.hc = conditional(self.depth > 0.001, self.depth, 0.001)
        self.aux = conditional(11.036*self.hc/self.ks > 1.001, 11.036*self.hc/self.ks, 1.001)
        self.qfc = 2/(ln(self.aux)/0.4)**2

        self.TOB = interpolate(1000*0.5*self.qfc*self.unorm, P1)

        # skin friction coefficient
        self.cfactor = interpolate(self.get_cfactor(), P1DG)

        # mu - ratio between skin friction and normal friction
        self.mu = interpolate(conditional(self.qfc > 0, self.cfactor/self.qfc, 0), P1DG)

        self.a = Constant((self.ks)/2)
        
        self.B = interpolate(conditional(self.a > self.depth, Constant(1.0), self.a/self.depth), P1DG)
        self.ustar = sqrt(0.5*self.qfc*self.unorm)
        self.exp1 = conditional((conditional((self.settling_velocity/(0.4*self.ustar)) - 1 > 0, (self.settling_velocity/(0.4*self.ustar)) - 1, -(self.settling_velocity/(0.4*self.ustar)) + 1)) > 10**(-4), conditional((self.settling_velocity/(0.4*self.ustar)) - 1 > 3, 3, (self.settling_velocity/(0.4*self.ustar))-1), 0)
        self.coefftest = conditional((conditional((self.settling_velocity/(0.4*self.ustar)) - 1 > 0, (self.settling_velocity/(0.4*self.ustar)) - 1, -(self.settling_velocity/(0.4*self.ustar)) + 1)) > 10**(-4), self.B*(1-self.B**self.exp1)/self.exp1, -self.B*ln(self.B))
        if self.wetting_and_drying:
            self.coeff = interpolate(conditional(conditional(self.coefftest>10**(-12), 1/self.coefftest, 10**12)>1, conditional(self.coefftest>10**(-12), 1/self.coefftest, 10**12), 1), P1DG)
        else:
            self.coeff = interpolate(conditional(self.coefftest > 0, 1/self.coefftest, 0), P1DG)
            
        # erosion flux - for vanrijn
        s0 = (conditional(1000*0.5*self.qfc*self.unorm*self.mu > 0, 1000*0.5*self.qfc*self.unorm*self.mu, 0) - self.taucr)/self.taucr
        self.ceq = interpolate(0.015*(self.average_size/self.a) * ((conditional(s0 < 0, 0, s0))**(1.5))/(self.dstar**0.3), P1DG)

        if self.conservative:
            self.tracer_init_value = Constant(self.depth.at([0, 0])*self.ceq.at([0, 0])/self.coeff.at([0, 0]))
        else:
            self.tracer_init_value = Constant(self.ceq.at([0, 0])/self.coeff.at([0, 0]))

        if tracer is None:
            self.tracer_init = self.set_tracer_init(self.P1DG)
        else:
            self.tracer_init = project(tracer, self.P1DG)

        self.depo, self.ero = self.set_source_tracer(P1DG, solver_obj=None, init=True)

        if self.conservative:
            if self.depth_integrated:
                self.depth_int_sink = interpolate(self.depo/self.depth, P1DG)
                self.depth_int_source = interpolate(self.ero, P1DG)
            else:
                self.sink = interpolate(self.depo/(self.depth**2), P1DG)
                self.source = interpolate(self.ero/self.depth, P1DG)
        else:
            if self.implicit_source:
                self.sink = interpolate(self.depo/self.depth, P1DG)
                self.source = interpolate(self.ero/self.depth, P1DG)
            else:
                if self.t_old.dat.data[:] == 0.0:
                    self.source = interpolate((-(self.depo*self.tracer_init) + self.ero)/self.depth, P1DG)
                    self.sink = None
                else:
                    self.source = interpolate((-(self.depo*tracer) + self.ero)/self.depth, P1DG)

        if self.t_old.dat.data[:] == 0.0:
            if self.conservative:
                self.qbsourcedepth = interpolate(-(self.depo*self.tracer_init/self.depth) + self.ero, P1DG)
            else:
                self.qbsourcedepth = interpolate(-(self.depo*self.tracer_init) + self.ero, P1DG)
        else:
            if self.conservative:
                self.qbsourcedepth = interpolate(-(self.depo*tracer/self.depth) + self.ero, P1DG)
            else:
                self.qbsourcedepth = interpolate(-(self.depo*tracer) + self.ero, P1DG)

        if self.convective_vel_flag:
            # correction factor to advection velocity in sediment concentration equation

            self.Bconv = conditional(self.depth > 1.1*self.ksp, self.ksp/self.depth, self.ksp/(1.1*self.ksp))
            self.Aconv = conditional(self.depth > 1.1*self.a, self.a/self.depth, self.a/(1.1*self.a))

            # take max of value calculated either by ksp or depth
            self.Amax = conditional(self.Aconv > self.Bconv, self.Aconv, self.Bconv)

            self.r1conv = 1 - (1/0.4)*conditional(self.settling_velocity/self.ustar < 1, self.settling_velocity/self.ustar, 1)

            self.Ione = conditional(self.r1conv > 10**(-8), (1 - self.Amax**self.r1conv)/self.r1conv, conditional(self.r1conv < - 10**(-8), (1 - self.Amax**self.r1conv)/self.r1conv, ln(self.Amax)))

            self.Itwo = conditional(self.r1conv > 10**(-8), -(self.Ione + (ln(self.Amax)*(self.Amax**self.r1conv)))/self.r1conv, conditional(self.r1conv < - 10**(-8), -(self.Ione + (ln(self.Amax)*(self.Amax**self.r1conv)))/self.r1conv, -0.5*ln(self.Amax)**2))

            self.alpha = -(self.Itwo - (ln(self.Amax) - ln(30))*self.Ione)/(self.Ione * ((ln(self.Amax) - ln(30)) + 1))

            # final correction factor
            self.corrective_velocity_factor = Function(self.P1DG).interpolate(conditional(conditional(self.alpha > 1, 1, self.alpha) < 0, 0, conditional(self.alpha > 1, 1, self.alpha)))

        else:
            self.corrective_velocity_factor = Function(self.P1DG).interpolate(Constant(1.0))

        self.z_n = Function(P1)
        self.z_n1 = Function(P1)
        self.v = TestFunction(P1)
        self.old_bathymetry_2d = interpolate(self.bathymetry, P1)

        # define bed gradient
        self.dzdx = interpolate(self.old_bathymetry_2d.dx(0), P1)
        self.dzdy = interpolate(self.old_bathymetry_2d.dx(1), P1)

    def set_up_bedload(self, mesh):
        P1 = FunctionSpace(mesh, "CG", 1)

        # calculate angle of flow
        self.calfa = interpolate(self.horizontal_velocity/sqrt(self.unorm), P1)
        self.salfa = interpolate(self.vertical_velocity/sqrt(self.unorm), P1)

        self.beta = Constant(1.3)

        self.surbeta2 = Constant(1/1.5)
        self.cparam = Constant((2650-1000)*9.81*self.average_size*(self.surbeta2**2))

        if self.slope_eff:
            # slope effect magnitude correction due to gravity where beta is a parameter normally set to 1.3
            self.slopecoef = interpolate(1 + self.beta*(self.dzdx*self.calfa + self.dzdy*self.salfa), P1)
        else:
            self.slopecoef = interpolate(Constant(1.0), P1)

        if self.angle_correction:
            # slope effect angle correction due to gravity
            tt1 = conditional(1000*0.5*self.qfc*self.unorm > 10**(-10), sqrt(self.cparam/(1000*0.5*self.qfc*self.unorm)), sqrt(self.cparam/(10**(-10))))
            # add on a factor of the bed gradient to the normal
            aa = self.salfa + tt1*self.dzdy
            bb = self.calfa + tt1*self.dzdx
            norm = conditional(sqrt(aa**2 + bb**2) > 10**(-10), sqrt(aa**2 + bb**2), 10**(-10))

        # implement meyer-peter-muller bedload transport formula
        thetaprime = self.mu*(1000*0.5*self.qfc*self.unorm)/((2650-1000)*9.81*self.average_size)

        # if velocity above a certain critical value then transport occurs
        self.phi = interpolate(conditional(thetaprime < self.thetacr, 0, 8*(thetaprime-self.thetacr)**1.5), P1)

        self.z_n = Function(P1)
        self.z_n1 = Function(P1)
        self.v = TestFunction(P1)
        self.n = FacetNormal(mesh)
        self.old_bathymetry_2d = Function(P1).interpolate(self.bathymetry)

    def update_key_hydro(self, solver_obj):

        self.old_bathymetry_2d.interpolate(solver_obj.fields.bathymetry_2d)
        self.z_n.assign(self.old_bathymetry_2d)

        self.uv1, self.eta = solver_obj.fields.solution_2d.split()
        self.u_cg.project(self.uv1)
        self.elev_cg.project(self.eta)

        # calculate gradient of bed (noting bathymetry is -bed)
        self.dzdx.interpolate(self.old_bathymetry_2d.dx(0))
        self.dzdy.interpolate(self.old_bathymetry_2d.dx(1))

        self.horizontal_velocity.interpolate(self.u_cg[0])
        self.vertical_velocity.interpolate(self.u_cg[1])

        # Update depth
        if self.wetting_and_drying:
            bathymetry_displacement = solver_obj.depth.wd_bathymetry_displacement
            self.depth.interpolate(self.elev_cg + bathymetry_displacement(self.eta) + self.bathymetry)
        else:
            self.depth.interpolate(self.elev_cg + self.old_bathymetry_2d)

        self.hc = conditional(self.depth > 0.001, self.depth, 0.001)
        self.aux = conditional(11.036*self.hc/self.ks > 1.001, 11.036*self.hc/self.ks, 1.001)
        self.qfc = 2/(ln(self.aux)/0.4)**2

        if self.friction == 'nikuradse':
            self.quadratic_drag_coefficient.interpolate(self.get_cfactor())

        self.cfactor.interpolate(self.get_cfactor())

        # mu - ratio between skin friction and normal friction
        self.mu.interpolate(conditional(self.qfc > 0, self.cfactor/self.qfc, 0))

        # bed shear stress
        self.unorm = ((self.horizontal_velocity**2) + (self.vertical_velocity**2))
        self.TOB.interpolate(1000*0.5*self.qfc*self.unorm)

        self.f = (((1-self.porosity)*(self.z_n1 - self.z_n)/(self.dt*self.morfac))*self.v)*dx

    def update_suspended(self, solver_obj):

        P1DG = solver_obj.function_spaces.P1DG_2d

        self.B.interpolate(conditional(self.a > self.depth, Constant(1.0), self.a/self.depth))
        self.ustar = sqrt(0.5*self.qfc*self.unorm)
        self.exp1 = conditional((conditional((self.settling_velocity/(0.4*self.ustar)) - 1 > 0, (self.settling_velocity/(0.4*self.ustar)) - 1, -(self.settling_velocity/(0.4*self.ustar)) + 1)) > 10**(-4), conditional((self.settling_velocity/(0.4*self.ustar)) - 1 > 3, 3, (self.settling_velocity/(0.4*self.ustar))-1), 0)
        self.coefftest = conditional((conditional((self.settling_velocity/(0.4*self.ustar)) - 1 > 0, (self.settling_velocity/(0.4*self.ustar)) - 1, -(self.settling_velocity/(0.4*self.ustar)) + 1)) > 10**(-4), self.B*(1-self.B**self.exp1)/self.exp1, -self.B*ln(self.B))
        self.coeff.interpolate(conditional(self.coefftest > 0, 1/self.coefftest, 0))

        # erosion flux - van rijn
        s0 = (conditional(1000*0.5*self.qfc*self.unorm*self.mu > 0, 1000*0.5*self.qfc*self.unorm*self.mu, 0) - self.taucr)/self.taucr
        self.ceq.interpolate(0.015*(self.average_size/self.a) * ((conditional(s0 < 0, 0, s0))**(1.5))/(self.dstar**0.3))

        if self.conservative:
            self.tracer_init_value.assign(self.depth.at([0, 0])*self.ceq.at([0, 0])/self.coeff.at([0, 0]))
        else:
            self.tracer_init_value.assign(self.ceq.at([0, 0])/self.coeff.at([0, 0]))

        self.depo, self.ero = self.set_source_tracer(P1DG, solver_obj)
                

        if self.conservative:
            if self.depth_integrated:
                self.depth_int_sink.interpolate(self.depo/self.depth)
                self.depth_int_source.interpolate(self.ero)
            else:
                self.sink.interpolate(self.depo/(self.depth**2))
                self.source.interpolate(self.ero/self.depth)
            self.qbsourcedepth.interpolate(-(self.depo*solver_obj.fields.tracer_2d/self.depth) + self.ero)
        else:
            if self.implicit_source:
                self.sink.interpolate(self.depo/self.depth)
                self.source.interpolate(self.ero/self.depth)
            else:
                self.source.interpolate((-(self.depo*solver_obj.fields.tracer_2d) + self.ero)/self.depth)
            self.qbsourcedepth.interpolate(-(self.depo*solver_obj.fields.tracer_2d) + self.ero)

        if self.convective_vel_flag:
            # correction factor to advection velocity in sediment concentration equation
            self.Bconv = conditional(self.depth > 1.1*self.ksp, self.ksp/self.depth, self.ksp/(1.1*self.ksp))
            self.Aconv = conditional(self.depth > 1.1*self.a, self.a/self.depth, self.a/(1.1*self.a))

            # take max of value calculated either by ksp or depth
            self.Amax = conditional(self.Aconv > self.Bconv, self.Aconv, self.Bconv)

            self.r1conv = 1 - (1/0.4)*conditional(self.settling_velocity/self.ustar < 1, self.settling_velocity/self.ustar, 1)

            self.Ione = conditional(self.r1conv > 10**(-8), (1 - self.Amax**self.r1conv)/self.r1conv, conditional(self.r1conv < - 10**(-8), (1 - self.Amax**self.r1conv)/self.r1conv, ln(self.Amax)))

            self.Itwo = conditional(self.r1conv > 10**(-8), -(self.Ione + (ln(self.Amax)*(self.Amax**self.r1conv)))/self.r1conv, conditional(self.r1conv < - 10**(-8), -(self.Ione + (ln(self.Amax)*(self.Amax**self.r1conv)))/self.r1conv, -0.5*ln(self.Amax)**2))

            self.alpha = (-(self.Itwo - (ln(self.Amax) - ln(30))*self.Ione)/(self.Ione * ((ln(self.Amax) - ln(30)) + 1)))

            # final correction factor
            self.corrective_velocity_factor.interpolate(conditional(conditional(self.alpha > 1, 1, self.alpha) < 0, 0, conditional(self.alpha > 1, 1, self.alpha)))

        self.f += - (self.qbsourcedepth * self.v)*dx

    def update_bedload(self, solver_obj):

        # calculate angle of flow
        self.calfa.interpolate(self.horizontal_velocity/sqrt(self.unorm))
        self.salfa.interpolate(self.vertical_velocity/sqrt(self.unorm))

        if self.slope_eff:
            # slope effect magnitude correction due to gravity where beta is a parameter normally set to 1.3
            # we use z_n1 and equals so that we can use an implicit method in Exner
            self.slopecoef = (1 + self.beta*(self.z_n1.dx(0)*self.calfa + self.z_n1.dx(1)*self.salfa))
        else:
            self.slopecoef = Constant(1.0)

        if self.angle_correction:
            # slope effect angle correction due to gravity
            tt1 = conditional(1000*0.5*self.qfc*self.unorm > 10**(-10), sqrt(self.cparam/(1000*0.5*self.qfc*self.unorm)), sqrt(self.cparam/(10**(-10))))
            # add on a factor of the bed gradient to the normal
            aa = self.salfa + tt1*self.dzdy
            bb = self.calfa + tt1*self.dzdx
            norm = conditional(sqrt(aa**2 + bb**2) > 10**(-10), sqrt(aa**2 + bb**2), 10**(-10))
            # we use z_n1 and equals so that we can use an implicit method in Exner
            calfamod = (self.calfa + (tt1*self.z_n1.dx(0)))/norm
            salfamod = (self.salfa + (tt1*self.z_n1.dx(1)))/norm

        # implement meyer-peter-muller bedload transport formula
        thetaprime = self.mu*(1000*0.5*self.qfc*self.unorm)/((2650-1000)*9.81*self.average_size)

        # if velocity above a certain critical value then transport occurs
        self.phi.interpolate(conditional(thetaprime < self.thetacr, 0, 8*(thetaprime-self.thetacr)**1.5))

        # bedload transport flux with magnitude correction
        self.qb_total = self.slopecoef*self.phi*sqrt(self.g*(2650/1000 - 1)*self.average_size**3)

        # formulate bedload transport flux with correct angle depending on corrections implemented
        if self.angle_correction:
            self.qbx = self.qb_total*calfamod
            self.qby = self.qb_total*salfamod
        else:
            self.qbx = self.qb_total*self.calfa
            self.qby = self.qb_total*self.salfa

        # add bedload transport to exner equation
        self.f += -(self.v*((self.qbx*self.n[0]) + (self.qby*self.n[1])))*ds(1) - (self.v*((self.qbx*self.n[0]) + (self.qby*self.n[1])))*ds(2) + (self.qbx*(self.v.dx(0)) + self.qby*(self.v.dx(1)))*dx
