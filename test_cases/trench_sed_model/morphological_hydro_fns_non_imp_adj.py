#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:05:02 2019

@author: mc4117
"""

import thetis_adjoint as th_adj
import time
import datetime
import numpy as np
import firedrake_adjoint as fire
import math
import os

def testing_derivative(rf, control_var, mod_var, mod):
    
    J_h = rf(control_var)
    print(J_h)
    
    der = rf.derivative()
    
    print(der)
    
    
    J_0 = rf(mod_var)
    print(J_0)


    print(der*mod)
    print(J_0-J_h)

    return J_0, J_h, der

def hydrodynamics_only(boundary_conditions_fn, mesh2d, bathymetry_2d, uv_init, elev_init, ks, average_size, dt, t_end, friction='nikuradse', friction_coef=0, fluc_bcs=False, viscosity=10**(-6), diffusivity=0.15):
    """
    Sets up a simulation with only hydrodynamics until a quasi-steady state when it can be used as an initial
    condition for the full morphological model. We update the bed friction at each time step.
    The actual run of the model are done outside the function

    Inputs:
    boundary_consditions_fn - function defining boundary conditions for problem
    mesh2d - define mesh working on
    bathymetry2d - define bathymetry of problem
    uv_init - initial velocity of problem
    elev_init - initial elevation of problem
    ks - bottom friction coefficient for quadratic drag coefficient
    average_size - average sediment size
    dt - timestep
    t_end - end time
    viscosity - viscosity of hydrodynamics. The default value should be 10**(-6)
    friction - choice of friction formulation - nikuradse and manning
    friction_coef - value of friction coefficient used in manning

    Outputs:
    solver_obj - solver which we need to run to solve the problem
    update_forcings_hydrodynamics - function defining the updates to the model performed at each timestep
    """
    def update_forcings_hydrodynamics(t_new):
        # update bed friction
        uv1, elev1 = solver_obj.fields.solution_2d.split()
        depth.project(elev1 + bathymetry_2d)
        # calculate skin friction coefficient
        cfactor.project(2*(0.4**2)/((th_adj.ln(11.036*depth/(ksp)))**2))

    # choose directory to output results
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    outputdir = 'outputs' + st

    # export interval in seconds
    if t_end < 40:
        t_export = 1
    else:
        t_export = np.round(t_end/40, 0)

    th_adj.print_output('Exporting to '+outputdir)

    # define function spaces
    V = th_adj.FunctionSpace(mesh2d, 'CG', 1)
    P1_2d = th_adj.FunctionSpace(mesh2d, 'DG', 1)

    # define parameters
    ksp = th_adj.Constant(3*average_size)

    # define depth
    depth = th_adj.Function(V).project(elev_init + bathymetry_2d)

    # define bed friction
    cfactor = th_adj.Function(P1_2d).project(2*(0.4**2)/((th_adj.ln(11.036*depth/(ksp)))**2))

    # set up solver
    solver_obj = th_adj.solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
    options = solver_obj.options
    options.simulation_export_time = t_export
    options.simulation_end_time = t_end
    options.output_directory = outputdir

    options.check_volume_conservation_2d = True
    options.fields_to_export = ['uv_2d', 'elev_2d']
    options.solve_tracer = False
    options.use_lax_friedrichs_tracer = False
    if friction == 'nikuradse':
        options.quadratic_drag_coefficient = cfactor
    elif friction == 'manning':
        if friction_coef == 0:
            friction_coef = 0.02
        options.manning_drag_coefficient = th_adj.Constant(friction_coef)
    else:
        print('Undefined friction')
    # set horizontal diffusivity parameter
    options.horizontal_diffusivity = th_adj.Constant(diffusivity)
    options.horizontal_viscosity = th_adj.Constant(viscosity)

    # crank-nicholson used to integrate in time system of ODEs resulting from application of galerkin FEM
    options.timestepper_type = 'CrankNicolson'
    options.timestepper_options.implicitness_theta = 1.0

    if not hasattr(options.timestepper_options, 'use_automatic_timestep'):
        options.timestep = dt

    # set boundary conditions

    swe_bnd, left_bnd_id, right_bnd_id, in_constant, out_constant, left_string, right_string = boundary_conditions_fn(bathymetry_2d, flag='hydro')

    for j in range(len(in_constant)):
        exec('constant_in' + str(j) + ' = th_adj.Constant(' + str(in_constant[j]) + ')', globals())

    str1 = '{'
    if len(left_string) > 0:
        for i in range(len(left_string)):
            str1 += "'" + str(left_string[i]) + "': constant_in" + str(i) + ","
        str1 = str1[0:len(str1)-1] + "}"
        exec('swe_bnd[left_bnd_id] = ' + str1)

    for k in range(len(out_constant)):
        exec('constant_out' + str(k) + '= th_adj.Constant(' + str(out_constant[k]) + ')', globals())

    str2 = '{'
    if len(right_string) > 0:
        for i in range(len(right_string)):
            str2 += "'" + str(right_string[i]) + "': constant_out" + str(i) + ","
        str2 = str2[0:len(str2)-1] + "}"
        exec('swe_bnd[right_bnd_id] = ' + str2)

    solver_obj.bnd_functions['shallow_water'] = swe_bnd

    solver_obj.assign_initial_conditions(uv=uv_init, elev=elev_init)
    return solver_obj, update_forcings_hydrodynamics


def morphological(boundary_conditions_fn, morfac, morfac_transport, suspendedload, convectivevel, bedload, angle_correction, slope_eff, seccurrent,
                  mesh2d, bathymetry_2d, input_dir, viscosity_hydro, ks_float, dt, final_time, average_size_float, diffusivity=None, 
                  beta_fn=None, surbeta2_fn=1/1.5, alpha_secc_fn=0.75, angle_fn=35, mesh_step_size=0.2, friction='nikuradse', friction_coef=0, d90=0, fluc_bcs=False, bed_form='meyer', sus_form='vanrijn', diffusivity=0.15):
    """
    Set up a full morphological model simulation using as an initial condition the results of a hydrodynamic only model.
    Inputs:
    boundary_consditions_fn - function defining boundary conditions for problem
    morfac - morphological scale factor
    morfac_transport - switch to turn on morphological component
    suspendedload - switch to turn on suspended sediment transport
    convectivevel - switch on convective velocity correction factor in sediment concentration equation
    bedload - switch to turn on bedload transport
    angle_correction - switch on slope effect angle correction
    slope_eff - switch on slope effect magnitude correction
    seccurrent - switch on secondary current for helical flow effect
    sediment_slide - switch on sediment slide mechanism to prevent steep slopes
    fluc_bcs - switch on fluctuating boundary conditions
    mesh2d - define mesh working on
    bathymetry2d - define bathymetry of problem
    input_dir - folder containing results of hydrodynamics model which are used as initial conditions here
    viscosity_hydro - viscosity value in hydrodynamic equations
    ks - bottom friction coefficient for quadratic drag coefficient
    average_size - average sediment size
    d90 - sediment size which 90% of the sediment are below
    dt - timestep
    final_time - end time
    beta_fn - magnitude slope effect parameter
    surbeta2_fn - angle correction slope effect parameter
    alpha_secc_fn - secondary current parameter
    angle_fn - maximum angle allowed by sediment slide mechanism
    mesh_step_size - NOT for defining mesh but to be used in the sediment slide mechanism
    bed_form - choice of bedload formula between 'meyer' (meyer-peter-muller) and 'soulsby' (soulsby-van-rijn)
    sus_form - choice of suspended load formula between 'vanrijn' (van Rijn 1984) and 'soulsby' (soulsby-van-rijn)
    friction - choice of friction formulation - nikuradse and manning
    friction_coef - value of friction coefficient used in manning

    Outputs:
    solver_obj - solver which we need to run to solve the problem
    update_forcings_hydrodynamics - function defining the updates to the model performed at each timestep
    diff_bathy - bedlevel evolution
    diff_bathy_file - bedlevel evolution file
    """
    if diffusivity == None:
        diffusivity_constant = th_adj.Constant(0.0)
        diffusivity_constant.assign(th_adj.AdjFloat(0.15))
    else:
        diffusivity_constant = th_adj.Constant(0.0)
        diffusivity_constant.assign(diffusivity)
        
    if beta_fn == None:
        beta = th_adj.Constant(0.0)
        beta.assign(th_adj.AdjFloat(1.3))
    else:
        beta = th_adj.Constant(0.0)
        beta.assign(beta_fn)        
        

    average_size = th_adj.Constant(0.0)
    average_size.assign(average_size_float)
    
    ks = th_adj.Constant(0.0)
    ks.assign(ks_float)    

    def update_forcings_tracer(t_new):

        # update bathymetry
        old_bathymetry_2d.project(bathymetry_2d)

        # extract new elevation and velocity and project onto CG space
        uv1, elev1 = solver_obj.fields.solution_2d.split()
        uv_cg.project(uv1)
        elev_cg.project(elev1)
        depth.project(elev_cg + old_bathymetry_2d)

        horizontal_velocity.project(uv_cg[0])
        vertical_velocity.project(uv_cg[1])

        # update boundary conditions if have fluctuating conditions
        if fluc_bcs:
            in_fn, out_fn = boundary_conditions_fn(orig_bathymetry, flag='morpho', morfac=morfac, t_new=t_new, state='update')
            for j in range(len(in_fn)):
                exec('constant_in' + str(j) + '.assign(' + str(in_fn[j]) + ')')

            for k in range(len(out_fn)):
                exec('constant_out' + str(k) + '.assign(' + str(out_fn[k]) + ')')

        # update bedfriction
        hc = th_adj.conditional(depth > 0.001, depth, 0.001)
        aux = th_adj.conditional(th_adj.ge(11.036*hc/ks_const, 1.001), 11.036*hc/ks_const, 1.001)
        qfc = 2/(th_adj.ln(aux)/0.4)**2

        # calculate skin friction coefficient
        hclip = th_adj.conditional(th_adj.ge(ksp, depth), ksp, depth)
        cfactor.project(th_adj.conditional(depth > ksp, 2*((2.5*th_adj.ln(11.036*depth/ksp))**(-2)), th_adj.Constant(0.0)))
        
        if suspendedload == False:
            t_old.assign(th_adj.AdjFloat(t_new))
        if morfac_transport == True:
            if t_old.dat.data[:] == t_new:
                z_n.project(old_bathymetry_2d)

                # mu - ratio between skin friction and normal friction
                mu.project(th_adj.conditional(qfc > 0, cfactor/qfc, th_adj.Constant(0.0)))

                # bed shear stress
                unorm = ((horizontal_velocity**2) + (vertical_velocity**2))
                TOB.project(1000*0.5*qfc*unorm)

                # calculate gradient of bed (noting bathymetry is -bed)
                dzdx.project(old_bathymetry_2d.dx(0))
                dzdy.project(old_bathymetry_2d.dx(1))

                # initialise exner equation if not already initialised in sediment slide
                f = 0

                if suspendedload:
                    # source term

                    # deposition flux - calculating coefficient to account for stronger conc at bed
                    B.project(th_adj.conditional(a > depth, th_adj.Constant(1.0), a/depth))
                    ustar = (th_adj.sqrt(0.5*qfc*unorm))
                    exp1 = (th_adj.conditional((th_adj.conditional((settling_velocity/(0.4*ustar)) - 1 > 0, (settling_velocity/(0.4*ustar)) - 1, -(settling_velocity/(0.4*ustar)) + 1)) > 10**(-4), th_adj.conditional((settling_velocity/(0.4*ustar)) - 1 > 3, 3, (settling_velocity/(0.4*ustar))-1), 0))
                    coefftest = (th_adj.conditional((th_adj.conditional((settling_velocity/(0.4*ustar)) - 1 > 0, (settling_velocity/(0.4*ustar)) - 1, -(settling_velocity/(0.4*ustar)) + 1)) > 10**(-4), B*(1-B**exp1)/exp1, -B*th_adj.ln(B)))
                    coeff.project(th_adj.conditional(coefftest > 0, 1/coefftest, 0))

                    if sus_form == 'vanrijn':
                        # erosion flux - above critical velocity bed is eroded
                        s0 = (th_adj.conditional(1000*0.5*qfc*unorm*mu > 0, 1000*0.5*qfc*unorm*mu, 0) - taucr)/taucr
                        ceq.project(0.015*(average_size/a) * ((th_adj.conditional(s0 < 0, 0, s0))**(1.5))/(dstar**0.3))
                    elif sus_form == 'soulsby':
                        ucr.assign(0.19*(average_size**0.1)*(th_adj.ln(4*depth/d90)/th_adj.ln(10)))
                        s0.assign(th_adj.conditional((th_adj.sqrt(unorm)-ucr)**2.4 > 0, (th_adj.sqrt(unorm)-ucr)**2.4, 0))
                        ceq.assign(ass*s0/depth)
                    else:
                        print('Unrecognised suspended sediment transport formula. Please choose "vanrijn" or "soulsby"')

                    # calculate depth-averaged source term for sediment concentration equation
                    depo = (settling_velocity * coeff)
                    ero = (settling_velocity * ceq)
                    source.project((-(depo*solver_obj.fields.tracer_2d) + ero)/depth)
                    # update sediment rate to ensure equilibrium at inflow
                    sediment_rate.project(ceq/coeff)
                    
                    if convectivevel:
                        # correction factor to advection velocity in sediment concentration equation
                        Bconv = (th_adj.conditional(depth > 1.1*ksp, ksp/depth, th_adj.Constant(1/1.1)))
                        Aconv = (th_adj.conditional(depth > 1.1*a, a/depth, th_adj.Constant(1/1.1)))

                        # take max of value calculated either by ksp or depth
                        Amax = (th_adj.conditional(Aconv > Bconv, Aconv, Bconv))

                        r1conv = (1 - (1/0.4)*th_adj.conditional(settling_velocity/ustar < 1, settling_velocity/ustar, 1))

                        Ione = (th_adj.conditional(r1conv > 10**(-8), (1 - Amax**r1conv)/r1conv, th_adj.conditional(r1conv < - 10**(-8), (1 - Amax**r1conv)/r1conv, th_adj.ln(Amax))))

                        Itwo = (th_adj.conditional(r1conv > 10**(-8), -(Ione + (th_adj.ln(Amax)*(Amax**r1conv)))/r1conv, th_adj.conditional(r1conv < - 10**(-8), -(Ione + (th_adj.ln(Amax)*(Amax**r1conv)))/r1conv, -0.5*th_adj.ln(Amax)**2)))

                        alpha = (-(Itwo - (th_adj.ln(Amax) - th_adj.ln(30))*Ione)/(Ione * ((th_adj.ln(Amax) - th_adj.ln(30)) + 1)))

                        # final correction factor
                        alphatest2.project(th_adj.conditional(th_adj.conditional(alpha > 1, 1, alpha) < 0, 0, th_adj.conditional(alpha > 1, 1, alpha)))

                if bedload:

                    # calculate angle of flow
                    calfa.project(horizontal_velocity/th_adj.sqrt(unorm))
                    salfa.project(vertical_velocity/th_adj.sqrt(unorm))

                    if slope_eff:
                        # slope effect magnitude correction due to gravity where beta is a parameter normally set to 1.3
                        # we use z_n1 and equals so that we can use an implicit method in Exner
                        slopecoef = (1 + beta*(z_n1.dx(0)*calfa + z_n1.dx(1)*salfa))
                    else:
                        slopecoef = th_adj.Constant(1.0)

                    if angle_correction:
                        # slope effect angle correction due to gravity
                        tt1 = (th_adj.conditional(1000*0.5*qfc*unorm > 10**(-10), th_adj.sqrt(cparam/(1000*0.5*qfc*unorm)), th_adj.sqrt(cparam/(10**(-10)))))
                        # add on a factor of the bed gradient to the normal
                        aa = (salfa + tt1*dzdy)
                        bb = (calfa + tt1*dzdx)
                        norm = (th_adj.conditional(th_adj.sqrt(aa**2 + bb**2) > 10**(-10), th_adj.sqrt(aa**2 + bb**2), 10**(-10)))
                        # we use z_n1 and equals so that we can use an implicit method in Exner
                        calfamod = (calfa + (tt1*z_n1.dx(0)))/norm
                        salfamod = (salfa + (tt1*z_n1.dx(1)))/norm

                    if seccurrent:
                        # accounts for helical flow effect in a curver channel

                        # again use z_n1 and equals so can use an implicit method in Exner
                        free_surface_dx = depth.dx(0) - z_n1.dx(0)
                        free_surface_dy = depth.dx(1) - z_n1.dx(1)

                        velocity_slide = (horizontal_velocity*free_surface_dy)-(vertical_velocity*free_surface_dx)

                        tandelta_factor.project(7*9.81*1000*depth*qfc/(2*alpha_secc*((horizontal_velocity**2) + (vertical_velocity**2))))

                        if angle_correction:
                            # if angle has already been corrected we must alter the corrected angle to obtain the corrected secondary current angle
                            t_1 = (TOB*slopecoef*calfamod) + (vertical_velocity*tandelta_factor*velocity_slide)
                            t_2 = (TOB*slopecoef*salfamod) - (horizontal_velocity*tandelta_factor*velocity_slide)
                        else:
                            t_1 = (TOB*slopecoef*calfa) + (vertical_velocity*tandelta_factor*velocity_slide)
                            t_2 = ((TOB*slopecoef*salfa) - (horizontal_velocity*tandelta_factor*velocity_slide))

                        # calculated to normalise the new angles
                        t4 = th_adj.sqrt((t_1**2) + (t_2**2))

                        # updated magnitude correction and angle corrections
                        slopecoef = t4/TOB

                        calfanew = t_1/t4
                        salfanew = t_2/t4

                    if bed_form == 'meyer':
                        # implement meyer-peter-muller bedload transport formula
                        thetaprime = (mu*(1000*0.5*qfc*unorm)/((2650-1000)*9.81*average_size))

                        # if velocity above a certain critical value then transport occurs
                        phi.project(th_adj.conditional(thetaprime < thetacr, 0, 8*(thetaprime-thetacr)**1.5))

                        # bedload transport flux with magnitude correction
                        qb_total = slopecoef*phi*th_adj.sqrt(g*(2650/1000 - 1)*average_size**3)
                    elif bed_form == 'soulsby':
                        abb.project(th_adj.conditional(depth >= average_size, 0.005*depth*((average_size/depth)**1.2)/coeff_soulsby, 0.005*depth/coeff_soulsby))
                        ucr_bed.project(th_adj.conditional(depth > d90, 0.19*(average_size**0.1)*(th_adj.ln(4*depth/d90))/(th_adj.ln(10)), 0.19*(average_size**0.1)*(th_adj.ln(4))/(th_adj.ln(10))))
                        s0_bed.project(th_adj.conditional((th_adj.sqrt(unorm)-ucr_bed)**2.4 > 0, (th_adj.sqrt(unorm)-ucr_bed)**2.4, 0))
                        qb_total = slopecoef*abb*s0_bed*th_adj.sqrt(unorm)
                    else:
                        print('Unrecognised bedload transport formula. Please choose "meyer" or "soulsby"')
                    # add time derivative to exner equation with a morphological scale factor
                    f += (((1-porosity)*(z_n1 - z_n)/(dt*morfac)) * v)*fire.firedrake.dx

                    # formulate bedload transport flux with correct angle depending on corrections implemented
                    if angle_correction and not seccurrent:
                        qbx = qb_total*calfamod
                        qby = qb_total*salfamod
                    elif seccurrent:
                        qbx = qb_total*calfanew
                        qby = qb_total*salfanew
                    else:
                        qbx = qb_total*calfa
                        qby = qb_total*salfa

                    # add bedload transport to exner equation
                    f += -(v*((qbx*n[0]) + (qby*n[1])))*fire.firedrake.ds(1) - (v*((qbx*n[0]) + (qby*n[1])))*fire.firedrake.ds(2) + (qbx*(v.dx(0)) + qby*(v.dx(1)))*fire.firedrake.dx

                else:
                    # if no bedload transport component initialise exner equation with time derivative
                    f = (((1-porosity)*(z_n1 - z_n)/(dt*morfac)) * v)*fire.firedrake.dx

                if suspendedload:
                    # add suspended sediment transport to exner equation multiplied by depth as the exner equation is not depth-averaged

                    qbsourcedepth.project(-(depo*solver_obj.fields.tracer_2d) + ero)
                    f += - (qbsourcedepth*v)*fire.firedrake.dx

                # solve exner equation using finite element methods
                fire.solve(f == 0, z_n1)

                # update bed
                bathymetry_2d.assign(z_n1)

                if round(t_new, 2) % t_export == 0:
                    # calculate difference between original bathymetry and new bathymetry
                    bathymetry_file.write(bathymetry_2d)
                    diff_bathy.project(-bathymetry_2d + orig_bathymetry)
                    diff_bathy_file.write(diff_bathy)
            t_old.assign(th_adj.AdjFloat(t_new))

    # choose directory to output results
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    outputdir = 'outputs' + st

    # final time of simulation
    t_end = final_time/morfac
    print(final_time)
    print(t_end)
    print(morfac)
    # export interval in seconds
    if t_end > 100:
        t_export = np.round(t_end/100, 0)
    else:
        t_export = 1
        
    t_old = th_adj.Constant(0.0)

    th_adj.print_output('Exporting to '+outputdir)
    x, y = th_adj.SpatialCoordinate(mesh2d)

    # define function spaces
    P1_2d = th_adj.FunctionSpace(mesh2d, 'DG', 1)
    V = th_adj.FunctionSpace(mesh2d, 'CG', 1)
    vector_cg = th_adj.VectorFunctionSpace(mesh2d, 'CG', 1)

    # define test functions on mesh
    v = fire.TestFunction(V)
    n = th_adj.FacetNormal(mesh2d)
    z_n1 = fire.Function(V, name="z^{n+1}")
    z_n = fire.Function(V, name="z^{n}")

    # define original bathymetry before bedlevel changes
    orig_bathymetry = th_adj.Function(V).project(bathymetry_2d)

    # calculate bed evolution
    diff_bathy = th_adj.Function(V).project(-bathymetry_2d + orig_bathymetry)

    # define output file for bed evolution
    diff_bathy_file = th_adj.File(outputdir + "/diff_bathy.pvd")
    diff_bathy_file.write(diff_bathy)

    # define parameters
    g = th_adj.Constant(9.81)
    porosity = th_adj.Constant(0.4
                               
    ks_const = th_adj.Function(V).assign(ks)
    ksp = th_adj.Function(V).project(3*average_size)
    a = th_adj.Function(V).project(ks_const/2)
    viscosity = th_adj.Constant(10**(-6))

    # angle correction slope effect parameters
    surbeta2 = th_adj.Constant(surbeta2_fn)
    cparam = th_adj.Constant((2650-1000)*9.81*average_size*(surbeta2**2))
    # secondary current parameter
    alpha_secc = th_adj.Constant(alpha_secc_fn)

    # calculate critical shields parameter thetacr
    R = th_adj.Constant(2650/1000 - 1)
    dstar = th_adj.Function(P1_2d).project(average_size*((g*R)/(viscosity**2))**(1/3))
    if min(dstar.dat.data[:] < 1):
        raise ValueError('ERROR: dstar value less than 1')

    thetacr = th_adj.Function(P1_2d).project(th_adj.conditional(th_adj.le(dstar, 4), 0.24*(dstar**(-1)), \
                   th_adj.conditional(th_adj.le(dstar, 10), 0.14*(dstar**(-0.64)),\
                   th_adj.conditional(th_adj.le(dstar, 20), 0.04*(dstar**(-0.1)), \
                   th_adj.conditional(th_adj.le(dstar, 150), 0.013*(dstar**(0.29)),\
                   th_adj.Constant(0.055))))))

    # critical bed shear stress
    taucr = th_adj.Function(P1_2d).project((2650-1000)*g*average_size*thetacr)
    
    # calculate settling velocity
    average_size_func = th_adj.Function(P1_2d).assign(average_size)
    #settling_velocity = th_adj.Function(P1_2d).project((10*viscosity/average_size_func)*(th_adj.sqrt(1 + 0.01*((R*g*(average_size_func**3))/(viscosity**2)))-1))

    settling_velocity = th_adj.Function(P1_2d).project(th_adj.conditional(average_size_func < 100*(10**(-6)), g*(average_size_func**2)*R/(18*viscosity),\
                                                   th_adj.conditional(th_adj.le(average_size_func, 1000*(10**(-6))), (10*viscosity/average_size_func)*(th_adj.sqrt(1 + 0.01*((R*g*(average_size_func**3))/(viscosity**2)))-1),\
                                                   1.1*th_adj.sqrt(g*average_size_func*R))))

    # initialise velocity, elevation and depth
    elev_init, uv_init = initialise_fields(mesh2d, input_dir, outputdir)

    uv_cg = th_adj.Function(vector_cg).project(uv_init)

    elev_cg = th_adj.Function(V).project(elev_init)

    depth = th_adj.Function(V).project(elev_cg + bathymetry_2d)

    old_bathymetry_2d = th_adj.Function(V).project(bathymetry_2d)

    horizontal_velocity = th_adj.Function(V).project(uv_cg[0])
    vertical_velocity = th_adj.Function(V).project(uv_cg[1])

    # define bed friction
    hc = th_adj.conditional(depth > 0.001, depth, 0.001)
    aux = th_adj.conditional(11.036*hc/ks_const > 1.001, 11.036*hc/ks_const, 1.001)
    qfc = 2/(th_adj.ln(aux)/0.4)**2
    # skin friction coefficient
    hclip = th_adj.conditional(ksp > depth, ksp, depth)
    cfactor = th_adj.Function(P1_2d).project(th_adj.conditional(depth > ksp, 2*((2.5*th_adj.ln(11.036*hclip/ksp))**(-2)), th_adj.Constant(0.0)))
    # mu - ratio between skin friction and normal friction
    mu = th_adj.Function(P1_2d).project(th_adj.conditional(qfc > 0, cfactor/qfc, 0))

    # calculate bed shear stress
    unorm = (horizontal_velocity**2) + (vertical_velocity**2)
    TOB = th_adj.Function(V).project(1000*0.5*qfc*unorm)

    # define bed gradient
    dzdx = th_adj.Function(V).project(old_bathymetry_2d.dx(0))
    dzdy = th_adj.Function(V).project(old_bathymetry_2d.dx(1))

    if suspendedload:
        # deposition flux - calculating coefficient to account for stronger conc at bed
        B = th_adj.Function(P1_2d).project(th_adj.conditional(a > depth, th_adj.Constant(1.0), a/depth))
        ustar = th_adj.sqrt(0.5*qfc*unorm)
        exp1 = th_adj.conditional((th_adj.conditional((settling_velocity/(0.4*ustar)) - 1 > 0, (settling_velocity/(0.4*ustar)) - 1, -(settling_velocity/(0.4*ustar)) + 1)) > 10**(-4), th_adj.conditional((settling_velocity/(0.4*ustar)) - 1 > 3, 3, (settling_velocity/(0.4*ustar))-1), 0)
        coefftest = th_adj.conditional((th_adj.conditional((settling_velocity/(0.4*ustar)) - 1 > 0, (settling_velocity/(0.4*ustar)) - 1, -(settling_velocity/(0.4*ustar)) + 1)) > 10**(-4), B*(1-B**exp1)/exp1, -B*th_adj.ln(B))
        coeff = th_adj.Function(P1_2d).project(th_adj.conditional(coefftest > 0, 1/coefftest, 0))

        if sus_form == 'vanrijn':
            # erosion flux - above critical velocity bed is eroded
            s0 = (th_adj.conditional(1000*0.5*qfc*unorm*mu > 0, 1000*0.5*qfc*unorm*mu, 0) - taucr)/taucr
            ceq = th_adj.Function(P1_2d).project(0.015*(average_size/a) * ((th_adj.conditional(s0 < 0, 0, s0))**(1.5))/(dstar**0.3))
        elif sus_form == 'soulsby':
            if d90 == 0:
                # if the value of d90 is unspecified set d90 = d50
                d90 = th_adj.Constant(average_size)
            else:
                d90 = th_adj.Constant(d90)
            coeff_soulsby = th_adj.Constant((R*g*average_size)**1.2)
            ass = th_adj.Constant(0.012*average_size*(dstar**(-0.6))/coeff_soulsby)
            ucr = th_adj.Function(P1_2d).project(0.19*(average_size**0.1)*(th_adj.ln(4*depth/d90))/(th_adj.ln(10)))
            s0 = th_adj.Function(P1_2d).project(th_adj.conditional((th_adj.sqrt(unorm)-ucr)**2.4 > 0, (th_adj.sqrt(unorm)-ucr)**2.4, 0))
            ceq = th_adj.Function(P1_2d).project(ass*s0/depth)
        else:
            print('Unrecognised suspended sediment transport formula. Please choose "vanrijn" or "soulsby"')

        # update sediment rate to ensure equilibrium at inflow

        sediment_rate = th_adj.Function(P1_2d).project(ceq/coeff)
        testtracer = th_adj.Function(P1_2d).project(ceq/coeff)

        # calculate depth-averaged source term for sediment concentration equation
        depo = settling_velocity * coeff
        ero = settling_velocity * ceq   
        source = th_adj.Function(P1_2d).projecct((-(depo*testtracer) + ero)/depth)

        # add suspended sediment transport to exner equation multiplied by depth as the exner equation is not depth-averaged
        qbsourcedepth = th_adj.Function(P1_2d).project(-(depo*testtracer) + ero)

        if convectivevel:
            # correction factor to advection velocity in sediment concentration equation

            Bconv = th_adj.conditional(depth > 1.1*ksp, ksp/depth, ksp/(1.1*ksp))
            Aconv = th_adj.conditional(depth > 1.1*a, a/depth, a/(1.1*a))

            # take max of value calculated either by ksp or depth
            Amax = th_adj.conditional(Aconv > Bconv, Aconv, Bconv)

            r1conv = 1 - (1/0.4)*th_adj.conditional(settling_velocity/ustar < 1, settling_velocity/ustar, 1)

            Ione = th_adj.conditional(r1conv > 10**(-8), (1 - Amax**r1conv)/r1conv, th_adj.conditional(r1conv < - 10**(-8), (1 - Amax**r1conv)/r1conv, th_adj.ln(Amax)))

            Itwo = th_adj.conditional(r1conv > 10**(-8), -(Ione + (th_adj.ln(Amax)*(Amax**r1conv)))/r1conv, th_adj.conditional(r1conv < - 10**(-8), -(Ione + (th_adj.ln(Amax)*(Amax**r1conv)))/r1conv, -0.5*th_adj.ln(Amax)**2))

            alpha = -(Itwo - (th_adj.ln(Amax) - th_adj.ln(30))*Ione)/(Ione * ((th_adj.ln(Amax) - th_adj.ln(30)) + 1))

            # final correction factor
            alphatest2 = th_adj.Function(P1_2d).project(th_adj.conditional(th_adj.conditional(alpha > 1, 1, alpha) < 0, 0, th_adj.conditional(alpha > 1, 1, alpha)))
        else:
            alphatest2 = th_adj.Function(P1_2d).project(th_adj.Constant(1.0))

    if bedload:
        # calculate angle of flow
        calfa = th_adj.Function(V).project(horizontal_velocity/th_adj.sqrt(unorm))
        salfa = th_adj.Function(V).project(vertical_velocity/th_adj.sqrt(unorm))

        if slope_eff:
            # slope effect magnitude correction due to gravity where beta is a parameter normally set to 1.3
            slopecoef = 1 + beta*(dzdx*calfa + dzdy*salfa)
        else:
            slopecoef = th_adj.Constant(1.0)

        if angle_correction:
            # slope effect angle correction due to gravity
            tt1 = th_adj.conditional(1000*0.5*qfc*unorm > 10**(-10), th_adj.sqrt(cparam/(1000*0.5*qfc*unorm)), th_adj.sqrt(cparam/(10**(-10))))
            # add on a factor of the bed gradient to the normal
            aa = salfa + tt1*dzdy
            bb = calfa + tt1*dzdx
            norm = th_adj.conditional(th_adj.sqrt(aa**2 + bb**2) > 10**(-10), th_adj.sqrt(aa**2 + bb**2), 10**(-10))

        if seccurrent:
            # accounts for helical flow effect in a curver channel
            free_surface_dx = th_adj.Function(V).project(elev_cg.dx(0))
            free_surface_dy = th_adj.Function(V).project(elev_cg.dx(1))

            velocity_slide = (horizontal_velocity*free_surface_dy)-(vertical_velocity*free_surface_dx)

            tandelta_factor = th_adj.Function(V).projet(7*9.81*1000*depth*qfc/(2*alpha_secc*((horizontal_velocity**2) + (vertical_velocity**2))))

            t_1 = (TOB*slopecoef*calfa) + (vertical_velocity*tandelta_factor*velocity_slide)
            t_2 = ((TOB*slopecoef*salfa) - (horizontal_velocity*tandelta_factor*velocity_slide))

            # calculated to normalise the new angles
            t4 = th_adj.sqrt((t_1**2) + (t_2**2))

            # updated magnitude correction and angle corrections
            slopecoef = t4/TOB

        if bed_form == 'meyer':
            # implement meyer-peter-muller bedload transport formula
            thetaprime = mu*(1000*0.5*qfc*unorm)/((2650-1000)*9.81*average_size)

            # if velocity above a certain critical value then transport occurs
            phi = th_adj.Function(V).project(th_adj.conditional(thetaprime < thetacr, 0, 8*(thetaprime-thetacr)**1.5))

        elif bed_form == 'soulsby':
            if d90 == 0:
                d90 = th_adj.Constant(average_size)
            coeff_soulsby = th_adj.Constant((R*g*average_size)**1.2)
            abb = th_adj.Function(P1_2d).project(th_adj.conditional(depth >= average_size, 0.005*depth*((average_size/depth)**1.2)/coeff_soulsby, 0.005*depth/coeff_soulsby))
            ucr_bed = th_adj.Function(P1_2d).project(th_adj.conditional(depth > d90, 0.19*(average_size**0.1)*(th_adj.ln(4*depth/d90))/(th_adj.ln(10)), 0.19*(average_size**0.1)*(th_adj.ln(4))/(th_adj.ln(10))))
            s0_bed = th_adj.Function(P1_2d).project(th_adj.conditional((th_adj.sqrt(unorm)-ucr_bed)**2.4 > 0, (th_adj.sqrt(unorm)-ucr_bed)**2.4, 0))
        else:
            print('Unrecognised bedload transport formula. Please choose "meyer" or "soulsby"')
    # set up solver
    solver_obj = th_adj.solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
    options = solver_obj.options
    options.simulation_export_time = t_export
    options.simulation_end_time = t_end
    options.output_directory = outputdir
    options.check_volume_conservation_2d = True
    if suspendedload:
        # switch on tracer calculation if using sediment transport component
        options.solve_tracer = True
        options.fields_to_export = ['uv_2d', 'elev_2d', 'tracer_2d', 'bathymetry_2d']
        options.tracer_advective_velocity_factor = alphatest2
        options.tracer_source_2d = source
    else:
        options.solve_tracer = False
        options.fields_to_export = ['uv_2d', 'elev_2d', 'bathymetry_2d']
    options.use_lax_friedrichs_tracer = False
    # set bed friction
    if friction == 'nikuradse':
        options.quadratic_drag_coefficient = cfactor
    elif friction == 'manning':
        if friction_coef == 0:
            friction_coef = 0.02
        options.manning_drag_coefficient = th_adj.Constant(friction_coef)
    else:
        print('Undefined friction')
    # set horizontal diffusivity parameter
    options.horizontal_diffusivity = th_adj.Constant(diffusivity)
    options.horizontal_viscosity = th_adj.Constant(viscosity_hydro)
    # crank-nicholson used to integrate in time system of ODEs resulting from application of galerkin FEM
    options.timestepper_type = 'CrankNicolson'
    options.timestepper_options.implicitness_theta = 1.0

    bathymetry_file = th_adj.File(outputdir + "/bathy.pvd")
    bathymetry_file.write(solver_obj.fields.bathymetry_2d)

    if not hasattr(options.timestepper_options, 'use_automatic_timestep'):
        options.timestep = dt

    
    # set boundary conditions
    swe_bnd, left_bnd_id, right_bnd_id, in_constant, out_constant, left_string, right_string = boundary_conditions_fn(orig_bathymetry, flag='morpho')

    for j in range(len(in_constant)):
        exec('constant_in' + str(j) + ' = th_adj.Constant(' + str(in_constant[j]) + ')', globals())

    str1 = '{'
    if len(left_string) > 0:
        for i in range(len(left_string)):
            str1 += "'" + str(left_string[i]) + "': constant_in" + str(i) + ","
        str1 = str1[0:len(str1)-1] + "}"
        exec('swe_bnd[left_bnd_id] = ' + str1)

    for k in range(len(out_constant)):
        exec('constant_out' + str(k) + '= th_adj.Constant(' + str(out_constant[k]) + ')', globals())

    str2 = '{'
    if len(right_string) > 0:
        for i in range(len(right_string)):
            str2 += "'" + str(right_string[i]) + "': constant_out" + str(i) + ","
        str2 = str2[0:len(str2)-1] + "}"
        exec('swe_bnd[right_bnd_id] = ' + str2)
        
    solver_obj.bnd_functions['shallow_water'] = swe_bnd
    
    if suspendedload:
        solver_obj.bnd_functions['tracer'] = {left_bnd_id: {'value': sediment_rate}}
        
        for i in solver_obj.bnd_functions['tracer'].keys():
            if i in solver_obj.bnd_functions['shallow_water'].keys():
                solver_obj.bnd_functions['tracer'][i].update(solver_obj.bnd_functions['shallow_water'][i])
        for i in solver_obj.bnd_functions['shallow_water'].keys():
            if i not in solver_obj.bnd_functions['tracer'].keys():
                solver_obj.bnd_functions['tracer'].update({i:solver_obj.bnd_functions['shallow_water'][i]})

        # set initial conditions
        solver_obj.assign_initial_conditions(uv=uv_init, elev=elev_init, tracer=testtracer)

    else:
        # set initial conditions
        solver_obj.assign_initial_conditions(uv=uv_init, elev=elev_init)

    solver_obj.iterate(update_forcings = update_forcings_tracer)
    return bathymetry_2d


def export_final_state(inputdir, uv, elev,):
    """
    Export fields to be used in a subsequent simulation
    """
    if not os.path.exists(inputdir):
        os.makedirs(inputdir)
    th_adj.print_output("Exporting fields for subsequent simulation")
    chk = th_adj.DumbCheckpoint(inputdir + "/velocity", mode=th_adj.FILE_CREATE)
    chk.store(uv, name="velocity")
    th_adj.File(inputdir + '/velocityout.pvd').write(uv)
    chk.close()
    chk = th_adj.DumbCheckpoint(inputdir + "/elevation", mode=th_adj.FILE_CREATE)
    chk.store(elev, name="elevation")
    th_adj.File(inputdir + '/elevationout.pvd').write(elev)
    chk.close()

def initialise_fields(mesh2d, inputdir, outputdir,):
    """
    Initialise simulation with results from a previous simulation
    """
    DG_2d = th_adj.FunctionSpace(mesh2d, 'DG', 1)
    # elevation
    with th_adj.timed_stage('initialising elevation'):
        chk = th_adj.DumbCheckpoint(inputdir + "/elevation", mode=th_adj.FILE_READ)
        elev_init = th_adj.Function(DG_2d, name="elevation")
        chk.load(elev_init)
        th_adj.File(outputdir + "/elevation_imported.pvd").write(elev_init)
        chk.close()
    # velocity
    with th_adj.timed_stage('initialising velocity'):
        chk = th_adj.DumbCheckpoint(inputdir + "/velocity", mode=th_adj.FILE_READ)
        V = th_adj.VectorFunctionSpace(mesh2d, 'DG', 1)
        uv_init = th_adj.Function(V, name="velocity")
        chk.load(uv_init)
        th_adj.File(outputdir + "/velocity_imported.pvd").write(uv_init)
        chk.close()
        return elev_init, uv_init,
