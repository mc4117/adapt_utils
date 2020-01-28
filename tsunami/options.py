from thetis import *

import scipy.interpolate as si
from scipy.io.netcdf import NetCDFFile
import numpy as np

from adapt_utils.swe.options import ShallowWaterOptions
from adapt_utils.tsunami.conversion import lonlat_to_utm, to_latlon


__all__ = ["TohokuOptions"]


class TohokuOptions(ShallowWaterOptions):
    # TODO: doc
    def __init__(self, approach='fixed_mesh', utm=True):
        super(TohokuOptions, self).__init__(approach=approach)
        self.utm = utm
        lon, lat, elev = self.read_bathymetry_file()
        lon_min = np.min(lon)
        lon_max = np.max(lon)
        lat_min = np.min(lat)
        lat_max = np.max(lat)
        self.default_mesh = RectangleMesh(40, 40, lon_max-lon_min, lat_max-lat_min)
        x, y = SpatialCoordinate(self.default_mesh)
        self.default_mesh.coordinates.interpolate(as_vector([x+lon_min, y+lat_min]))
        if utm:
            self.default_mesh.coordinates.interpolate(as_vector(lonlat_to_utm(y, x, force_zone_number=54)))
        P1 = FunctionSpace(self.default_mesh, "CG", 1)
        self.set_bathymetry(P1, dat=(lon, lat, elev))
        self.set_bathymetry(P1)
        self.set_initial_condition(P1)

    def read_bathymetry_file(self, km=False):
        nc = NetCDFFile('tohoku.nc', mmap=False)
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:-1]
        elev = nc.variables['elevation'][:-1,:]/1000 if km else nc.variables['elevation'][:-1,:]
        nc.close()
        return lon, lat, elev

    def read_surface_file(self):
        nc = NetCDFFile('surf_zeroed.nc', mmap=False)
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        elev = nc.variables['z'][:,:]
        nc.close()
        return lon, lat, elev

    def set_bathymetry(self, fs, dat=None):
        """Initial surface data courtesy of GEBCO."""
        assert fs.ufl_element().degree() == 1 and fs.ufl_element().family() == 'Lagrange'
        x0, y0, elev = dat or self.read_bathymetry_file()
        if self.utm:
            x0, y0 = lonlat_to_utm(y0, x0, force_zone_number=54)
        bathy_interp = si.RectBivariateSpline(y0, x0, elev)
        self.bathymetry = Function(fs, name="Bathymetry")
        self.print_debug("Interpolating bathymetry...")
        msg = "Coordinates ({:.1f}, {:.1f}) Bathymetry {:.3f} km"
        for i in range(fs.mesh().num_vertices()):
            xy = fs.mesh().coordinates.dat.data[i] 
            self.bathymetry.dat.data[i] = bathy_interp(xy[1], xy[0])
            self.print_debug(msg.format(xy[0], xy[1], self.bathymetry.dat.data[i]/1000))
        self.print_debug("Done!")
        return self.bathymetry

    def set_initial_condition(self, fs):
        """Initial suface data courtesy of Saito."""
        assert fs.ufl_element().degree() == 1 and fs.ufl_element().family() == 'Lagrange'
        x0, y0, elev = self.read_surface_file()
        if self.utm:
            x0, y0 = lonlat_to_utm(y0, x0, force_zone_number=54)
        surf_interp = si.RectBivariateSpline(y0, x0, elev)
        self.initial_value = Function(fs, name="Initial free surface")
        self.print_debug("Interpolating initial surface...")
        msg = "Coordinates ({:.1f}, {:.1f}) Surface {:.3f} m"
        for i in range(fs.mesh().num_vertices()):
            xy = fs.mesh().coordinates.dat.data[i] 
            self.initial_value.dat.data[i] = surf_interp(xy[1], xy[0])
            self.print_debug(msg.format(xy[0], xy[1], self.initial_value.dat.data[i]))
        self.print_debug("Done!")
        return self.initial_value
