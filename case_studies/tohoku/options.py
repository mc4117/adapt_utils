from thetis import *
from thetis.configuration import *

import netCDF4
import matplotlib.pyplot as plt

from adapt_utils.swe.tsunami.options import TsunamiOptions
from adapt_utils.swe.tsunami.conversion import from_latlon


__all__ = ["TohokuOptions"]


class TohokuOptions(TsunamiOptions):
    """
    Setup for model of the Tohoku tsunami which struck the east coast of Japan in 2011, leading to
    the meltdown of Daiichi nuclear power plant, Fukushima.

    Data sources:
      * Bathymetry data extracted from GEBCO.
      * Initial free surface elevation field generated by inversion on tide gauge data by
        [Saito et al.].

    [Saito et al.] T. Saito, Y. Ito, D. Inazu, R. Hino, "Tsunami source of the 2011 Tohoku‐Oki
                   earthquake, Japan: Inversion analysis based on dispersive tsunami simulations",
                   Geophysical Research Letters (2011), 38(7).
    """
    def __init__(self, offset=475, **kwargs):
        self.force_zone_number = 54
        self.offset = offset
        super(TohokuOptions, self).__init__(**kwargs)

        # Timestepping: export once per minute for 25 minutes
        self.dt_per_export = 12
        self.dt_per_remesh = 12
        self.end_time = 1500.0
        # self.end_time = 3600.0

        # Gauge locations  # TODO: remove timeseries
        self.gauges["P02"] = {"lonlat": (142.5016, 38.5002),
                              "data": [0.00, 0.07, 0.12, 0.46, 0.85, 1.20, 1.55, 1.90, 2.25, 2.50,
                                       2.80, 3.10, 3.90, 4.80, 4.46, 2.25, -0.45, -0.17, -1.60,
                                       -0.82, -0.44, -0.26, -0.08, 0.13, 0.42, 0.71]}
        self.gauges["P06"] = {"lonlat": (142.5838, 38.6340),
                              "data": [0.00, 0.10, 0.30, 0.65, 1.05, 1.35, 1.65, 1.95, 2.25, 2.55,
                                       2.90, 3.50, 4.50, 4.85, 3.90, 1.55, -0.35, -1.05, -0.65,
                                       -0.30, -0.15, 0.05, 0.18, 0.35, 0.53, 0.74]}
        self.gauges["801"] = {"lonlat": (141.6856, 38.2325)}
        self.gauges["802"] = {"lonlat": (142.0969, 39.2586)}
        self.gauges["803"] = {"lonlat": (141.8944, 38.8578)}
        self.gauges["804"] = {"lonlat": (142.1867, 39.6272)}
        self.gauges["806"] = {"lonlat": (141.1856, 36.9714)}

        # Coastal locations of interest, including major cities and nuclear power plants
        self.locations_of_interest["Fukushima Daiichi"] = {"lonlat": (141.0281, 37.4213)}
        self.locations_of_interest["Onagawa"] = {"lonlat": (141.5008, 38.3995)}
        self.locations_of_interest["Fukushima Daini"] = {"lonlat": (141.0249, 37.3166)}
        self.locations_of_interest["Tokai"] = {"lonlat": (140.6067, 36.4664)}
        self.locations_of_interest["Hamaoka"] = {"lonlat": (138.1433, 34.6229)}
        self.locations_of_interest["Tohoku"] = {"lonlat": (141.3903, 41.1800)}
        self.locations_of_interest["Tokyo"] = {"lonlat": (139.6917, 35.6895)}

        # Convert coordinates to UTM and create timeseries array
        for loc in (self.gauges, self.locations_of_interest):
            for l in loc:
                loc[l]["timeseries"] = []
                if self.utm:
                    lon, lat = loc[l]["lonlat"]
                    loc[l]["utm"] = from_latlon(lat, lon, force_zone_number=54)
                    loc[l]["coords"] = loc[l]["utm"]
                else:
                    loc[l]["coords"] = loc[l]["lonlat"]

    def annotate_plot(self, axes, coords=None, gauges=False):
        """
        Annotate a plot on axes `axes` in coordinate system `coords` with all gauges or locations of
        interest, as determined by the Boolean kwarg `gauges`.
        """
        coords = coords or "lonlat"
        try:
            assert coords in ("lonlat", "utm")
        except AssertionError:
            raise ValueError("Coordinate system {:s} not recognised.".format(coords))
        dat = self.gauges if gauges else self.locations_of_interest
        for loc in dat:
            x, y = dat[loc][coords]
            xytext = (x + 0.3, y)
            color = "indigo"
            if loc == "P02":
                color = "navy"
                xytext = (x + 0.5, y - 0.4)
            elif loc == "P06":
                color = "navy"
                xytext = (x + 0.5, y + 0.2)
            elif "80" in loc:
                color = "darkgreen"
                xytext = (x - 0.8, y)
            elif loc == "Fukushima Daini":
                continue
            elif loc == "Fukushima Daiichi":
                loc = "Fukushima"
            elif loc in ("Tokyo", "Hamaoka"):
                xytext = (x + 0.3, y-0.6)
            ha = "center" if gauges else "left"
            axes.annotate(loc, xy=(x, y), xytext=xytext, fontsize=10, color=color, ha=ha)
            circle = plt.Circle((x, y), 0.1, color=color)
            axes.add_patch(circle)


    def read_bathymetry_file(self, km=False):
        abspath = os.path.realpath(__file__)
        fname = abspath.replace('options.py', 'resources/bathymetry.nc')
        nc = netCDF4.Dataset(fname, 'r')
        o = self.offset
        lon = nc.variables['lon'][o:]
        lat = nc.variables['lat'][:]
        rescale = 1000.0 if km else 1.0
        elev = nc.variables['elevation'][:, o:]/rescale
        nc.close()
        return lon, lat, elev

    def read_surface_file(self, zeroed=True):
        fname = 'resources/surf'
        if zeroed:
            fname = '_'.join([fname, 'zeroed'])
        fname += '.nc'
        abspath = os.path.realpath(__file__)
        fname = abspath.replace('options.py', 'resources/bathymetry.nc')
        nc = netCDF4.Dataset(fname, 'r')
        lon = nc.variables['lon' if zeroed else 'x'][:]
        lat = nc.variables['lat' if zeroed else 'y'][:]
        elev = nc.variables['z'][:, :]
        nc.close()
        return lon, lat, elev

    def set_boundary_conditions(self, fs):
        self.boundary_conditions = {}
        return self.boundary_conditions

    def set_qoi_kernel(self, solver_obj):
        pass
        # raise NotImplementedError  # TODO