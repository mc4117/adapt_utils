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
