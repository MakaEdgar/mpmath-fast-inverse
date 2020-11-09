from .mpinv_wrapper import *
from .mpmath_wrapper import *

print("WARNING: current mpinv version works only with 50 mpmath digits. Include \"mpmath.mp.dps = 50\" to your code")
import mpmath as mp; mp.mp.dps = 50
