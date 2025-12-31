from warnings import warn
from . import number_v

def __getattr__(name):
    if name == 'number_ighv':
        warn("number_ighv is deprecated. Please use number_v instead.", DeprecationWarning, stacklevel=2)
        return number_v


