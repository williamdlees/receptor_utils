from warnings import warn
from . import number_v
from . import aux_formats
from . import simple_bio_seq
from . import sequence_alignment

def __getattr__(name):
    if name == 'number_ighv':
        warn("number_ighv is deprecated. Please use number_v instead.", DeprecationWarning, stacklevel=2)
        return number_v
    raise AttributeError(f"module {__name__} has no attribute {name}")
