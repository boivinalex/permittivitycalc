import os
if os.environ.get('DISPLAY','') == '':
    print('No display found. Using non-interactive Agg backend')
    import matplotlib
    matplotlib.use('Agg')

from .sparam_data import AirlineData, multiple_meas, run_default, run_example
from . import permittivity_plot as pplot
from .helper_functions import get_METAS_data, perm_compare