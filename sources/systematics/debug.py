import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import nuSQuIDS as nsq
import nuflux
# import seaborn as sns

from params import *
import chisq as chi

print(chi.get_chisq(theta23, m31, 2))