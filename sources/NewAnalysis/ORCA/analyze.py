import numpy as np
import ORCA.digitalizer as dgt
import ORCA.MCGenerator as gen
import pandas as pd
from matplotlib import pyplot as plt
from ORCA.params import *
import ORCA.util
tracks = dgt.Digitalizer(input_track, input_scale)
cascades = dgt.Digitalizer(input_cascade, input_scale)

tracks.set_palette(0, -3, 100)
cascades.set_palette(0, -3, 100)

tracks.digitalize(22, 22)
cascades.digitalize(22, 22)

D_tracks = tracks.extracted
D_tracks[11][14] = 0.1 # only manual hardcode part
D_cascades = cascades.extracted

tracks.fit()
cascades.fit()

fakegen = gen.Generator(input_MC, tracks.gaussians, cascades.gaussians, x_bins)
res_e_true, res_e_reco, res_zen_true, res_zen_reco, res_W_mc = fakegen.generate()