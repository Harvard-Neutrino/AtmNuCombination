import numpy as np
import digitalizer as dgt
import MCGenerator as gen
import pandas as pd
from matplotlib import pyplot as plt
from params import *
import util

tracks = dgt.Digitalizer(input_track, input_scale)
cascades = dgt.Digitalizer(input_cascade, input_scale)

tracks.set_palette(0, -3, 20)
cascades.set_palette(0, -3, 20)

tracks.digitalize(22, 22)
cascades.digitalize(22, 22)

D_tracks = tracks.extracted
D_cascades = cascades.extracted

tracks.fit()
cascades.fit()

tracksgen = gen.Generator(input_MC, tracks.gaussians, cascades.gaussians, x_bins)
res_true, res_reco = tracksgen.generate(1, 1)
