import plotting
import sensitivity
from params import *


# plotting.distribution_plots()
plotting.plot_t23_chi_raw_profile("t23_raw_profile_new_cascade", 0)
plotting.plot_t23_chi_raw_profile("t23_raw_profile_new_track", 1)
plotting.plot_t23_chi_raw_profile("t23_raw_profile_new_all", 2)
# plotting.plot_m31_chi_raw_profile()
# plotting.plot_contour_chi()
# plotting.plot_t23_chi_raw_profile_all_top(savename = "t23_chi_sq_profile_raw_all_top")
# plotting.plot_resolution()

# For debug purposes
# print("chi-sq is", sensitivity.get_chisq(theta23, m31, top = 1))
# sensitivity.get_rated_weight_truth(top = 1)
# sensitivity.get_energy_bins(theta23, m31, top = 1)
# print(sensitivity.get_t23_chi_profile(top = 2))

