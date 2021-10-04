import plotting
import sensitivity
from params import *

# top = 0
# plotting.distribution_plots()

# plotting.plot_t23_chi_raw_profile(savename = "t23_chi_sq_profile_raw_cascade", top = 0)
# plotting.plot_t23_chi_raw_profile_all_top(savename = "t23_chi_sq_profile_raw_all_top_CORRECT")

#plotting.plot_m31_chi_raw_profile(savename = "m31_chi_sq_profile_raw_CORRECT", top = top)
# plotting.plot_m31_chi_raw_profile_all_top(savename = "m31_chi_sq_profile_raw_all_top_CORRECT")

# plotting.plot_t23_min_chi_profile(savename = "t23_chi_sq_min_profile_cascade", top = 0)
#plotting.plot_t23_min_chi_profile_all_top(savename = "t23_min_chi_sq_profile_all_top")

plotting.plot_contour_chi(savename = "chi_sq_contour_all_coarse_large_m")
# plotting.plot_contour_chi(savename = "chi_sq_contour_all_top_fine")
# plotting.plot_contour_chi(savename = "chi_sq_contour_all", top = 2)
