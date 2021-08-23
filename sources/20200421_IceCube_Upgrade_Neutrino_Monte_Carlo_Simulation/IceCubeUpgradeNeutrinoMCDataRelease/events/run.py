import plotting
import sensitivity
from params import *


# plotting.distribution_plots()
#print("this is output file for both plot")
#plotting.plot_t23_chi_raw_profile("t23_raw_profile_new_both", 2)
# plotting.plot_t23_chi_raw_profile("t23_raw_profile_new_track", 1)
# plotting.plot_t23_chi_raw_profile("t23_raw_profile_new_all", 2)
# plotting.plot_m31_chi_raw_profile()
# plotting.plot_contour_chi()
# plotting.plot_t23_chi_raw_profile_all_top(savename = "t23_chi_sq_profile_raw_all_top")
# plotting.plot_resolution()

# For debug purposes
top = 2
print("for both")
print("first we get the literature value:")
rate_weight_truth, energy_hist_truth, energy_bins_truth = sensitivity.get_rated_weight_truth(top)
print("chi-sq is", sensitivity.get_chisq(theta23, m31, energy_hist_truth, top))
#sensitivity.get_rated_weight_truth(top = 1)
#sensitivity.get_energy_bins(theta23, m31, top = 1)
print(sensitivity.get_t23_chi_profile(top))

