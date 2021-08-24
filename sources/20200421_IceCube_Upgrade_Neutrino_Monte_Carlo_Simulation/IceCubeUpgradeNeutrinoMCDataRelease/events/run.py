import plotting
import sensitivity
from params import *

# set the topology. This will be integrated into shell interface later
top = 0

# set the truth energy bins
rate_weight_truth, energy_hist_truth, energy_bins_truth = sensitivity.get_rated_weight_truth(top)

# plotting.distribution_plots()
#print("this is output file for both plot")
#plotting.plot_t23_chi_raw_profile("t23_raw_profile_new_both", 2)
# plotting.plot_t23_chi_raw_profile("t23_raw_profile_new_track", 1)
# plotting.plot_t23_chi_raw_profile("t23_raw_profile_new_all", 2)
# plotting.plot_m31_chi_raw_profile()
# plotting.plot_contour_chi()
# plotting.plot_t23_chi_raw_profile_all_top(savename = "t23_chi_sq_profile_raw_all_top")
# plotting.plot_resolution()\





# plotting.plot_t23_chi_raw_profile(energy_hist_truth, savename = "t23_chi_sq_profile_raw_new_new", top = 0)
sensitivity.stupidity(energy_hist_truth, top)

# For debug purposes
#top = 0
#print("for cascade")
#print("first we get the literature value:")
# print("first define the truth values, should get into get_rated_weight_truth only")
#rate_weight_truth, energy_hist_truth, energy_bins_truth = sensitivity.get_rated_weight_truth(top)
#sensitivity.get_rated_weight_truth(top = 1)
#sensitivity.get_energy_bins(theta23, m31, top = 1)
print("now running same thing with get chi profile")
print(sensitivity.get_t23_chi_profile(energy_hist_truth, top = top))

# print("now let us go into PROFILE function")
#sensitivity.get_t23_chi_profile(energy_hist_truth, m31 = m31, top = 0)

# print("see if this is the same as just calling the get_chisq function")
# sensitivity.get_chisq(t23l[0], m31, energy_hist_truth, top)
# 
# print("now with the old script")
# sensitivity.get_t23_chi_profile_old(m31 = m31, top = 0)
