"""
   GBM_exp.py
   Author: Ginny Cunningham
   Date: December 11, 2017
   For a given magnitude and time of a GRB, calculate the expected magnitude at a later time assuming a power law decay. 
   Usage:  python GBM_exp.py [Initial_Magnitude] [Age of Burst] 
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
print " "

#####Inputs
inputs = sys.argv
m0 = float(inputs[1]) #Initial magnitude (read from email alert)
t0 = float(inputs[2]) #Time since burst (read from email alert) [s]
#delta_t = float(inputs[3]) #Expected delta_t from t0 to time of observation with SEDm [s]
delta_t = 600. #Use a placeholder value of 10 minutes for now [s]
t_obs = t0+delta_t #Expected time of observation [s]
#print m0,t0,t_obs

#####Power law index
gamma = 1.

#####Calculate expected magnitude
tend = 10*t_obs #End of observation run [s]
t = np.linspace(t0, tend, 1000) #Range of times to observe over [s] 
m = 2.5*gamma*np.log10(t/t0)+m0 #final magnitude
m_exp = 2.5*gamma*np.log10((t_obs)/t0)+m0 #Expected magnitude at t_obs=t0+delta_t
print "Expected Magnitude %s s after initial UVOT observations: %.2f" %(delta_t, m_exp)

#####Examples of exposure times for various magnitudes
if 10 < m_exp <= 11:
    exptime = 120
elif 11 < m_exp <= 12:
    exptime = 240
elif 12 < m_exp <= 13:
    exptime = 360
elif 13 < m_exp <= 14:
    exptime = 500
else:
    exptime = -1
    print "Exposure Time not Within Expected Range."
print "Recommended Exposure time: %s seconds" %exptime

#####Plotting
plt.semilogx(t,m)
plt.gca().invert_yaxis()
p1 = plt.scatter(t0, m0, c='g', label="Initial Magnitude")
p2 = plt.scatter(t0+delta_t, m_exp, c='r', label="Expected Magnitude at Observation")
plt.legend()
plt.xlabel("Time since Trigger [s]")
plt.ylabel("Magnitude")
plt.show()

print " "


