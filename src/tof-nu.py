#Time of flight from periapsis give nu

import tof_tool as toolbox
import numpy as np
import math as m
import matplotlib.pyplot as plt

#Problem Inputs
r0 = np.array([1,0,0])
v0 = np.array([0,1.1,0])
nu = 6*m.pi
mu = 1

#Calculate orbital parameters
a,p,e,i,ascend,peri,SME,conic = toolbox.rv_orbparams(r0[0], r0[1], r0[2], v0[0], v0[1], v0[2], 1)

full_orb = (nu - (nu % (2*m.pi))) / (2*m.pi)
print("full orb ", full_orb)


nu_set = np.linspace(0, 4*m.pi, 50)
rem_nu_set = []
Eset = []
tof_set = []
for nu in nu_set:
    rem_nu_set.append(toolbox.rev_and_rem(nu)[1])
    Ecurr = toolbox.tr_to_ecc_anomaly(e,nu)
    Eset.append(Ecurr)
    tof_set.append(toolbox.tof(a,mu, nu, Ecurr,e))



# fig, ax = plt.subplots()
# ax.plot(nu_set, Eset)
# ax.plot(nu_set, tof_set)
# ax.set(xlabel="True Anomaly", ylabel="Eccentric Anomaly/ToF", title="Nu vs. E")
# ax.grid()
# plt.show()





fig, ax = plt.subplots()
ax.plot(tof_set, rem_nu_set)
ax.set(xlabel='time (DU)', ylabel='True Anomaly (rad)', title='True Anomaly vs. Time')
ax.grid()
plt.show()
