"""
Created on Fri Dec 23 13:46:31 2022

@author: ronan
"""

import json
import numpy as np
import matplotlib.pyplot as plt
file = open("sb_radar.json")
vel = []
pos = []
data = json.load(file)
for i in data["data"]:
    if i[4] == "Hz":
        vel.append(300*float(i[3])/float(i[5]))
    else :
        pos.append(float(i[3])*300/2)

file.close()
plt.hist(vel,bins=40)
plt.show()
plt.hist(pos,bins=40)
plt.show()

velplt = []
posplt = []

for i in vel:
    if i>0.5:
        print(i)
    else:
        velplt.append(i)
for i in pos:
    if i>3000:
        print(i)
    else:
        posplt.append(i)
        
plt.hist(velplt,bins=60)
plt.xlim(0,0.5)
plt.xlabel("velocity (m/s)")
plt.ylabel("number of occurrences")
plt.savefig("velocity_uncertainties.pdf",dpi=600)
plt.show()
print("mean velocity uncertainty",np.mean(velplt))
plt.hist(posplt,bins=60)
plt.xlim(0,3000)
plt.xlabel("distance (m)")
plt.ylabel("number of occurrences")
plt.savefig("position_uncertainties.pdf",dpi=600)
print("mean position uncertainty",np.mean(posplt))
plt.show()