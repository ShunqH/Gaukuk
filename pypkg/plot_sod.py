import numpy as np 
import matplotlib.pyplot as plt
from read_gaukuk import ReadGaukuk

DEN = 0 
MTX = 1 
MTY = 2 
MTZ = 3 
ENG = 4 
gamma = 1.4 
gm1Rec = 1.0 / (gamma - 1)

frameID = 0
filename = "../bin/prim_" + str(frameID).zfill(5)
frame = ReadGaukuk(filename)
dx = frame.xc[1] - frame.xc[0]
xlist = np.linspace(frame.xc[0]-dx, frame.xc[-1]+dx, frame.lenx)

k = int(frame.nz/2)
j = int(frame.ny/2)
# k = 3
# j = 3

fig = plt.figure(figsize=(12,8))
plt.subplots_adjust(hspace=0.02, wspace=0.2)
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)

for frameID in range(0,21,5):
    filename = "../bin/prim_" + str(frameID).zfill(5)
    frame = ReadGaukuk(filename)
    rho = frame.data[DEN,k,j,:]
    vx = frame.data[MTX,k,j,:]
    vy = frame.data[MTY,k,j,:]
    vz = frame.data[MTZ,k,j,:]
    pre = frame.data[ENG,k,j,:] 
    eng = pre * gm1Rec + 0.5 * rho * (vx*vx + vy*vy + vz*vz)
    eng = pre * gm1Rec / rho 
    ax1.plot(xlist, rho, label = f"t = {frame.t:.2f}")
    ax2.plot(xlist, pre)
    ax3.plot(xlist, vx)
    ax4.plot(xlist, eng)


ax1.set_xlim(-0.5, 0.5)
ax1.set_ylim(0, 1.1)
ax1.grid()
ax1.legend(fontsize = 15)
ax1.set_xlabel(r"$x$", fontsize = 15)
ax1.set_ylabel(r"$\rho$", fontsize = 15)

ax2.set_xlim(-0.5, 0.5)
ax2.set_ylim(0, 1.1)
ax2.grid()
ax2.set_xlabel(r"$x$", fontsize = 15)
ax2.set_ylabel(r"$P$", fontsize = 15)

ax3.set_xlim(-0.5, 0.5)
ax3.set_ylim(0, 1.1)
ax3.grid()
ax3.set_xlabel(r"$x$", fontsize = 15)
ax3.set_ylabel(r"$v_x$", fontsize = 15)

ax4.set_xlim(-0.5, 0.5)
ax4.set_ylim(1.5, 2.98)
ax4.grid()
ax4.set_xlabel(r"$x$", fontsize = 15)
ax4.set_ylabel("Specific Internal Energy", fontsize = 15)

plt.savefig("sod.png", bbox_inches='tight', dpi=300)
# plt.show()
plt.close()