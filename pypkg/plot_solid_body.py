import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from read_gaukuk import ReadGaukuk
import os 

DEN = 0 
MTX = 1 
MTY = 2 
MTZ = 3 
ENG = 4 
gamma = 1.4 
gm1Rec = 1.0 / (gamma - 1)
totalFrames = 51

tlist = np.zeros(totalFrames) 
englist = np.zeros(totalFrames)
sPath = "/Users/shunq/Documents/CPP/cal_gaukuk/solid_body/3Dm=4/"
tag = "solid"

for frameID in range(31, totalFrames):
    filename = "../bin/cons_" + str(frameID).zfill(5)
    savename = sPath + tag + "_" + str(frameID).zfill(5) + ".png"
    frame = ReadGaukuk(filename)

    k = int((frame.lenz)/2)
    kl = frame.nGhost
    kr = frame.nGhost + frame.nz 
    jl = frame.nGhost
    jr = frame.nGhost + frame.ny 
    il = frame.nGhost
    ir = frame.nGhost + frame.nx
    rho = frame.data[DEN,k,jl:jr,il:ir]
    mtx = frame.data[MTX,k,jl:jr,il:ir]
    mty = frame.data[MTY,k,jl:jr,il:ir]
    mtz = frame.data[MTZ,k,jl:jr,il:ir]
    eng = frame.data[ENG,k,jl:jr,il:ir]

    x = frame.xc
    y = frame.yc
    rhoMin = 1e-2
    rhoMax = 5

    # print(il,ir,jl,jr)
    X, Y = np.meshgrid(x, y)

    # cmap='viridis'
    # cmap='magma'
    cmap='inferno'
    # cmap='cividis'
    # cmap='cubehelix'

    plt.figure(figsize=(10, 8))
    fig, ax = plt.subplots(1, 1)
    plt.pcolormesh(X, Y, rho, cmap=cmap, shading='gouraud',
                    # vmin=rhoMin, vmax=rhoMax, 
                    norm=LogNorm(vmin=rhoMin, vmax=rhoMax,), 
                    )
    circle = plt.Circle((-3, 0), 1, color='black', fill=True)
    ax.add_patch(circle)
    plt.colorbar(label='Density')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title("t = " + f"{(frame.t):.2f}", fontsize = 18)
    plt.axis('equal')   
    plt.savefig(savename, bbox_inches='tight', dpi=300)
    # plt.show()
    plt.close()
    print("frame = ", frameID, end="\r")
    del frame


os.system("ffmpeg -y -framerate 10 -i " + sPath + tag + "_%05d.png -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -c:v libx264 -pix_fmt yuv420p -crf 23 " + sPath + tag + ".mp4")

