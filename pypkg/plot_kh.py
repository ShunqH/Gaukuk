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

for frameID in range(21):
    filename = "../bin/prim_" + str(frameID).zfill(5)
    savename = "kh_" + str(frameID).zfill(5) + ".png"
    frame = ReadGaukuk(filename)

    k = int(frame.nz/2)
    jl = frame.nGhost
    jr = frame.nGhost + frame.ny 
    il = frame.nGhost
    ir = frame.nGhost + frame.nx
    rho = frame.data[DEN,k,jl:jr,il:ir]
    vx = frame.data[MTX,k,jl:jr,il:ir]
    vy = frame.data[MTY,k,jl:jr,il:ir]
    vz = frame.data[MTZ,k,jl:jr,il:ir]
    pre = frame.data[ENG,k,jl:jr,il:ir]

    x = frame.xc
    y = frame.yc

    # print(il,ir,jl,jr)
    X, Y = np.meshgrid(x, y)

    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X, Y, rho, cmap='viridis', shading='auto')
    plt.colorbar(label='Density')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Density')
    plt.axis('equal')   
    plt.savefig(savename, bbox_inches='tight', dpi=300)
    # plt.show()
    plt.close()

    del frame