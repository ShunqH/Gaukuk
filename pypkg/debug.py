import numpy as np 

def ReadTArray(filename):
    with open(filename, "rb") as f:
        type= np.fromfile(f, dtype=np.int32, count=1)[0]
        size = np.fromfile(f, dtype=np.int32, count=1)[0] 
        dtype=np.float64
        if type== 8:
            dtype=np.float64
        elif type== 4:
            dtype=np.float32
        data = np.fromfile(f, dtype=dtype, count=size)
    return data

nVar = 5 
nx = 128
ny = 16
nz = 16 
nGhost = 1
lenx = nx + 2 * nGhost + 1
leny = ny + 2 * nGhost
lenz = nz + 2 * nGhost
lenArr = lenx * leny * lenz

path = "../bin/"
frameID = 3
filename = path + "flx1_"+str(frameID).zfill(5)
test = ReadTArray(filename)
test = test.reshape((nVar, lenz, leny, lenx))
print(test[0,9,5,:])

