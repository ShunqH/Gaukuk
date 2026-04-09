import numpy as np

class ReadGaukuk:
    def __init__(self, filename, isCons = True, dtype=np.float64):
        self.filename = filename
        self.dtype = dtype

        with open(filename, "rb") as f:
            RealType = np.fromfile(f, dtype=np.int32, count=1)[0]
            if RealType == 8:
                dtype = np.float64
            elif RealType == 4:
                dtype = np.float32

            # frame information
            self.t = np.fromfile(f, dtype=dtype, count=1)[0]

            # header
            self.nvar = np.fromfile(f, dtype=np.int32, count=1)[0]
            self.nx = np.fromfile(f, dtype=np.int32, count=1)[0]
            self.ny = np.fromfile(f, dtype=np.int32, count=1)[0]
            self.nz = np.fromfile(f, dtype=np.int32, count=1)[0]
            self.nGhost = np.fromfile(f, dtype=np.int32, count=1)[0]
            self.lenx = self.nx + 2 * self.nGhost 
            self.leny = self.ny + 2 * self.nGhost
            self.lenz = self.nz + 2 * self.nGhost
            self.lenArr = self.lenz * self.leny * self.lenx 

            # mesh
            self.xc = np.fromfile(f, dtype=dtype, count=self.nx)
            self.yc = np.fromfile(f, dtype=dtype, count=self.ny)
            self.zc = np.fromfile(f, dtype=dtype, count=self.nz)

            # data
            self.dataLen = np.fromfile(f, dtype=np.int32, count=1)[0] 
            if (self.dataLen != self.nvar * self.lenArr ):
                print("data size error, reshape might fail")
            self.data = np.fromfile(f, dtype=dtype, count=self.dataLen)
            self.data = self.data.reshape(
                    (self.nvar, self.lenz, self.leny, self.lenx)
                )