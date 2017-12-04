import os

nCoresList = [1, 4]
sizesList = [32, 64, 128, 256, 512, 1024]

for nCores in nCoresList:
    for sz in sizesList:
        for v in range(3):
            fName = "threads" + str(nCores) + "size" + str(sz) + "v" + str(v) + ".out"
            os.system("mpisubmit.bg -n 1" + " -stdout " + fName +
                      " prog -- " + str(sz) + " " + str(nCores))
