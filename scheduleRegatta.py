import os

nCoresList = [1, 4, 16]
sizesList = [16, 32, 64, 128, 256, 512, 1024, 2048]


homepath = "/home/edu-cmc-skpod17-325-05/"

for nCores in nCoresList:
    for sz in sizesList:
        for v in range(3):
            fName = "threads" + str(nCores) + "size" + str(sz) + "v" + str(v) + ".out"
            os.system("ompsubmit -n " + "16" + " -w 15:00" + " --stdout " + fName + " prog -- " + str(sz) + " " + str(nCores))
