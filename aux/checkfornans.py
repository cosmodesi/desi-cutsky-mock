import numpy as np
import fitsio
import glob
from combiconvert_rdz2xyz import LC
import argparse



parser = argparse.ArgumentParser()
parser.add_argument("--ngc_sgc", type=str, help="NGC or SGC preferred rotation")
args = parser.parse_args()

files = glob.glob("//global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CubicBox/ELG/z1.100/EZmock_B2000G512Z1.1N24000470_b0.345d1.45r40c0.05_seed56/*gz")
print(len(files))

# list_ = ["x", "y", "z", "vx", "vy", "vz"]
list_ = ["vz"]

data = {"x": np.empty(0),"y": np.empty(0),"z": np.empty(0)}

for i, file_ in enumerate(files):
    fits = fitsio.FITS(file_, "r")

    for col in list_:
        index_ = np.isnan(fits[1][col][:])
        idx = np.arange(len(fits[1][col][:]))

        names, counts = np.unique(index_, return_counts=True)
    
        if True in names:
            print(i, "/", len(files), ":", col, names, counts)
            # print(fits[1]["x"][idx[index_]], fits[1]["y"][idx[index_]], fits[1]["z"][idx[index_]], fits[1]["vx"][idx[index_]], fits[1]["vy"][idx[index_]], fits[1]["vz"][idx[index_]])
            data["x"] = np.append(data["x"], fits[1]["x"][idx[index_]][0])
            data["y"] = np.append(data["y"], fits[1]["y"][idx[index_]][0])
            data["z"] = np.append(data["z"], fits[1]["z"][idx[index_]][0])
    # print()
    fits.close()

print(data)


convert = LC("./EZmock/config/config_EZmock_ELG.ini", args)

shellnums = convert.compute_shellnums()

for shellnum in shellnums:
    convert.convert_xyz2rdz(data, shellnum)

print("NGC")

fits = fitsio.FITS("/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky/ELG/z1.100/cutsky_ELG_z1.100_EZmock_B2000G512Z1.1N24000470_b0.345d1.45r40c0.05_seed56.fits", "r")

index_ = np.isnan(fits[1]["Z"][:])
idx = np.arange(len(fits[1]["Z"][:]))


print(np.unique(index_, return_counts=True))
print(fits[1]["RA"][idx[index_]], fits[1]["DEC"][idx[index_]], fits[1]["Z"][idx[index_]])
fits.close()
