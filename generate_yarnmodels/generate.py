import time
import multiprocessing as mp
import shutil
import os
import HYLC.PatternSolver as ps
from HYLC.parallel import run_parallel
import HYLC.values
import numpy as np

# some settings
max_steps = 2000  # max solver steps
max_time = 60  # -1 # max solver time spent
ffmt = "%.10g"  # output format


def make_settings(f):
    """initialize solver and material settings by name"""
    s = ps.SimulationSettings()
    s.pypfile = f
    HYLC.values.set_material_settings(s, s.pypfile)
    HYLC.values.set_simulation_settings(s, s.pypfile)

    if("stocktiny" in f or "stocktiniest" in f):
        s.material.yE = 1e6
        s.material.G = 4e5
        s.material.gamma = 1
        s.material.kc = 1e2
        s.material.density = 1.2e3
        s.extpx = 1.6
        s.extpy = 1.8
    return s


def getNverts(f):
    """build a dummy simulation to get number of vertices in the pattern"""
    s = make_settings(f)
    sim = ps.Simulation(s)
    N = len(sim.getQ())
    del sim
    return N


def compute_dQrd(strains, f):
    s = make_settings(f)
    s.strains = strains
    sim = ps.Simulation(s)
    Q0 = sim.getQ()

    time0 = time.time()
    for i in range(max_steps):
        if sim.isFinished():
            break
        if max_time > 0:
            elapsed = time.time() - time0
            if (elapsed > max_time):
                break
        sim.step()

    # Qd = sim.getQdeformed()

    # if all(strains[i] == 0 for i in range(3,6)) and False:
    #     # pure inplane deformation, no need for newton iteration as Jacobian is constant S
    #     S = sim.getSurfaceS()
    #     Qrd = np.copy(Q0)
    #     Qrd[:,3] = Qd[:,3] # copy theta
    #     utilde = sim.getGT()[:,:3] # for R==I, utilde=Rg=g
    #     Qrd[:,:3] += np.einsum('ij,nj->ni',np.linalg.inv(S),utilde) # Xrd = X0 + Sinv utilde
    # else:

    # bending deformation with non-constant Jacobian, and non-trivial inverse mapping deformed->ref
    # do newton iteration on reference space such that it maps to the correct deformed space
    # ie. newton solve on xbar(Xrd) = xdef, for Xrd
    Qrd = sim.computeQrd(eps_r=0.01, max_steps=10)

    del sim
    return Qrd - Q0


# note this needs to be declared outside, for multiprocessing...
def do_thing(args):
    sx, sa, sy, pypfile, i0, i1, i2 = args
    strains = np.zeros(6,)
    strains[0] = sx
    strains[1] = sa
    strains[2] = sy
    # print(strains)
    D = compute_dQrd(strains, pypfile)
    return i0, i1, i2, D


def generate_data(name, pypfile, rgesxsy, rgesa):
    Nv = getNverts(pypfile)

    SX = rgesxsy
    SA = rgesa
    SY = rgesxsy
    Nsims = len(SX)*len(SA)*len(SY)
    # allocate data = [Ntotal x 4], for each vertex and deformation
    YD = np.empty((Nv*Nsims, 4), dtype=np.float32)

    os.makedirs(name+"/sxsasy/", exist_ok=True)
    shutil.copyfile(pypfile, name+"/pyp")

    # write header about deformation
    with open(name+"/sxsasy/axes.txt", 'w') as thefile:
        for ARR in [SX, SA, SY]:
            for i, val in enumerate(ARR):
                thefile.write(ffmt % val)
                if i < len(ARR)-1:
                    thefile.write(" ")
            thefile.write("\n")

    def args():
        for i2 in range(len(SY)):
            for i1 in range(len(SA)):
                for i0 in range(len(SX)):
                    yield SX[i0], SA[i1], SY[i2], pypfile, i0, i1, i2

    def data_ix(vix, isx, isa, isy):
        return isx + (len(SX)) * isa + (len(SX) * len(SA)) * isy + (len(SX) * len(SA) * len(SY)) * vix

    print("Starting.")
    arggen = args()
    p = mp.Pool()
    r = p.imap_unordered(do_thing, arggen)
    count = 0
    for i0, i1, i2, D in r:  # non-parallel processing of unordered incoming results
        print("  %05d/%05d: sx=%.2f sa=%.2f sy=%.2f" %
              (count, Nsims, SX[i0], SA[i1], SY[i2]))
        count += 1
        for vix in range(len(D)):
            YD[data_ix(vix, i0, i1, i2), :] = D[vix, :]

    # serialize folder/sxsasy/data.npy
    np.save(name + "/sxsasy/data.npy", YD, allow_pickle=False)

    print("Done.")


# note this needs to be declared outside, for multiprocessing...
def do_thing_bending(args):
    benddir, val, pypfile, i0 = args
    strains = np.zeros(6,)
    if benddir == "x":
        strains[3] = val
    else:
        strains[5] = val
    D = compute_dQrd(strains, pypfile)
    return i0, D


def generate_bending_data(name, pypfile, IIX, IIY):
    Nv = getNverts(pypfile)

    for benddir, bendarr in [
        ["x", IIX],
        ["y", IIY]
    ]:
        Nsims = len(bendarr)
        # allocate data = [Ntotal x 4], for each vertex and deformation
        YD = np.empty((Nv*Nsims, 4), dtype=np.float32)

        os.makedirs(name+"/bend"+benddir+"/", exist_ok=True)
        # shutil.copyfile(pypfile, name+"/pyp") # assume already existing from sxsasy data

        # write header about deformation
        with open(name+"/bend"+benddir+"/axes.txt", 'w') as thefile:
            for i, val in enumerate(bendarr):
                thefile.write(ffmt % val)
                if i < len(bendarr)-1:
                    thefile.write(" ")
            thefile.write("\n")

        def args():
            for i0 in range(len(bendarr)):
                yield benddir, bendarr[i0], pypfile, i0

        def data_ix(vix, i):
            return i + len(bendarr) * vix

        print("Starting.")
        arggen = args()
        p = mp.Pool()
        r = p.imap_unordered(do_thing_bending, arggen)
        count = 0
        for i0, D in r:  # non-parallel processing of unordered incoming results
            print("  %05d/%05d: II%s=%.2f" %
                  (count, Nsims, benddir, bendarr[i0]))
            count += 1
            for vix in range(len(D)):
                YD[data_ix(vix, i0), :] = D[vix, :]

        # serialize folder/bend/data.npy
        np.save(name+"/bend"+benddir+"/data.npy", YD, allow_pickle=False)

    print("Done.")


def do_thing_bending4D(args):
    sx, sa, sy, val, benddir, pypfile, i0, i1, i2, i3 = args
    strains = np.zeros(6,)
    strains[:3] = sx, sa, sy
    if benddir == "x":
        strains[3] = val
    else:
        strains[5] = val
    D = compute_dQrd(strains, pypfile)
    return i0, i1, i2, i3, D

def generate_bending4D_data(name, pypfile, SX,SA,SY,IIX,IIY):
    Nv = getNverts(pypfile)

    for benddir, bendarr in [
        ["x", IIX],
        ["y", IIY]
    ]:
        Nsims = len(bendarr)*len(SX)*len(SA)*len(SY)
        # allocate data = [Ntotal x 4], for each vertex and deformation
        YD = np.empty((Nv*Nsims, 4), dtype=np.float32)

        os.makedirs(name+"/bend4D"+benddir+"/", exist_ok=True)
        shutil.copyfile(pypfile, name+"/pyp") # assume already existing from sxsasy data

        # write header about deformation
        with open(name+"/bend4D"+benddir+"/axes.txt", 'w') as thefile:
            for ARR in [SX, SA, SY, bendarr]:
                for i, val in enumerate(ARR):
                    thefile.write(ffmt % val)
                    if i < len(ARR)-1:
                        thefile.write(" ")
                thefile.write("\n")

        def args():
            for i3 in range(len(bendarr)):
                for i2 in range(len(SY)):
                    for i1 in range(len(SA)):
                        for i0 in range(len(SX)):
                            yield SX[i0], SA[i1], SY[i2], bendarr[i3], benddir, pypfile, i0, i1, i2, i3

        def data_ix(vix, isx, isa, isy, ibend):
            return isx + (len(SX)) * isa + (len(SX) * len(SA)) * isy + (len(SX) * len(SA) * len(SY)) * ibend + (len(SX) * len(SA) * len(SY) * len(bendarr)) * vix

        print("Starting.")
        arggen = args()
        p = mp.Pool()
        r = p.imap_unordered(do_thing_bending4D, arggen)
        count = 0
        for i0, i1, i2, i3, D in r:  # non-parallel processing of unordered incoming results
            print("  %05d/%05d: sx=%.2f sa=%.2f sy=%.2f II%s=%.2f" %
                (count, Nsims, SX[i0], SA[i1], SY[i2], benddir, bendarr[i3]))
            count += 1
            for vix in range(len(D)):
                YD[data_ix(vix, i0, i1, i2, i3), :] = D[vix, :]

        # serialize folder/bend/data.npy
        np.save(name+"/bend4D"+benddir+"/data.npy", YD, allow_pickle=False)

    print("Done.")
