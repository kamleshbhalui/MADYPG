import generate
import numpy as np

nsamples = 31
for name, pypfile, bendx, bendy in [
    ['stock', "pyp/hylc2020/stockinette.pyp", [-200, 200], [-230, 230]],
    ['rib', "pyp/hylc2020/cartridge_belt_rib.pyp", [-100, 100], [-80,  80]],
    # ['honey', "pyp/hylc2020/slip_stitch_honeycomb.pyp", [-120, 120], [-120, 120]],
    # ['basket', "pyp/hylc2020/basket.pyp", [-170, 170], [-170, 170]],
    # ['satin', "pyp/hylc2020/satin.pyp", [-170, 170], [-170, 170]],
]:
    IIX = np.linspace(bendx[0],bendx[1], nsamples)
    IIX -= IIX[np.abs(IIX).argmin()]
    IIY = np.linspace(bendy[0],bendy[1], nsamples)
    IIY -= IIY[np.abs(IIY).argmin()]
    print("Generating data for", name, "with:\n",IIX, IIY)
    generate.generate_bending_data("model_%s_bend" % (name), pypfile, IIX, IIY)
