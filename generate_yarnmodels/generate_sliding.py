import generate
import numpy as np
import time


rge5xy = [-0.2,0.0,0.33,0.67,1.0]
rge5a = [-0.5,-0.25,0.0,0.25,0.5]
rge9xy = np.concatenate([np.linspace(-0.2,0.0,9 - 9*4//6), np.linspace(0.0,1.0,9*4//6 + 1)[1:]])
rge9a = np.linspace(-0.7,0.7,9)
rge15xy = np.concatenate([np.linspace(-0.2,0.0,15 - 15*4//6), np.linspace(0.0,1.0,15*4//6 + 1)[1:]])
rge15a = np.linspace(-0.7,0.7,15)
rge31xy = np.concatenate([np.linspace(-0.2,0.0,31 - 31*4//6), np.linspace(0.0,1.0,31*4//6 + 1)[1:]])
rge31a = np.linspace(-0.7,0.7,31)


for name, pypfile in [
    # ['stocktiniest', "pyp/stocktiniest.pyp"],
    ['stock', "pyp/hylc2020/stockinette.pyp"],
    # ['rib', "pyp/hylc2020/cartridge_belt_rib.pyp"],
    # ['honey', "pyp/hylc2020/slip_stitch_honeycomb.pyp"],
    # ['basket', "pyp/hylc2020/basket.pyp"],
    # ['satin', "pyp/hylc2020/satin.pyp"],
]:
  for tstname,rgesxsy,rgesa in [
    ['15slide',rge15xy,rge15a],
  ]:
    print("Pattern:",name)
    print("Samples: %d x %d x %d (Total: %d)" % (len(rgesxsy),len(rgesa),len(rgesxsy), len(rgesxsy)*len(rgesxsy)*len(rgesa)))
    generate.generate_data("model_%s_%s" % (name,tstname), pypfile, rgesxsy,
                          rgesa, rgesxsy, disable_slide_constr=True)