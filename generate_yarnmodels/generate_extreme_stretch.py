import generate
import numpy as np
import time


SX = np.concatenate([np.linspace(-0.2,0.0,15 - 15*4//6), np.linspace(0.0,3.0,15*4//6 + 1)[1:]])
SA = np.array([0.0,0.01]) # dummy
SY = np.array([0.0,0.01]) # dummy

timings = []
timingsfile = 'timings_%s.txt' % time.strftime("%Y%m%d-%H%M%S")

for name, pypfile in [
    # ['stocktiniest', "pyp/stocktiniest.pyp"],
    ['stock', "pyp/hylc2020/stockinette.pyp"],
    ['rib', "pyp/hylc2020/cartridge_belt_rib.pyp"],
    # ['honey', "pyp/hylc2020/slip_stitch_honeycomb.pyp"],
    ['basket', "pyp/hylc2020/basket.pyp"],
    # ['satin', "pyp/hylc2020/satin.pyp"],
]:
  tstname = "stretchx"
  print("Pattern:",name)
  print("Samples: %d x %d x %d (Total: %d)" % (len(SX),len(SA),len(SY), len(SX)*len(SA)*len(SY)))
  t0 = time.perf_counter()
  generate.generate_data("model_%s_%s" % (name,tstname), pypfile, SX,
                        SA, SY, stretch_test=True)
  dt = time.perf_counter() - t0
  timings.append(["model_%s_%s" % (name,tstname), dt])
  with open(timingsfile, 'a') as f1:
    f1.write("%s %s\n" % ("model_%s_%s" % (name,tstname), dt))
