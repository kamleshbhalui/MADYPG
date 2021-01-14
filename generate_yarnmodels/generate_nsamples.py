import generate
import numpy as np

# maxish = 1.0
# minish = -0.25
# rge = (np.linspace(0, 1, 24)**2 * (maxish - minish) + minish)
# rge = rge - rge[np.argmin(np.abs(rge))]
# rgexy = np.round(rge, 3)
# rgea = np.round(np.linspace(-0.7, 0.7, 23),  3)
# print(rgexy)
# print(rgea)
# print(len(rgexy)*len(rgexy)*len(rgea))


rge5xy = [-0.2,0.0,0.33,0.67,1.0]
rge5a = [-0.5,-0.25,0.0,0.25,0.5]
rge9xy = np.concatenate([np.linspace(-0.2,0.0,9 - 9*4//6), np.linspace(0.0,1.0,9*4//6 + 1)[1:]])
rge9a = np.linspace(-0.7,0.7,9)
rge15xy = np.concatenate([np.linspace(-0.2,0.0,15 - 15*4//6), np.linspace(0.0,1.0,15*4//6 + 1)[1:]])
rge15a = np.linspace(-0.7,0.7,15)
rge31xy = np.concatenate([np.linspace(-0.2,0.0,31 - 31*4//6), np.linspace(0.0,1.0,31*4//6 + 1)[1:]])
rge31a = np.linspace(-0.7,0.7,31)

# print(rge9xy)
# print(rge9a)
# print(rge15xy)
# print(rge15a)
# print(rge31xy)
# print(rge31a)

for name, pypfile in [
    # ['stocktiniest', "pyp/stocktiniest.pyp"],
    ['stock', "pyp/hylc2020/stockinette.pyp"],
    ['rib', "pyp/hylc2020/cartridge_belt_rib.pyp"],
    # ['honey', "pyp/hylc2020/slip_stitch_honeycomb.pyp"],
    # ['basket', "pyp/hylc2020/basket.pyp"],
    # ['satin', "pyp/hylc2020/satin.pyp"],
]:
  for tstname,rgesxsy,rgesa in [
    ['5',rge5xy,rge5a],
    ['9',rge9xy,rge9a],
    ['15',rge15xy,rge15a],
    ['31',rge31xy,rge31a],
  ]:
      print("Test %s will have %d total samples"%(tstname,len(rgesxsy)*len(rgesxsy)*len(rgesa)))
      generate.generate_data("model_%s_%s" % (name,tstname), pypfile, rgesxsy,
                            rgesa, rgesxsy)
