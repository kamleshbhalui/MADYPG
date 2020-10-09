
import sys
import os
import subprocess

assert(len(sys.argv) == 4)

executable = sys.argv[1]
parentfolder = sys.argv[2]
outfolder = sys.argv[3]

for subfolder in os.listdir(parentfolder):
  if not os.path.isdir(os.path.join(parentfolder,subfolder)):
    continue

  infolder = os.path.join(parentfolder,subfolder)
  outfile = os.path.join(outfolder,subfolder + ".bin")
  os.makedirs(outfolder,exist_ok=True)
  print(subfolder, "->", outfile)
  print([executable, infolder, outfile])
  subprocess.run([executable, infolder, outfile])

# TODO ALTERNATIVELY: RECURSIVE FOLDER ITERATION, AND IF A FOLDER CONTAINS ANY .obj THEN DONT GO FURTHER BUT CONVERT