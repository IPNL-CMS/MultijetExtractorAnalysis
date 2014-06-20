#! /usr/bin/env python
# Launch crab for every datasets in mc_signal_datasets.list

import os, shutil, subprocess
from optparse import OptionParser

isCastor = os.system("uname -n | grep cern &> /dev/null") == 0

parser = OptionParser()
parser.add_option("-p", "--path", dest="path", type="string", help="where to store crab folders")

(options, args) = parser.parse_args()

if options.path is None or not os.path.isdir(options.path):
  parser.error("you must specify a valid path")

crabFolders = [name for name in os.listdir(options.path) if os.path.isdir(os.path.join(options.path, name)) and name.startswith("crab_MC_")]

for crabFolder in crabFolders:
  dataset = crabFolder.rstrip("/").replace("crab_MC_", "")
  # Remove date
  dataset_name = dataset
  print("Processing %s" % dataset)
  outputName = "MULTIJET_%s.list" % (dataset)
  fullPath = "%s/%s" % (options.path, crabFolder)
  if os.path.exists(outputName):
    print("'%s' already exists. Skipping." % outputName)
    continue

  p = subprocess.Popen(["crabOutputList.py", "--check", "--xrootd", fullPath], stdout=subprocess.PIPE)
  dpmFiles = [line.replace("\n", "") for line in p.stdout.readlines()]
  p.wait()

  if p.returncode != 0:
    print("Error: can't merge for %s because crabOutputList was not successfull" % dataset)
    continue

  file = open(outputName, "w")
  for f in dpmFiles:
    file.write(f + "\n")
  file.close()

    #singleLineFiles = ""
    #for f in dpmFiles:
    #  singleLineFiles = "%s%s " % (singleLineFiles, f)

    #os.system("hadd %s %s" % (outputName, singleLineFiles))
    #print("hadd %s %s" % (outputName, singleLineFiles))
