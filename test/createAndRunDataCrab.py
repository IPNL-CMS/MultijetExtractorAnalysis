#! /usr/bin/env python

import os, datetime, pwd

from optparse import OptionParser
parser = OptionParser()
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
(options, args) = parser.parse_args()

datasets = [
["/Jet/apequegn-Jet_Run2012A-22Jan2013_04Dec13-v1-2b71cc75519c435218af9c37f5843477/USER", "Jet_Run2012A-22Jan2013", "FT53_V21A_AN6"],
["/JetHT/apequegn-JetHT_Run2012B-22Jan2013_04Dec13-v1-2b71cc75519c435218af9c37f5843477/USER", "JetHT_Run2012B-22Jan2013", "FT53_V21A_AN6"],
["/JetHT/apequegn-JetHT_Run2012C-22Jan2013_04Dec13-v1-2b71cc75519c435218af9c37f5843477/USER", "JetHT_Run2012C-22Jan2013", "FT53_V21A_AN6"],
["/JetHT/apequegn-JetHT_Run2012D-22Jan2013_04Dec13-v1-2b71cc75519c435218af9c37f5843477/USER", "JetHT_Run2012D-22Jan2013", "FT53_V21A_AN6"]
    ]

# Get email address
email = "%s@ipnl.in2p3.fr" % (pwd.getpwuid(os.getuid()).pw_name)

d = datetime.datetime.now().strftime("%d%b")

version = 1

print("Creating configs for crab. Today is %s, you are %s and it's version %d" % (d, email, version))
print("")

for dataset in datasets:
  dataset_path = dataset[0]
  dataset_name = dataset[1]  + "_woPU_pt25_eta50_PUtightOnly"
  dataset_globaltag = dataset[2]

  ui_working_dir = ("crab_data_%s_%s") % (dataset_name, d)
  output_file = "crab_data_%s_%s.cfg" % (dataset_name, d)
  output_dir = ("Extracted_step2/data/%s/%s" % (d, dataset_name))

  python_config = "Extractor_MULTIJET_data.py";
#  if "Electron" in dataset_path:
#    python_config = "Extractor_MTT_semie.py"
#  else:
#    python_config = "Extractor_MTT_semimu.py"

  print("Creating config file for '%s'" % (dataset_path))
  print("\tName: %s" % dataset_name)
  print("\tGlobal tag: %s" % dataset_globaltag)
  print("\tOutput directory: %s" % output_dir)
  print("")
    

  os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@email@#%s#g\" -e \"s#@globaltag@#%s#g\" -e \"s#@pset@#%s#g\" -e \"s#@outputdir@#%s#\" crab_data.cfg.template.ipnl > %s" % (dataset_path, ui_working_dir, email, dataset_globaltag, python_config, output_dir, output_file))

  if options.run:
    cmd = "crab -create -submit -cfg %s" % (output_file)
    os.system(cmd);
