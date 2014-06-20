#! /usr/bin/env python

import os, copy, datetime, pwd, re

from optparse import OptionParser
parser = OptionParser()
parser.add_option("", "--create-cfg", action="store_true", dest="create_cfg", default=False, help="create config files for crab")
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
parser.add_option("", "--status", action="store_true", dest="status", default=False, help="run crab -status")
parser.add_option("", "--get", action="store_true", dest="get", default=False, help="run crab -get")
parser.add_option("", "--resubmit", action="store_true", dest="resubmit", default=False, help="run crab -resubmit bad")
parser.add_option("", "--submit", action="store_true", dest="submit", default=False, help="run crab -submit all")
parser.add_option("", "--kill", action="store_true", dest="kill", default=False, help="run crab -kill all")
(options, args) = parser.parse_args()

if options.run:
  options.create_cfg = True

datasets = [
["/Jet/apequegn-Jet_Run2012A-22Jan2013_09Apr14-v1-829b2d54640b0ff2d246edb621ac7ffd/USER" , "Jet_Run2012A-22Jan2013", "FT53_V21A_AN6"],
["/JetHT/apequegn-JetHT_Run2012B-22Jan2013_09Apr14-v1-9f561dc7a7b870f3c507ce92c9c7f820/USER", "JetHT_Run2012B-22Jan2013", "FT53_V21A_AN6"],
["/JetHT/apequegn-JetHT_Run2012C-22Jan2013_09Apr14-v1-9f561dc7a7b870f3c507ce92c9c7f820/USER", "JetHT_Run2012C-22Jan2013", "FT53_V21A_AN6"],
["/JetHT/apequegn-JetHT_Run2012D-22Jan2013_09Apr14-v1-9f561dc7a7b870f3c507ce92c9c7f820/USER", "JetHT_Run2012D-22Jan2013", "FT53_V21A_AN6"],
["/JetMon/apequegn-JetMon_Run2012B-22Jan2013_09Apr14-v1-3ba8c59f861a2218f0d0749eda03cb49/USER", "JetMon_Run2012B-22Jan2013", "FT53_V21A_AN6"],
["/JetMon/apequegn-JetMon_Run2012C-22Jan2013_09Apr14-v1-3ba8c59f861a2218f0d0749eda03cb49/USER", "JetMon_Run2012C-22Jan2013", "FT53_V21A_AN6"],
["/JetMon/apequegn-JetMon_Run2012D-22Jan2013_09Apr14-v1-3ba8c59f861a2218f0d0749eda03cb49/USER", "JetMon_Run2012D-22Jan2013", "FT53_V21A_AN6"]
    ]

# Get email address
email = "%s@ipnl.in2p3.fr" % (pwd.getpwuid(os.getuid()).pw_name)

d = datetime.datetime.now().strftime("%d%b")

version = 1

print("Creating configs for crab. Today is %s, you are %s and it's version %d" % (d, email, version))
print("")

for dataset in datasets:
  dataset_path = dataset[0]
  dataset_name = dataset[1]  + "_woPU_pt30_eta50_puJetIdT"
  dataset_globaltag = dataset[2]

  ui_working_dir = ("crab_data_%s_%s") % (dataset_name, d)
  #ui_working_dir = ("crab_data_%s_16Jun") % (dataset_name)
  output_file = "crab_data_%s_%s.cfg" % (dataset_name, d)
  output_dir = ("Extracted_step2/data/%s/%s" % (d, dataset_name))
  #output_dir = ("Extracted_step2/data/16Jun/%s" % (dataset_name))

  python_config = "Extractor_MULTIJET_data.py";
#  if "Electron" in dataset_path:
#    python_config = "Extractor_MTT_semie.py"
#  else:
#    python_config = "Extractor_MTT_semimu.py"

  if options.create_cfg:
    print("Creating config file for '%s'" % (dataset_path))
    print("\tName: %s" % dataset_name)
    print("\tGlobal tag: %s" % dataset_globaltag)
    print("\tOutput directory: %s" % output_dir)
    print("")
    
    os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@email@#%s#g\" -e \"s#@globaltag@#%s#g\" -e \"s#@pset@#%s#g\" -e \"s#@outputdir@#%s#\" crab_data.cfg.template.ipnl > %s" % (dataset_path, ui_working_dir, email, dataset_globaltag, python_config, output_dir, output_file))

  if options.run:
    cmd = "crab -create -submit -cfg %s" % (output_file)
    os.system(cmd)

  if options.status:
    cmd = "crab -status -c %s" % (ui_working_dir)
    os.system(cmd)

  if options.get:
    cmd = "crab -get -c %s" % (ui_working_dir)
    os.system(cmd)

  if options.resubmit:
    cmd = "crab -resubmit bad -c %s" % (ui_working_dir)
    os.system(cmd)

  if options.submit:
    cmd = "crab -submit all -c %s" % (ui_working_dir)
    os.system(cmd) 

  if options.kill:
    cmd = "crab -kill all -c %s" % (ui_working_dir)
    os.system(cmd) 
