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
(options, args) = parser.parse_args()

if options.run:
  options.create_cfg = True

datasets = [
    # Multijet analysis
#    ["/QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-1000to1400_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-1000to1400_pythia"],
#    ["/QCD_Pt-1000_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-1000toInf_herwig_START53_V7C_22Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-1000toInf_herwig"],
#    ["/QCD_Pt-120to170_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-120to170_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-120to170_pythia"],
#    ["/QCD_Pt_120to170_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-120to170_herwig_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-120to170_herwig"],
#    ["/QCD_Pt-300to470_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-300to470_herwig_START53_V7C_22Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-300to470_herwig"],
#    ["/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-300to470_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-300to470_pythia"],
#    ["/QCD_Pt-170to300_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-170to300_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-170to300_pythia"],
#    ["/QCD_Pt_170to300_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-170to300_herwig_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-170to300_herwig"],
#    ["/QCD_Pt-470to600_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-470to600_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-470to600_pythia"],
#    ["/QCD_Pt-600to800_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-600to800_herwig_START53_V7C_22Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-600to800_herwig"],
#    ["/QCD_Pt-600to800_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-600to800_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-600to800_pythia"],
#    ["/QCD_Pt-470to600_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-470to600_herwig_START53_V7C_22Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-470to600_herwig"],
#    ["/QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-800to1000_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-800to1000_pythia"],
#    ["/QCD_Pt-800to1000_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-800to1000_herwig_START53_V7C_22Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-800to1000_herwig"],
    ["/QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6/apequegn-QCD_HT-1000ToInf_START53_V7A_04Dec13-v1-869f7bc494a6202e5b9ae91c5af1c97d/USER", "QCD_HT-1000ToInf"],
    ["/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6/apequegn-QCD_HT-500To1000_START53_V7A_04Dec13-v1-869f7bc494a6202e5b9ae91c5af1c97d/USER ", "QCD_HT-500To1000"]
    
    ]

# Get email address
email = "%s@ipnl.in2p3.fr" % (pwd.getpwuid(os.getuid()).pw_name)

d = datetime.datetime.now().strftime("%d%b%y")

version = 1

print("Creating configs for crab. Today is %s, you are %s and it's version %d" % (d, email, version))
print("")

for dataset in datasets:

  dataset_name = dataset[1] + "_woPU_pt40_eta50"
  dataset_path = dataset[0]
  dataset_size = -1
  if len(dataset) > 2:
    dataset_size = dataset[2]
  #dataset_globaltag = re.search('START\d{0,2}_V\d[A-Z]?', dataset_path).group(0)

  #publish_name = "%s_%s_%s-v%d" % (dataset_name, dataset_globaltag, d, version)
  output_file = "crab_MC_%s_%s.cfg" % (dataset_name, d)
  ui_working_dir = ("crab_MC_%s_%s") % (dataset_name, d)

  if options.create_cfg:
    output_dir = ("Extracted_step2/MC/%s/%s" % (d, dataset_name))
  
  python_config = "Extractor_MULTIJET_MC.py";

  print("Creating config file for '%s'" % (dataset_path))
  print("\tName: %s" % dataset_name)
  print("\tOutput directory: %s" % output_dir)
  print("")

  os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@email@#%s#g\" -e \"s#@pset@#%s#g\" -e \"s#@outputdir@#%s#\" crab_MC.cfg.template.ipnl > %s" % (dataset_path, ui_working_dir, email, python_config, output_dir, output_file))

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
