#! /usr/bin/env python

import os, copy, datetime, pwd, re, sys
from CRABClient.UserUtilities import getUsernameFromSiteDB

def check_output(*popenargs, **kwargs):
    import subprocess
    r"""Run command with arguments and return its output as a byte string.
 
    Backported from Python 2.7 as it's implemented as pure python on stdlib.
 
    >>> check_output(['/usr/bin/python', '--version'])
    Python 2.6.2
    """
    process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    if retcode:
        cmd = kwargs.get("args")
        if cmd is None:
            cmd = popenargs[0]
        error = subprocess.CalledProcessError(retcode, cmd)
        error.output = output
        raise error
    return output

def getGitTag():
    if "tag" not in getGitTag.__dict__:
        getGitTag.tag = check_output(["git", "describe", "--tags"]).rstrip('\n')

    return getGitTag.tag

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-j", "--process", action="store", dest="cores", type="int", default=1, help="Number of core to use for launching")
parser.add_option("", "--status", action="store_true", dest="status", default=False, help="run crab -status")
parser.add_option("", "--submit", action="store_true", dest="submit", default=False, help="run crab -submit all")
parser.add_option("", "--kill", action="store_true", dest="kill", default=False, help="kill crab jobs")
parser.add_option("", "--resubmit", action="store_true", dest="resubmit", default=False, help="Resubmit crab jobs. If no job ids are given, all (and only) jobs in status failed will be resubmitted")
parser.add_option("", "--get", action="store_true", dest="get", default=False, help="getlog crab jobs")
parser.add_option("", "--getlist", action="store_true", dest="getlist", default=False, help="getlist of crab jobs")
parser.add_option("", "--report", action="store_true", dest="report", default=False, help="report crab jobs")

(options, args) = parser.parse_args()

datasets = [
#["/Jet/apequegn-Jet_Run2012A-22Jan2013_09Apr14-v1-829b2d54640b0ff2d246edb621ac7ffd/USER" , "Jet_Run2012A-22Jan2013", "FT53_V21A_AN6"],
#["/JetHT/apequegn-JetHT_Run2012B-22Jan2013_09Apr14-v1-9f561dc7a7b870f3c507ce92c9c7f820/USER", "JetHT_Run2012B-22Jan2013", "FT53_V21A_AN6"],
#["/JetHT/apequegn-JetHT_Run2012C-22Jan2013_09Apr14-v1-9f561dc7a7b870f3c507ce92c9c7f820/USER", "JetHT_Run2012C-22Jan2013", "FT53_V21A_AN6"],
#["/JetHT/apequegn-JetHT_Run2012D-22Jan2013_09Apr14-v1-9f561dc7a7b870f3c507ce92c9c7f820/USER", "JetHT_Run2012D-22Jan2013", "FT53_V21A_AN6"],
#["/JetMon/apequegn-JetMon_Run2012B-22Jan2013_09Apr14-v1-3ba8c59f861a2218f0d0749eda03cb49/USER", "JetMon_Run2012B-22Jan2013", "FT53_V21A_AN6"],
#["/JetMon/apequegn-JetMon_Run2012C-22Jan2013_09Apr14-v1-3ba8c59f861a2218f0d0749eda03cb49/USER", "JetMon_Run2012C-22Jan2013", "FT53_V21A_AN6"],
#["/JetMon/apequegn-JetMon_Run2012D-22Jan2013_09Apr14-v1-3ba8c59f861a2218f0d0749eda03cb49/USER", "JetMon_Run2012D-22Jan2013", "FT53_V21A_AN6"]
#["/JetHT/apequegn-JetHT_7_4_X_RunD_25Jun15-v1-339c8c37e62a7df1890df9cb21e1059e/USER", "JetHT_7_4_X_RunD_25Jun15", "PHYS14_25_V2"]
#["/JetHT/Run2015B-PromptReco-v1/MINIAOD", "JetHT_Run2015B-PromptReco_miniAOD", "74X_dataRun2_Prompt_v0"]    
["/JetHT/Run2015C-PromptReco-v1/MINIAOD", "JetHT_Run2015C-PromptReco_miniAOD", "74X_dataRun2_v2"]    
    ]

# Get email address
email = "%s@ipnl.in2p3.fr" % (pwd.getpwuid(os.getuid()).pw_name)

d = datetime.datetime.now().strftime("%d%b%y")

version = 1

print("Creating configs for crab. Today is %s, you are %s and it's version %d" % (d, email, version))
print("")

#golden_json_file = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251642_13TeV_PromptReco_Collisions15_JSON.txt"
#golden_json_file = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt"
#golden_json_file = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-254349_13TeV_PromptReco_Collisions15_JSON_v2.txt"
golden_json_file = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt"

#def processDataset(dataset):
for dataset in datasets:
    #dataset_name = dataset[1] + "_woPU_pt30_eta50_puJetIdT"
    dataset_name = dataset[1] + "_woPU_pt10_eta50_notRmPUJets_usingGT"
    dataset_path = dataset[0]
    dataset_size = -1
    if len(dataset) > 2:
        dataset_size = dataset[2]
    #dataset_globaltag = re.search('START\d{0,2}_V\d[A-Z]?', dataset_path).group(0)
    dataset_globaltag = dataset[2]

    #publish_name = "%s_%s_%s-v%d" % (dataset_name, dataset_globaltag, d, version)
    output_file = (os.path.join("tasks","crab_data_%s_%s.py" % (dataset_name, d)))
    ui_working_dir = ("crab_data_%s_%s") % (dataset_name, d)

    try:
        getUsernameFromSiteDB()
    except Exception as e:
        print("Error when trying to find CERN username")
        print e
        sys.exit
    #output_dir = ("/store/user/%s/Multijet/Extracted/dC/Run2/extractor_%s/%s/%s" % (getUsernameFromSiteDB(), getGitTag(), d, dataset_name))
    output_dir = ("/store/user/%s/Multijet/Extracted/data/Run2/%s/%s" % (getUsernameFromSiteDB(), d, dataset_name))

    if options.submit:
        print("Creating config file for '%s'" % (dataset_path))
        print("\tName: %s" % dataset_name)
        print("\tGlobal tag: %s" % dataset_globaltag)
        print("\tOutput directory: %s" % output_dir)
        print("")

        if "/USER" in dataset_path:
            dbs_url = "phys03"
        else:
            dbs_url = "global"

        os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@remote_dir@#%s#g\" -e \"s#@email@#%s#g\" -e \"s#@globaltag@#%s#g\" -e \"s#@dbs_url@#%s#g\" -e \"s#@golden_json_file@#%s#g\" crab_data_template.py > %s" % (dataset_path, ui_working_dir, output_dir, email, dataset_globaltag, dbs_url, golden_json_file, output_file))


    if options.submit:
        cmd = "crab submit -c %s" % (output_file)
        #print cmd
        os.system(cmd)

    correct_ui_working_dir = "crab_%s" % ui_working_dir
    ui_working_dir_area = (os.path.join("tasks", correct_ui_working_dir)) 

    if options.status:
        cmd = "crab status -d %s" % ui_working_dir_area
        os.system(cmd)

    if options.resubmit:
        cmd = "crab resubmit -d %s" % ui_working_dir_area
        os.system(cmd)

    if options.get:
        cmd = "crab getlog -d %s" % ui_working_dir_area
        os.system(cmd)

    if options.report:
        cmd = "crab report -d %s" % ui_working_dir_area
        os.system(cmd)

    if options.kill:
        cmd = "crab kill -d %s" % ui_working_dir_area
        os.system(cmd)

    if options.getlist:
        cmd = "crab getoutput --xrootd -d %s" % ui_working_dir_area
        os.system(cmd)

#import multiprocessing
#pool = multiprocessing.Pool(options.cores)
#pool.map(processDataset, datasets)

