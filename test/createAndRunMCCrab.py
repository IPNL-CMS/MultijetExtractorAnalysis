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
    # Multijet analysis
#    ["/QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-1000to1400_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-1000to1400_pythia"],
    #["/QCD_Pt-1000_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-1000toInf_herwig_START53_V7C_22Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-1000toInf_herwig"],
    #["/QCD_Pt-120to170_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-120to170_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-120to170_pythia"],
    #["/QCD_Pt_120to170_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-120to170_herwig_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-120to170_herwig"],
    #["/QCD_Pt-300to470_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-300to470_herwig_START53_V7C_22Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-300to470_herwig"],
    #["/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-300to470_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-300to470_pythia"],
    #["/QCD_Pt-170to300_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-170to300_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-170to300_pythia"],
    #["/QCD_Pt_170to300_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-170to300_herwig_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-170to300_herwig"],
    #["/QCD_Pt-470to600_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-470to600_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-470to600_pythia"],
    #["/QCD_Pt-600to800_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-600to800_herwig_START53_V7C_22Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-600to800_herwig"],
    #["/QCD_Pt-600to800_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-600to800_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-600to800_pythia"],
    #["/QCD_Pt-470to600_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-470to600_herwig_START53_V7C_22Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-470to600_herwig"],
    #["/QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6/apequegn-QCD_Pt-800to1000_pythia_START53_V7A_20Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-800to1000_pythia"],
    #["/QCD_Pt-800to1000_CTEQ6L1_8TeV_herwig6/apequegn-QCD_Pt-800to1000_herwig_START53_V7C_22Dec13-v1-c5f9c59e100f883a59cec8d8908af608/USER", "QCD_Pt-800to1000_herwig"],
    #["/QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6/apequegn-QCD_HT-1000ToInf_START53_V7A_04Dec13-v1-869f7bc494a6202e5b9ae91c5af1c97d/USER", "QCD_HT-1000ToInf"],
    #["/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6/apequegn-QCD_HT-500To1000_START53_V7A_04Dec13-v1-869f7bc494a6202e5b9ae91c5af1c97d/USER ", "QCD_HT-500To1000"]

    #["/QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM", "QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_miniAOD"],
    #["/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISpring15DR74-Asympt50nsRaw_MCRUN2_74_V9A-v3/MINIAODSIM", "QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_bx50_miniAOD"],
    ["/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v3/MINIAODSIM", "QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_bx25_miniAOD"],
#    ["/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/apequegn-QCD_Pt-15to7000_Flat_13TeV_50ns_pythia8_25Jun15-v1-949f4a919651fc0956a5c1f2e3affdf2/USER", "QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_50ns_pythia8"],
    #["/QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6/apequegn-QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_25Jun15-v1-949f4a919651fc0956a5c1f2e3affdf2/USER", "QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6"],
    #["/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/apequegn-QCD_Pt-15to7000_Flat_13TeV_25ns_pythia8_25Jun15-v1-949f4a919651fc0956a5c1f2e3affdf2/USER", "QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_25ns_pythia8"],

    ["/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM", "QCD_HT100to200"],
    ["/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM", "QCD_HT200to300"],
    ["/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM", "QCD_HT300to500"],
    ["/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM", "QCD_HT500to700"],
    ["/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM", "QCD_HT700to1000"],
    ["/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM", "QCD_HT1000to1500"],
    ["/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM", "QCD_HT1500to2000"],
    ["/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM", "QCD_HT2000toInf"],
    
    ]

# Get email address
email = "%s@ipnl.in2p3.fr" % (pwd.getpwuid(os.getuid()).pw_name)

d = datetime.datetime.now().strftime("%d%b%y")

version = 1

print("Creating configs for crab. Today is %s, you are %s and it's version %d" % (d, email, version))
print("")

#def processDataset(dataset):
for dataset in datasets:
    #dataset_name = dataset[1] + "_woPU_pt30_eta50_puJetIdT"
    dataset_name = dataset[1] + "_woPU_pt10_eta50_notRmPUJets"
    dataset_path = dataset[0]
    dataset_size = -1
    if len(dataset) > 2:
        dataset_size = dataset[2]
    #dataset_globaltag = re.search('START\d{0,2}_V\d[A-Z]?', dataset_path).group(0)

    #publish_name = "%s_%s_%s-v%d" % (dataset_name, dataset_globaltag, d, version)
    output_file = (os.path.join("tasks","crab_MC_%s_%s.py" % (dataset_name, d)))
    ui_working_dir = ("crab_MC_%s_%s") % (dataset_name, d)

    try:
        getUsernameFromSiteDB()
    except Exception as e:
        print("Error when trying to find CERN username")
        print e
        sys.exit
    #output_dir = ("/store/user/%s/Multijet/Extracted/MC/Run2/extractor_%s/%s/%s" % (getUsernameFromSiteDB(), getGitTag(), d, dataset_name))
    output_dir = ("/store/user/%s/Multijet/Extracted/MC/Run2/%s/%s" % (getUsernameFromSiteDB(), d, dataset_name))

    python_config = "Extractor_MULTIJET_MC.py";

    if options.submit:
        print("Creating config file for '%s'" % (dataset_path))
        print("\tName: %s" % dataset_name)
        print("\tOutput directory: %s" % output_dir)
        print("")

        if "/USER" in dataset_path:
            dbs_url = "phys03"
        else:
            dbs_url = "global"

        os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@remote_dir@#%s#g\" -e \"s#@email@#%s#g\" -e \"s#@dbs_url@#%s#g\" crab_MC_template.py > %s" % (dataset_path, ui_working_dir, output_dir, email, dbs_url, output_file))


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
