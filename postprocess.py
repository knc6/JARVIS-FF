
import os
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
import sys,traceback
import urllib2
import zipfile 
import urllib
from BeautifulSoup import BeautifulSoup
import json
MAPI_KEY = os.environ.get("MAPI_KEY", "")
from pymatgen.matproj.rest import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from math import ceil
from collections import OrderedDict

import matplotlib
matplotlib.use('Agg')
import glob,json
from pymatgen.core.composition import Composition
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter

from mpinterfaces import get_struct_from_mp
from mpinterfaces.lammps import CalibrateLammps
from mpinterfaces.utils import *
import os, zipfile  ,shutil
from os.path import join
import httplib
import json

base = "/jar"
headers = {"Accept": "application/json"}
def ZipDir(inputDir, outputZip):
    '''Zip up a directory and preserve symlinks and empty directories'''
    zipOut = zipfile.ZipFile(outputZip, 'w', compression=zipfile.ZIP_DEFLATED)

    rootLen = len(os.path.dirname(inputDir))
    def _ArchiveDirectory(parentDirectory):
        #contents = os.listdir(parentDirectory)
        contents = ['init.mod', 'potential.mod', 'in.elastic', 'data',  'log.lammps', 'restart.equil' ]
        #store empty directories
        if not contents:
            #http://www.velocityreviews.com/forums/t318840-add-empty-directory-using-zipfile.html
            archiveRoot = parentDirectory[rootLen:].replace('\\', '/').lstrip('/')
            zipInfo = zipfile.ZipInfo(archiveRoot+'/')
            zipOut.writestr(zipInfo, '')
        for item in contents:
            fullPath = os.path.join(parentDirectory, item)
            if os.path.isdir(fullPath) and not os.path.islink(fullPath):
                _ArchiveDirectory(fullPath)
            else: 
                archiveRoot = fullPath[rootLen:].replace('\\', '/').lstrip('/')
                if os.path.islink(fullPath):
                    # http://www.mail-archive.com/python-list@python.org/msg34223.html
                    zipInfo = zipfile.ZipInfo(archiveRoot)
                    zipInfo.create_system = 3
                    # long type of hex val of '0xA1ED0000L',
                    # say, symlink attr magic...
                    zipInfo.external_attr = 2716663808L
                    zipOut.writestr(zipInfo, os.readlink(fullPath))
                else:
                    zipOut.write(fullPath, archiveRoot, zipfile.ZIP_DEFLATED)
    _ArchiveDirectory(inputDir)

    zipOut.close()
        
    
    

def push_case(data={}):
    conn = httplib.HTTPConnection("", 5900)
    conn.request("POST","%s/push/case"%base, json.dumps(data), headers)
    response = conn.getresponse()
    _data = response.read()
    conn.close()
    return _data



def analyz_loge(log='log.lammps'):
    import sys
    """ Analyzes log.lammps file, at present  energy/atom data extraction is implemented"""
    en=0
    press=0
    c11=0
    c22=0
    c33=0
    c44=0
    c55=0
    c66=0
    c12=0
    c13=0
    c23=0
    c14=0
    c15=0
    c16=0
    c14=0
    c24=0
    c25=0
    c26=0
    c34=0
    c35=0
    c36=0
    c45=0
    c46=0
    c56=0
    try:
       logfile = open(log, "r")
       lines = logfile.read().splitlines()
       for i, line in enumerate(lines):
           if 'Loop time of' in line:
               toten=float(lines[i-1].split()[12])
               press=float(lines[i-1].split()[2])
               press=float(press)*0.0001
               en=float(lines[i-1].split()[12])/float(lines[i-1].split()[17])
               break
       logfile.close()
    except:
       pass
    try:
       logfile = open(log, "r")
       lines = logfile.read().splitlines()
       for i, line in enumerate(lines):
           if 'print "Elastic Constant C11all = ${C11all} ${cunits}"' in line:
               c11=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C22all = ${C22all} ${cunits}"' in line:
               c22=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C33all = ${C33all} ${cunits}"' in line:
               c33=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C12all = ${C12all} ${cunits}"' in line:
               c12=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C13all = ${C13all} ${cunits}"' in line:
               c13=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C23all = ${C23all} ${cunits}"' in line:
               c23=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C44all = ${C44all} ${cunits}"' in line:
               c44=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C55all = ${C55all} ${cunits}"' in line:
               c55=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C66all = ${C66all} ${cunits}"' in line:
               c66=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C14all = ${C14all} ${cunits}"' in line:
               c14=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C15all = ${C15all} ${cunits}"' in line:
               c15=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C16all = ${C16all} ${cunits}"' in line:
               c16=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C24all = ${C24all} ${cunits}"' in line:
               c24=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C25all = ${C25all} ${cunits}"' in line:
               c25=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C26all = ${C26all} ${cunits}"' in line:
               c26=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C34all = ${C34all} ${cunits}"' in line:
               c34=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C35all = ${C35all} ${cunits}"' in line:
               c35=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C36all = ${C36all} ${cunits}"' in line:
               c36=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C45all = ${C45all} ${cunits}"' in line:
               c45=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C46all = ${C46all} ${cunits}"' in line:
               c46=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
           if 'print "Elastic Constant C56all = ${C56all} ${cunits}"' in line:
               c56=((str((lines[i+1])).split("=")[1]).split("GPa"))[0]
       logfile.close()
    except:
       pass
    return round(en,2),round(press,2),float(toten),round(float(c11),1),round(float(c22),1),round(float(c33),1),round(float(c12),1),round(float(c13),1),round(float(c23),1),round(float(c44),1),round(float(c55),1),round(float(c66),1),round(float(c14),1),round(float(c15),1),round(float(c16),1),round(float(c24),1),round(float(c25),1),round(float(c26),1),round(float(c34),1),round(float(c35),1),round(float(c36),1),round(float(c45),1),round(float(c46),1),round(float(c56),1)

#print en,c11,c22

# all the info/warnings/outputs redirected to the log file: 
# lammps_Al_O.log
logger = get_logger('lammps_Al_O')

#def analyz_log(log='log.lammps'):
#    """ Analyzes log.lammps file, at present  energy/atom data extraction is implemented"""
#    logfile = open(log, "r")
#    lines = logfile.read().splitlines()
#    try:
#       for i, line in enumerate(lines):
#           if 'Loop time of' in line:
#               en=float(lines[i-1].split()[12])/float(lines[i-1].split()[17])
#       logfile.close()
#    except:
#       en=0
#    return en


els=Element#all_symbols()
#els=['Ti']
allsymb=Element #all_symbols()
files=[]
info=[]
infoarr=[]
memory=[]
data_mem={}
id=0
_cases = []
       #for nist in glob.glob('lammps_Mishin-Ni-Al-Co-2013.eam.alloy.nist'):
for nist in glob.glob('*nist'):
    #for el in els:
        #try:
           ff=((nist.split("lammps_")[1]).split(".nist")[0])
           mp_json=str(ff)+str("_data.json")
           step1_json=str(ff)+str("_step1.json")
           all_jobs = jobs_from_file(step1_json)
           file=mp_json
           with open(file, 'r') as f:
                data = json.load(f)
           for mp_id,job in zip(data,all_jobs):
               try:
                   log=str(job.job_dir)+'/log.lammps'
               #en=analyz_log(log=log)
                   (en,press,toten,c11,c22,c33,c12,c13,c23,c44,c55,c66,c14,c15,c16,c24,c25,c26,c34,c35,c36,c45,c46,c56)=analyz_loge(log)
                   try:  
                       if job.vis.mplmp.parameters['units'] :
                          en=float(en)*0.043 
                          toten=float(toten)*0.043 
                          #c11=float(c11)*0.000101325 
                          #c22=float(c22)*0.000101325 
                          #c33=float(c33)*0.000101325 
                          #c12=float(c12)*0.000101325 
                          #c13=float(c13)*0.000101325 
                          #c23=float(c23)*0.000101325 
                          #c44=float(c44)*0.000101325 
                          #c55=float(c55)*0.000101325 
                          #c66=float(c66)*0.000101325 
                          #c14=float(c14)*0.000101325 
                          #c16=float(c16)*0.000101325 
                          #c24=float(c24)*0.000101325 
                          #c25=float(c25)*0.000101325 
                          #c26=float(c26)*0.000101325 
                          #c34=float(c34)*0.000101325 
                          #c35=float(c35)*0.000101325 
                          #c36=float(c36)*0.000101325 
                          #c45=float(c45)*0.000101325 
                          #c46=float(c46)*0.000101325 
                          #c56=float(c56)*0.000101325 
                   except:
                          pass
                   set_el=(job.vis.mplmp.structure.symbol_set)
                   el_list=sorted(list(set_el)) 
                   #el_list=sorted(set_el,key=lambda student: student[0]) 
                   search= ('-'.join([item for item in el_list]))
                   #search=sorted(search.split("-"))
                   #search=search.join("-")

           #print mp_id['mp_id'],'   ',(job.vis.mplmp.structure.symbol_set),'    ',search,'    ',en,'   ',ff
               #info=(search, mp_id['mp_id'],mp_id['pf'],en,ff,c11,c22,c33,c12,c13,c23,c44,c66)
                   id=id+1
                   case=str("Calc-")+str(id)
                   Bv=float(c11+c22+c33)/9.0+float(c12+c13+c23)*2.0/9.0
                   Gv=(float(c11+c22+c33)-float(c12+c13+c23))/15.0+float(c44+c55+c66)/5.0
                   Bv=round(float(Bv),1)
                   Gv=round(float(Gv),1)
                   info=(search, case,mp_id['mp_id'],(str(''.join(['%s%s' % (k,int(v)) for k,v in (mp_id['ucf']).iteritems()]))),en,c11,c22,c33,c12,c13,c23,c44,c55,c66,Bv,Gv,round(mp_id['ehull'],2),ff)
   
                   target_file=str("/scratch/lfs/kamal/JARVIS/All2/COMBINED/")+str(job.job_dir)
                   dest=str("/scratch/lfs/kamal/JARVIS/All2/STORE/")+str(case)+str(".zip")
                   ZipDir(target_file,dest)
                   f = zipfile.ZipFile(str("/scratch/lfs/kamal/JARVIS/All2/STORE/")+str(ff)+str(".zip"), 'w')

                   f.write(str(ff))
                   print case
                   _case = {}
                   _case['case-number'] = case
                   _case['search'] = str(search)
                   _case['structure'] = (str(''.join(['%s%s' % (k,int(v)) for k,v in (mp_id['ucf']).iteritems()])))
                   _case['energy'] = str(en)
                   _case['ehull'] = str(mp_id['ehull'])
                   _case['totenergy'] = str(toten)
                   _case['composition'] = str(job.vis.mplmp.structure.composition)
#                   print type( job.vis.mplmp.structure.composition),job.vis.mplmp.structure.composition.to_json()

                   _case['forcefield'] = ff
                   _case['mpid'] = mp_id['mp_id']
                   _case['elastic_matrix'] = {}
                   _case['elastic_matrix']['c1'] = [0 for i in range(6)]
                   _case['elastic_matrix']['c2']	= [0 for i in range(6)]
                   _case['elastic_matrix']['c3']	= [0 for i in range(6)]
                   _case['elastic_matrix']['c4']	= [0 for i in range(6)]
                   _case['elastic_matrix']['c5']	= [0 for i in range(6)]
                   _case['elastic_matrix']['c6']	= [0 for i in range(6)]
                   

                   _case['elastic_matrix']['c1'][0] = c11
                   _case['elastic_matrix']['c1'][1] = c12
                   _case['elastic_matrix']['c1'][2] = c13
                   _case['elastic_matrix']['c1'][3] = c14
                   _case['elastic_matrix']['c1'][4] = c15
                   _case['elastic_matrix']['c1'][5] = c16
                   _case['elastic_matrix']['c2'][0] = c12
                   _case['elastic_matrix']['c2'][1] = c22
                   _case['elastic_matrix']['c2'][2] = c23
                   _case['elastic_matrix']['c2'][3] = c24
                   _case['elastic_matrix']['c2'][4] = c25
                   _case['elastic_matrix']['c2'][5] = c26
                   _case['elastic_matrix']['c3'][0] = c13
                   _case['elastic_matrix']['c3'][1] = c23
                   _case['elastic_matrix']['c3'][2] = c33
                   _case['elastic_matrix']['c3'][3] = c34
                   _case['elastic_matrix']['c3'][4] = c35
                   _case['elastic_matrix']['c3'][5] = c36
                   _case['elastic_matrix']['c4'][0] = c14
                   _case['elastic_matrix']['c4'][1] = c24
                   _case['elastic_matrix']['c4'][2] = c34
                   _case['elastic_matrix']['c4'][3] = c44
                   _case['elastic_matrix']['c4'][4] = c45
                   _case['elastic_matrix']['c4'][5] = c46
                   _case['elastic_matrix']['c5'][0] = c15
                   _case['elastic_matrix']['c5'][1] = c25
                   _case['elastic_matrix']['c5'][2] = c35
                   _case['elastic_matrix']['c5'][3] = c45
                   _case['elastic_matrix']['c5'][4] = c55
                   _case['elastic_matrix']['c5'][5] = c56
                   _case['elastic_matrix']['c6'][0] = c16
                   _case['elastic_matrix']['c6'][1] = c26
                   _case['elastic_matrix']['c6'][2] = c36
                   _case['elastic_matrix']['c6'][3] = c46
                   _case['elastic_matrix']['c6'][4] = c56
                   _case['elastic_matrix']['c6'][5] = c66

                   _case['Bv'] = str(Bv)
                   _case['Gv'] = str(Gv)
                   _cases.append(_case)
                   #data_mem[case]={"search":search,"mpid":mp_id['mp_id'],"formula":job.vis.mplmp.structure.formula,"ff":ff,"en":en,"c11":c11,"c22":c22,"c33":c33,"c12":c12,"c13":c13,"c23":c23,"c44":c44,"c55":c55,"c66":c66,"Bv":Bv,"Gv":Gv,"ehull":mp_id['ehull']}
                   infoarr.append(info)
                   #print push_case(_case)
               except:
                      #import traceback
                      #print traceback.print_exc()
                      #print job.job_dir


                      pass
        #except:
        #      pass
Intfarray=[]
Intfarray=sorted(infoarr,key=lambda student: student[0])
#print Intfarray[27][4]
count=0
keys=[]
f=open('EAM_alloy16.csv','w')
line= str('key')+str(',')+str('symbol')+str(',')+str('period')+str(',')+str('group')+'\n'
f.write(line)
new=0
for i in range(0,len(Intfarray)):
    if Intfarray[i][0] not in keys:
        keys.append(Intfarray[i][0])
        new=new+1
        count=0
    elif Intfarray[i][0]  in keys:
        count=count+1
    for j in range(0,len(Intfarray[0])):
        #print Intfarray[i][j],count+1,j+1
        line= str(new)+str(',')+str(Intfarray[i][j])+str(',')+str(count+1)+str(',')+str(j+1)+'\n'
        f.write(line)
print count,new,id
search_string='Case-209'
for element,content in  data_mem.iteritems():
    if search_string==element:
       print content['ff']
#       print json.dumps(content, sort_keys=True,indent=4, separators=(',', ': '))
#
#
json_file='data.json'
f_json=open(json_file,'w')
#f_json.write(json.dumps(data_mem))
f_json.write(json.dumps(_cases))
f_json.close()




