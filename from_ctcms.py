# Copyright (c) Kamal Choudhary
import os
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element #,all_symbols
import sys
import urllib2
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

from pymatgen.core.composition import Composition
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter

from mpinterfaces import get_struct_from_mp
from mpinterfaces.MP_lammps import CalibrateLammps
from mpinterfaces.utils import *

# all the info/warnings/outputs redirected to the log file: 
# lammps_Al_O.log
logger = get_logger('lammps_Al_O')




#els=all_symbols()
els=Element #['Al']
allsymb=Element #all_symbols()
files=[]
for el in els:
    try:
       #os.makedirs(str(os.getcwd())+'/'+str(el)+str('.dir'))
       #os.chdir(str(os.getcwd())+'/'+str(el)+str('.dir'))
       
       #change directory name here
       os.makedirs(str('/scratch/lfs/kamal/JARVIS')+'/'+str(el)+str('.dir')) 
       os.chdir(str('/scratch/lfs/kamal/JARVIS')+'/'+str(el)+str('.dir'))
       print(os.getcwd() )
       link=str("http://www.ctcms.nist.gov/potentials/")+str(el)+str(".html")
       page = urllib2.urlopen(link).read()
       soup = BeautifulSoup(page)
       soup.prettify()
       for anchor in soup.findAll('a', href=True):
           li=anchor['href'].strip()
           if  li.startswith("./Download/"):
               handle= anchor['href']
               handle= handle.replace("./","")
               
               if handle.endswith(".alloy"):
                  down_link=str("http://www.ctcms.nist.gov/potentials/")+str(handle)
                  list_el=handle.replace("Download/","").split("-")
                  arr_alloy=[el]
                  structures=[]
                  for e in list_el:
                      if  e in allsymb and e not in arr_alloy:
                          arr_alloy.append(e)
                  print el,down_link,"arr_al=",arr_alloy,"handle=",handle,"joined=",'-'.join(arr_alloy)
                  ff=anchor['href'].replace("./Download/","").split("/")[-1]
                  n_atom_types=len(arr_alloy)
                  #print el,down_link,arr_alloy,handle,'-'.join(arr_alloy)
                  alloy_file=urllib.URLopener()
                  alloy_file.retrieve(down_link,ff)
                  #os.makedirs(str('/scratch/lfs/kamal/JARVIS')+'/'+str(el)+str('.dir')+'/'+str(ff))
                  #os.chdir(str('/scratch/lfs/kamal/JARVIS')+'/'+str(el)+str('.dir')+'/'+str(ff))
                  print(os.getcwd() )

                  with MPRester(MAPI_KEY) as m:
                       data = m.get_entries_in_chemsys(arr_alloy,inc_structure='final', property_data=["material_id","icsd_id"])
                       structures=[]
                       icsd_arr=[]
                       mp_arr=[]
                       for d in data:
                           x=d.data["material_id"]
                           icsd=d.data["icsd_id"]
                           structure = m.get_structure_by_material_id(x)
                           structures.append(structure)
                           icsd_arr.append(icsd)
                           mp_arr.append(x)


                  scell_size = 1 #needs to be changed to 12 for COMB
                  for s in structures:
                      SpacegroupAnalyzer(s).get_conventional_standard_structure()
                      a, b, c = s.lattice.abc
                      s.make_supercell([ceil(scell_size/a)+1,
                                        ceil(scell_size/b)+1,
                                        ceil(scell_size/c)+1])





                  parameters = {'atom_style': 'charge',
                  'minimize':'1.0e-13  1.0e-20  1000  10000',
                  'fix':['fix_nve all nve',
                     '1 all box/relax aniso 0.0 vmax 0.001'] }
                  # list of pair styles
                  #pair_coeff_files = []
                  pair_styles = ['eam/alloy']
                  #pair_coeff_files.append([os.path.join(os.getcwd(), ff)])
                  # list of pair coefficient files
                  pair_coeff_files = [os.path.join(os.getcwd(), ff)]
                  # this file must be in the folder where this script is run
                  pair_coeff_files = [os.path.join(os.getcwd(), ff)]
                  turn_knobs = OrderedDict([('STRUCTURES', structures),('PAIR_STYLE', pair_styles),('PAIR_COEFF', pair_coeff_files)] )
                  # job directory and run settings
                  job_dir = str('lammps')+str('_')+str(ff)
                  nprocs = 4
                  nnodes = 1
                  walltime = '01:15:00'
                  mem = 500
                  job_bin = '/home/km468/Software/lammps/src/lmp_ufhpc < inp'
                  qadapter, job_cmd = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,walltime=walltime,job_bin=job_bin, mem=mem)
                  checkpoint_files = []
                  file=str(ff)+'_data.json'
                  f=open(file,'w')
                  f.write(json.dumps([dict(mp_id=mp) for mp in mp_arr]))
                  f.close()
                 
                 # f=open(mp_file,'w')
                 # for mp in  mp_arr:
                 #     f.write(json.dumps([dict(mp_id=mp)]))
                 # f.close()
                  checkpoint_files = []
                  chkpt_file = str(ff)+'_step1.json'


                  # setup calibration jobs and run
                  cal = CalibrateLammps(parameters, turn_knobs=turn_knobs,
                          qadapter=qadapter, job_cmd=job_cmd,
                          job_dir=job_dir, is_matrix=True,
                          checkpoint_file=chkpt_file,
                          cal_logger=logger)
                  cal.setup()
                  cal.run()
                  checkpoint_files.append(chkpt_file)
                  #os.chdir('../')





               elif handle.endswith(".fs"):
                  down_link=str("http://www.ctcms.nist.gov/potentials/")+str(handle)
                  list_el=handle.replace("Download/","").split("-")
                  arr_fs=[el]
                  for e in list_el:
                      if  e in allsymb and e not in arr_fs:
                          arr_fs.append(e)
##                  print el,down_link,arr_fs
               elif handle.endswith(".set"):
                  down_link=str("http://www.ctcms.nist.gov/potentials/")+str(handle)
                  list_el=handle.replace("Download/","").split("-")
                  arr_set=[el]
                  for e in list_el:
                      if  e in allsymb and e not in arr_set:
                          arr_set.append(e)
##                  print el,down_link,arr_set
#               elif handle.endswith(".dat"):
#                  print ""
#               elif handle.endswith(".pdf"):
#                  print ""
#               elif handle.endswith(".gz"):
#                  print ""
#               elif handle.endswith(".plt"):
#                  print ""
#               else:
#                  files.append(handle)
       #os.chdir("../")
       os.chdir('../')
    except:
       #print el,"None"
       pass

#testfile=urllib.URLopener()
#testfile.retrieve("http://www.ctcms.nist.gov/potentials/Download/U-Mo-Xe-SKS13/U_Mo_Xe.2013.eam.alloy","ff")

#print files
