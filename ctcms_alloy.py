# Copyright (c) Kamal Choudhary
import os
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
import sys,shutil
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
import glob


for file in glob.glob("*.alloy"):
    try:
        ff=str(file)
        element_ff=[]
        f=open(ff,'r')
        list_el=[]
        lines=f.readlines()
        content=(lines[3]).split(" ")
    #content=(lines[3]).split("' '|\n|\r\n")
        for val in content:
        
            if val != '' and val !='\n' and val !='\r\n':
               list_el.append(val)
        for i in range(0,len(list_el)):
             if i!=0:
                 element_ff.append(list_el[i])    
#    print ff,' ',element_ff
        with MPRester(MAPI_KEY) as m:
             data = m.get_entries_in_chemsys(element_ff,inc_structure='final', property_data=["unit_cell_formula","material_id","icsd_id","spacegroup","energy_per_atom","formation_energy_per_atom","pretty_formula","band_gap","total_magnetization","e_above_hull"])
             if (len(element_ff)>1):
                 try:
                     entries = m.get_entries_in_chemsys(element_ff)
                     pd = PhaseDiagram(entries)
                     plotter = PDPlotter(pd, show_unstable=True)
                     image=str(ff)+str("_DFT")+str(".jpg")
                     plotter.write_image(image)
                 except:
                     pass
             structures=[]
             structures_cvn=[]
             icsd_arr=[]
             mp_arr=[]
             sg_arr=[]
             enp_arr=[]
             fenp_arr=[]
             pf_arr=[]
             ucf_arr=[]
             bg_arr=[]
             tm_arr=[]
             ehull_arr=[]
             for d in data:
                 x=d.data["material_id"]
                 sg=d.data["spacegroup"]
                 enp=d.data["energy_per_atom"]
                 fenp=d.data["formation_energy_per_atom"]
                 pf=d.data["pretty_formula"]
                 ucf=d.data["unit_cell_formula"]
                 bg=d.data["band_gap"]
                 tm=d.data["total_magnetization"]
                 ehull=d.data["e_above_hull"]
                 icsd=d.data["icsd_id"]
                 structure = m.get_structure_by_material_id(x)
                 structures.append(structure)
                 icsd_arr.append(icsd)
                 mp_arr.append(x)
                 sg_arr.append(sg)
                 enp_arr.append(enp)
                 fenp_arr.append(fenp)
                 pf_arr.append(pf)
                 bg_arr.append(bg)
                 tm_arr.append(tm)
                 ucf_arr.append(ucf)
                 ehull_arr.append(ehull)


                 scell_size = 1
             for s in structures:
                 svn=SpacegroupAnalyzer(s).get_conventional_standard_structure()
                 a, b, c = s.lattice.abc
                 svn.make_supercell([ceil(scell_size/a)+1,
                                   ceil(scell_size/b)+1,
                                   ceil(scell_size/c)+1])
                 structures_cvn.append(svn)
             parameters = {'atom_style': 'charge' ,'control_file':'/home/kamal/inelast.mod'}
             pair_styles = ['eam/alloy']
             pair_coeff_files = [os.path.join(os.getcwd(), ff)]
             file=str(ff)+'_data.json'
             f=open(file,'w')
             f.write(json.dumps([dict(ucf=ucf,mp_id=mp,icsd=icsd,sg=sg,enp=enp,fenp=fenp,pf=pf,bg=bg,tm=tm,ehull=ehull) for ucf,mp,icsd,sg,enp,fenp,pf,bg,tm,ehull in zip(ucf_arr,mp_arr,icsd_arr,sg_arr,enp_arr,fenp_arr,pf_arr,bg_arr,tm_arr,ehull_arr)]))
             f.close()
             turn_knobs = OrderedDict(
             [  
             ('STRUCTURES', structures_cvn),
             ('PAIR_STYLE', pair_styles),
             ('PAIR_COEFF', pair_coeff_files)
              ] )
    # job directory and run settings
    #job_dir = 'lammps_job'
             job_dir = str('lammps')+str('_')+str(ff)+str('.nist')
             nprocs = 4
             nnodes = 1
             walltime = '00:30:00'
             mem = 500
             job_bin = 'lmp_ufhpc.openmpi  < in.elastic'
             qadapter, job_cmd = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,
                                      walltime=walltime,
                                      job_bin=job_bin, mem=mem)
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
    #plotter = PDPlotter(pd, show_unstable=True)
             print ("ff=",ff)
         #    return checkpoint_files

             all_jobs = jobs_from_file(chkpt_file)
             entries = []
             for job in all_jobs:
                 comp = job.vis.mplmp.structure.composition
                 energy =job.final_energy
                 entries.append(PDEntry(comp,energy))
             try:
                pd = PhaseDiagram(entries)
                plotter = PDPlotter(pd, show_unstable=True)
                image_ff=str(ff)+str("_CLASS")+str(".jpg")
                plotter.write_image(image_ff)
             except:
                  pass
  
    except:
          pass     ##        pass
