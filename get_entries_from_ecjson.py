# Copyright (c) Kamal Choudhary
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import json
from pymatgen.core.composition import Composition
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter
import numpy as np
##data = json.loads(open('data.json').read())
count=0
ref=open('RESULTS','w')
ref1=open('RESULTS1','w')
mp_arr=[]
ff_arr=[]
def get_entries_from_json(key=None,value=None,out=None,ff=None,json_f=None,PD=False):
    output=[[]]
 
    file=json_f #FF data
    with open(file, 'r') as f:
         data = json.load(f)
    with open('ec.json', 'r') as f1:
         data1 = json.load(f1)
    for c in data1:
        y = {}
        y['material_id'] = str(c['material_id'])
        y['K_Voigt'] = str(c['K_Voigt'])
        y['G_Voigt'] = str(c['G_Voigt'])
        for d in data:
                x = {}
                options=['mpid','ehull','Bv','Gv','totenergy','energy','case-number','elastic_matrix','forcefield','composition']
                x['mpid'] = str(d['mpid'])
                x['ehull'] = str(d['ehull'])
                x['Bv'] = str(d['Bv'])
                x['Gv'] = str(d['Gv'])
                x['totenergy'] = str(d['totenergy'])
                x['energy'] = str(d['energy'])
                x['case-number'] = str(d['case-number'])
                x['elastic_matrix'] = str(d['elastic_matrix'])
                x['forcefield'] = str(d['forcefield'])
                x['composition'] = str(d['composition'])
                if x['mpid'] ==y['material_id'] and  ff==None and PD == False:
                     print (x['mpid'],x['Bv'],y['K_Voigt'],x['Gv'],y['G_Voigt'],x['forcefield'],x['ehull'])
                     line= str(x['mpid'])+' '+str(x['Bv'])+ ' '+str(y['K_Voigt'])+' '+str(x['Gv'])+' '+str(y['G_Voigt'])+' '+str(x['forcefield'])+str(x['ehull'])+'\n'
                     ref.write(line)
                     if (float(x['Bv']))>0.0:
                         ff_arr.append(float(x['Bv']))
                         mp_arr.append(float(y['K_Voigt']))
                     if float(x['ehull'])==0.0:
                        line= str(x['mpid'])+' '+str(x['Bv'])+ ' '+str(y['K_Voigt'])+' '+str(x['Gv'])+' '+str(y['G_Voigt'])+' '+str(x['forcefield'])+str(x['ehull'])+'\n'
                        ref1.write(line)
    plt.plot(mp_arr, ff_arr, 'bo')           
    plt.xlabel('Bv-MP (GPa)')
    plt.ylabel('Bv-FF (GPa)')
    t1 = np.arange(0.0, 800.0, 100)
    plt.plot(t1,t1,'g-')
    plt.xlim(0,700)
    plt.ylim(0,700)
    plt.savefig('Bv.svg')
##PD

##INDIVIDUAL
#output=get_entries_from_json(key='mpid',value='mp-13',out='energy',ff='ffield.reax.Fe_O_C_H',json_f='data.json')
get_entries_from_json(out='Bv',json_f='data.json')
##
