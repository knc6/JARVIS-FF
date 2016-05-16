# Copyright (c) Kamal Choudhary
import matplotlib
matplotlib.use('Agg')


import json
from pymatgen.core.composition import Composition
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter

##data = json.loads(open('data.json').read())
count=0
def get_entries_from_json(key=None,value=None,out=None,ff=None,json_f=None,PD=False):
    output=[[]]
 
    file=json_f
    with open(file, 'r') as f:
         data = json.load(f)
    if key=='search':
       output=[]
       try:
          old=value.split('-')
          el_list=sorted(old)
          value=('-'.join([item for item in el_list]))
       except:
          pass
    for d in data:
            x = {}
            options=['mpid','ehull','Bv','Gv','totenergy','energy','case-number','elastic_matrix','forcefield','composition']
            x[key] = str(d[key])
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
            #if key  != 'search':
            if x[key] ==value and  x['forcefield']==ff and PD == False:
                 output= x[out]
            if x[key] ==value and  ff==None and PD == False:
                 print (x[out],x['forcefield'])
            if x[key] ==value and  x['forcefield']==ff and PD==True:
                    output.append(PDEntry(Composition(x['composition']),float(x['totenergy'])))
            if key=='search'  and len(value.split('-')) >1 and PD ==True :
               els=value.split('-')
               for el in els:
                  if x[key] ==el and  x['forcefield']==ff:
                      output.append(PDEntry(Composition(x['composition']),float(x['totenergy'])))  
               for el1 in els:
                   for el2 in els:
                      if x[key] ==str(el1)+str('-')+str(el2) and  x['forcefield']==ff:
                          output.append(PDEntry(Composition(x['composition']),float(x['totenergy'])))  
                      if x[key] ==str(el2)+str('-')+str(el1) and  x['forcefield']==ff:
                          output.append(PDEntry(Composition(x['composition']),float(x['totenergy'])))  
               
    return output
                       

##Used for phase diagram
entries=get_entries_from_json(key='search',value='Cu-H-O',out='energy',ff='ffield.CuOCH.comb3',json_f='data.json',PD=True)
#entries=get_entries_from_json(key='search',value='Fe-H-O',out='energy',ff='ffield.reax.Fe_O_C_H',json_f='data.json',PD=True)
#entries=get_entries_from_json(key='search',value='Al-Ni',out='energy',ff='Mishin-Ni-Al-2009.eam.alloy',json_f='data.json',PD=True)
pd = PhaseDiagram(entries)
plotter = PDPlotter(pd, show_unstable=False)
#plotter = PDPlotter(pd, show_unstable=True)
#print plotter.pd_plot_data
print plotter.get_plot
name=str('Cu-H-Ocomb')+str('_LMP_phasediagram.png')
plotter.write_image(name,image_format="png")
##

##INDIVIDUAL
#output=get_entries_from_json(key='mpid',value='mp-13',out='energy',ff='ffield.reax.Fe_O_C_H',json_f='data.json')
get_entries_from_json(key='mpid',value='mp-134',out='Bv',json_f='data.json')
##
