# JARVIS-FF

JARVIS-FF contains the scripts for making repository for evaluating classi cal force field available at http://www.ctcms.nist.gov/~knc6/periodic.html.
ctcms_alloy.py is the scripts for setting up MPinterface for running LAMMPS calculation. Similar scripts have been used for other force-fields such as AIREBO,COMB, ReaxFF and so on.
data.json is the classical force-field data for matrerial properties such as elastic constnats, energies.
from_ctcms.py is the script to get force-fields from ctcms website.
get_entries_from_json.py and get_entries_from_ecjson.py is the script to compare materials project and JARVIS-FF data such as energetics, convex hull plot.
postprocess.py is the post-procesing script after running ctcms_alloy.py
Note lammps.py in MPinterfaces has been changed to MP_lammps.py
