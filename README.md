# JARVIS-FF

JARVIS-FF contains the scripts for making repository for evaluating classical force field available at http://www.ctcms.nist.gov/~knc6/periodic.html.
1) ctcms_alloy.py is the scripts for setting up MPinterface (https://github.com/JARVIS-Unifies/MPInterfaces) and  running LAMMPS calculation. Similar scripts have been used for other force-fields such as AIREBO,COMB, ReaxFF and so on.
2) data.json is the classical force-field data for matrerial properties such as elastic constnats, energies.
3) from_ctcms.py is the script to get force-fields from ctcms interatomic potential repository website (http://www.ctcms.nist.gov/potentials/).
4) get_entries_from_json.py and get_entries_from_ecjson.py is the script to compare materials project and JARVIS-FF data such as energetics, convex hull plot. 
5) postprocess.py is the post-procesing script after running ctcms_alloy.py
Note lammps.py in MPinterfaces has been changed to MP_lammps.py
