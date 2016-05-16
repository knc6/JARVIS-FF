# JARVIS-FF

JARVIS-FF contains the scripts for making repository for evaluating classical force field available at http://www.ctcms.nist.gov/~knc6/periodic.html.
1) from_ctcms.py is the script to get force-fields from ctcms interatomic potential repository website (http://www.ctcms.nist.gov/potentials/)
2) ctcms_alloy.py is the scripts for setting up MPinterface (https://github.com/JARVIS-Unifies/MPInterfaces) and  running LAMMPS calculation. Similar scripts have been used for other force-fields such as AIREBO,COMB, ReaxFF and so on.
2) postprocess.py is the post-procesing script after running ctcms_alloy.py
3) data.json is the classical force-field data for matrerial properties such as elastic constnats, energies.
4) get_entries_from_json.py and get_entries_from_ecjson.py is the script to compare materials project (MP) and JARVIS-FF data such as energetics, convex hull plot. 
Note lammps.py in MPinterfaces has been changed to MP_lammps.py
