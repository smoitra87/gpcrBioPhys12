""" 

Deal with issues regarding gremlin run
"""

import glob,os

#----------------------------------------------------------------------
# Find missing results and rerun those experiments

# Get the list of poutputb files

done_files = glob.glob("out_perm*")
range_res = (1,348)

pre_str = 'out_perm_GPCR_ranga_'
offset = len(pre_str)

p = map(lambda x : x[offset:],map(lambda x :os.path.splitext(x)[0],\
	done_files))
p = sorted(map(lambda x : int(x),p))
set_done = set(p)
set_notdone = set(range(range_res[0],range_res[1]+1)) - set_done
for id in set_notdone : 
	print id
