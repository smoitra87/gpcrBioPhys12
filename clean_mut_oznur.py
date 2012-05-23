"""
Cleans oznur file
""" 
ozf = "data/mut_oznur_raw.dat"
ozf2 = "data/mut_oznur.dat"
with open(ozf) as fin :
	 lines = fin.readlines();

lines = filter(lambda x : len(x) > 0,map(str.strip,lines))

with open(ozf2,'w') as fout :
	fout.write('\n'.join(lines))


