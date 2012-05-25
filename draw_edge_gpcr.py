"""
Reads from edge files to draw edges on the GPCR structure

"""

import os,sys

import pymol
pymol.finish_launching()
from pymol import cmd

import glob;


DATA_DIR='data'
edgef = 'common_GPCR_ranga_adrn_pa.dat'
edgepath = os.path.join(DATA_DIR,edgef);
outfp = os.path.join(DATA_DIR,os.path.splitext(edgef)[0]+'.png')


#----------------------------------------------------------------------
# Draw Edges

# Clear previous stuff
cmd.delete('all');

# Load pdb
pdbfile = os.path.join(DATA_DIR,'1F88.pdb');
cmd.load(pdbfile)

# prettify
cmd.hide('all')
cmd.show('ribbon','1F88 and chain A')
cmd.bg_color('white')
cmd.set_view (\
    '-0.737545371,   -0.547816157,   -0.394868791,\
     0.542158365,   -0.131744832,   -0.829878628,\
     0.402601302,   -0.826157331,    0.394175380,\
    -0.000016227,    0.000207499, -201.692947388,\
    52.060752869,   10.034209251,   -1.057026982,\
    73.561737061,  329.858551025,    0.000000000 ')

cmd.spectrum(palette = 'green_red', selection = '1F88 and chain A')
cmd.select('retinal','1F88 and het and chain A and resn ret')
cmd.show('spheres','retinal')
cmd.color('magenta','retinal')

# constants
cutoff = 2

# Load the edges from file
with open(edgepath,'r') as fin : 
    data = fin.readlines();
edges = [line.strip().split() for line in data]


for n,e in  enumerate(edges) : 
    t = cmd.dist('dist'+str(n),'1F88//A/'+e[0]+'/CA',\
		'1F88//A/'+e[1]+'/CA')

    if t < cutoff  : 
        cmd.delete('dist'+str(n)) 
    else : 
        cmd.color('blue','dist'+str(n))
        cmd.hide('labels','dist'+str(n))
        cmd.show('spheres','1F88//A/'+e[0]+'/CA')
        cmd.show('spheres','1F88//A/'+e[1]+'/CA')
        cmd.color('yellow','1F88//A/'+e[0]+'/CA')
        cmd.color('yellow','1F88//A/'+e[1]+'/CA')

cmd.set('dash_radius','0.2');
cmd.set('dash_gap','0.0');
cmd.set('sphere_scale','0.3')
cmd.set('sphere_scale','0.8','retinal');
cmd.ray()
print('*'*5+'Image:%s'%os.path.abspath(outfp)+'*'*5)
cmd.png(os.path.abspath(outfp))

print '-------Completed Call--------';

