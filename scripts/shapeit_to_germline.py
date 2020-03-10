# remake shapeit -> germline script

from numpy import *
from scipy.interpolate import interp1d
from sys import argv

# use, e.g.:
# python shapeit_to_germline.py shapeitoutputprefix geneticmapfile
# python shapeit_to_germline.py Ethiopians_chr21_SHAPEITphased_duoHMM genetic_map_chr21_combined_b37.txt

prefix = argv[-2]
hapsfile = prefix + '.haps'
samplefile = prefix + '.sample'

hapmapfile = argv[-1]

famdata = [line.strip().split() for line in file(samplefile)][2:]

# read in haplotypes
print 'reading in haps data...'
hapsdata = [line.strip().split() for line in file(hapsfile)]
print 'splitting haps data'
mapinfo = array([line[:3] for line in hapsdata])
variants = array([line[3:5] for line in hapsdata])
alleles = array([line[5:] for line in hapsdata],dtype=int)


#hapsdata = array([line.strip().split()[:3] for line in file(hapsfile)])
bps = array( mapinfo[:,2] , dtype=int)

# read in hapmap file
mapdata = array([line.strip().split() for line in file(hapmapfile)][2:],dtype=float)
happos = mapdata[:,0]
hapcms = mapdata[:,2]

# now create interpolation function
f = interp1d(happos,hapcms,fill_value='extrapolate')

bpcms = f(bps)

# make map file
outmap = file(prefix + '.map' , 'w')
print 'writing mapfile to %s' % (outmap)
for i, line in enumerate(mapinfo):
	outmap.write('%s\t%s\t%s\t%s\n' % (line[0],line[1],bpcms[i],line[2]))


# make ped file. More annoying for sure!
# brute force for now, hopefully it runs okay (not huge files)
outped = file(prefix + '.ped' , 'w')
print 'writing pedfile to %s' % (outped)
alleles_acgt = array([variants[i][alleles[i]] for i in range(len(variants))])
for i, line in enumerate(famdata):
	if i % 20 == 0:
		print 'on individual %s' % (i)
	# zip together haplotypes 2*i and 2*i+1 so they're on the same line
	plinkalleles = array(ndarray.flatten(array(zip(alleles_acgt[:,2*i],alleles_acgt[:,2*i+1] ) ) ), dtype='S1')
	outped.write('%s %s\n' % ( ' '.join(line[:2] + line[3:]) , ' '.join(plinkalleles)))

outmap.close()
outped.close()
	
print 'done'
