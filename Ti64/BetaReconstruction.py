import os
import numpy as np
import matplotlib.pyplot as plt
from mtex import *

# Set the directory path where the data is located
adress = r'C:\Users\mrash\Desktop\EBSD-Ti64\Processing Parameters\ODFs\Max20\\'

# Set the sample ID
sample = '53'

# Load the EBSD data
ebsd = load_ebsd(adress + sample + '.ctf')

# Calculate the Burgers orientation between Beta and Alpha phase
beta2alpha = orientation.Burgers(ebsd('Beta').CS,ebsd('Alpha').CS)

# Calculate grains from the indexed EBSD data
grains,ebsd['indexed'].grainId = calcGrains(ebsd['indexed'],'threshold',1.5*degree)

# Reconstruct the parent grains from the indexed EBSD data
job = parentGrainReconstructor(ebsd, grains)
job.p2c = beta2alpha
job.calcTPVotes('minFit',2.5*degree,'maxFit',5*degree)
job.calcParentFromVote()

# Calculate grain boundary votes and reconstruct parent grains
for k in range(1,4):
    job.calcGBVotes('p2c','threshold',k * 2.5*degree)
    job.calcParentFromVote()

# Merge similar grains and inclusions
job.mergeSimilar('threshold',5*degree)
job.mergeInclusions()

# Get the parent EBSD data
parentEBSD = job.ebsd

# Set the IPF color key and direction
ipfKey = ipfColorKey(ebsd('Beta'))
ipfKey.inversePoleFigureDirection = vector3d.Z

# Plot the Alpha phase with its orientations
color = ipfKey.orientation2color(parentEBSD('Beta').orientations)
plot(ebsd('Alpha'),ebsd('Alpha').orientations,'figSize','large')

# Overlay the parent grains' boundaries
hold('on')
parentGrains = smooth(job.grains,5)
plot(parentGrains.boundary,'lineWidth',3)
hold('off')

# Calculate and plot the Beta phase ODF
odf = calcDensity(parentEBSD('Beta').orientations)
plotSection(odf)
mtexColorbar()

# Calculate and plot the Beta phase texture index
psi = calcKernel(parentEBSD('Beta').orientations)
psi = calcKernel(parentGrains('Beta').meanOrientation)
odf = calcDensity(parentEBSD('Beta').orientations,'kernel',psi)
h = [Miller(0,0,0,2,odf.CS)]
plotPDF(odf,h,'antipodal','silent')

# Calculate and plot the Alpha phase ODF
odf = calcDensity(ebsd('Alpha').orientations)
plotSection(odf, 'phi2', 0*degree)
caxis([0,30])
mtexColorbar()

# Save the Alpha phase ODF plot as a TIFF image file
fullFileName = os.path.join(adress, sample + '-AlphaODF-Zero.tiff')
imagewd = getframe(gcf())
imwrite(imagewd.cdata, fullFileName, 'Compression', 'none', 'Resolution', 300)

# Close all plots
plt.close('all')

# Calculate and plot the texture index of the Alpha phase
textureindex(odf)
psi = calcKernel(grains('Alpha').meanOrientation)
odf = calcDensity(ebsd('Alpha').orientations,'kernel',psi)
h = [Miller(0,1,2,odf.CS)]
plotPDF(odf,h)