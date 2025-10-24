# Name: plotfig.py
#
# Authors: F. Tarissan
# Purpose: plot a distribution
#####################################

import sys
from pylab import *

def load_file(fdin):
        lx=[]
        ly=[]
        for line in fdin:
                l=line.strip().split()
                lx.append(float(l[0]))
                ly.append(float(l[1]))

        return (lx, ly)

# Main loop
if __name__ == '__main__':

        outfig="distribution.pdf"
        
        # properties of figure
        fig=figure(figsize=(12,9))
        subplots_adjust(left = 0.1, bottom = 0.11, right = 0.98, top = 0.9, wspace = 0.5, hspace = 0.17)

        sp=fig.add_subplot(1,1,1)

        f=open(sys.argv[1],'r')
        (lx,ly) = load_file(f)
        f.close()

        sp.plot(lx,ly,'+', markersize=12)

        sp.set_title("Distribution", fontsize=26)
        sp.set_xlabel("Degree",fontsize=17)
        sp.set_ylabel("Nb nomdes with degree = x",fontsize=17)

        savefig(outfig, format="pdf")
