# Name: sandbox.py
#
# Authors: Fabien Tarissan
# Purpose: graph analysis tools
############################################################################

import sys
import numpy


# Set an instance of graph from infile (list of links)
def loadGraphLinks(infile,verbose=True):
        """read and format graph from file descriptor
            remove self-loop and multi-edges
        Args:
           infile: descriptor of input stream 
        Format: adjacency list
           each line is a link == couple of id nodes
        Returns:
           dict : id -> set of neighbors
        """

        if verbose:
                sys.stderr.write("Load undirected graph\n   Remove self-loops and multi-edges\n") 

        graph={}
        
        for line in infile:
                [u,v]= line.strip().split()
                if u!=v:
                        graph.setdefault(u,set()).add(v)
                        graph.setdefault(v,set()).add(u)

        return graph


# to compute clustering coefficient
def graphCC(graph,verbose=True):
        """Compute the clustering coefficient

	Args:
	    graph : input graph
	Format:
	    dict : id -> set of neighbors
	Returns:
	    dict : id -> clustering coefficient of id
        """

        if verbose:
                sys.stderr.write("Compute clustering coefficient\n") 
                
        cc={}
        for u in graph:
                neigh=graph[u]
                l=len(neigh)
                if l>1:
                        nb_v=l*(l-1)/2
                        nb_tri=0
                        for v in neigh:
                                nb_tri+=len(graph[u].intersection(graph[v]))
                        nb_tri/=2
                        cc[u]=nb_tri/nb_v
                else:
                        cc[u]=-1
        return (cc)


def printGraphStat(graph,outfile, verbose=True):
        """Print global  properties of simple (undirected) graphs

	Args:
	    graph: dict of nodes
	    outfile: descriptor of output stream
	Format:
	    graph is a dict : id -> set of neighbors
        """
        if verbose:
                sys.stderr.write("Print global properties of simple (undirected) graph\n") 

        nb_nodes=len(graph)
        deglist=[len(graph[x]) for x in graph]
        nb_edges=sum(deglist)/2
        nb_edges=0
        for x in graph:
                nb_edges+=len(graph[x])
        nb_edges//=2
        max_deg=max(deglist)
        avdeg=numpy.mean(deglist)
        nb_posslink=nb_nodes*(nb_nodes-1)/2
        density=float(nb_edges)/nb_posslink
        
        cc=graphCC(graph)
        cclist=[cc[x] for x in graph if cc[x]!=-1]
        avcc=numpy.mean(cclist)
        
        
        print("Number of nodes:",nb_nodes, file=outfile)
        print("Number of edges:",nb_edges, file=outfile)
        print("Density:",density, file=outfile)	
        print("Average degree:",avdeg, file=outfile)
        print("Highest degree:",max_deg, file=outfile)
        print("Clustering coefficient:",avcc, file=outfile)


def printDistribution(dist,outfile, verbose=True):
        """Print global  properties of simple (undirected) graphs

	Args:
	    dist: dict of distribution
	    outfile: descriptor of output stream
	Format:
	    graph is a dict : value -> nb of times the value appears
        """
        if verbose:
                sys.stderr.write("Print distribution\n") 


        for v in dist:
                print(v,dist[v], file=outfile)

        
############
# Main programm #
############

# Main loop
if __name__ == '__main__':
        # handling errors
        if len(sys.argv)!=2: # wrong argument number
                sys.stderr.write("Usage: %s infile\n"%sys.argv[0])
                sys.exit(1)


        #load graph 
        infile=open(sys.argv[1],"r")
        graph = loadGraphLinks(infile)
        infile.close()
        
        #compute and print properties
        outfile=open("properties", "w")
        printGraphStat(graph,outfile) 
        outfile.close()

        # computed and print degree distribution
        distdeg={}
        for x in graph:
                d=len(graph[x])
                if d not in distdeg:
                        distdeg[d] = 0
                distdeg[d] = distdeg[d] + 1
                
        outfile=open("degdist.txt", "w")
        printDistribution(distdeg,outfile)
        outfile.close()



        
