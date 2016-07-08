import sys
import argparse

parser = argparse.ArgumentParser(description="Detect blocks of overlapping reads using a gaussian distribution approach.")
parser.add_argument("file1", help="input file 1")
parser.add_argument("file2", help="input file 2")


args = parser.parse_args()

def get_cluster(f):
    cluster = {}
    cluster["header"] = ""
    cluster["blocks"] = []
    for line in f:
        if line.startswith(">cluster"):
            if cluster["header"]:
                yield cluster
                cluster["blocks"] = []
                cluster["header"] = line
            else:
                cluster["header"] = line
        else:
            cluster["blocks"].append(line.split()[2:])
                
f1 = open(args.file1) 
f2 = open(args.file2) 
clusters1 = get_cluster(f1)
clusters2 = get_cluster(f2)
while True:
    cluster1 = clusters1.next()
    cluster2 = clusters2.next()
    if not cluster1["header"] == cluster2["header"]:
        print "Different clusters found."
        print cluster1
        sys.exit(1)
            
    for block1 in cluster1["blocks"]:
        found = False
        for block2 in cluster2["blocks"]:
            if block2 == block1:
                found = True
        if not found:
            print "block not found."
            print block1
            sys.exit(1)
    #print "All good with clusters: " + cluster1["header"]

sys.exit(0)
        


        
