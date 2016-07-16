import time
import argparse
import numpy as np
from blist import blist
from scipy.stats import norm

# User Variables
parser = argparse.ArgumentParser(description="Detect blocks of overlapping reads using a gaussian distribution approach.")
parser.add_argument("file", help="input file")
parser.add_argument("--minblockheight", type=int, default=2, help="minimum block height (default: 2)")
parser.add_argument("--sizescale", type=float, default=0.2, help="scale stddev for a single read (default: 0.2)")
parser.add_argument("--merge", type=int, default=0, help="merge reads with almost similar means (default: 0)")
parser.add_argument("--distance", type=int, default=30, help="minimum distance between two clusters (default: 30)")
parser.add_argument("--minClusterHeight", type=int, default=10, help="minimum cluster height (default: 10)")
parser.add_argument("--printout", type=int, default=1, help="print out: (1) blocks (2) reads (default: 1)")
parser.add_argument("--tagFilter", type=float, default=0.0, help="skip tags with expression smaller than this value (default: 0.0)")
parser.add_argument("--type", type=int, default=1, help="file format of input file (default: 1):\n(1)bed\n(2)segemehl-output")
parser.add_argument("--sep", default="\t", help="separator (default: '\t')")

try:
    args = parser.parse_args()
except:
    print('\nThere was an error while parsing the arguments.\n')
    parser.print_help()

# Global Variables
clusterStart = -1
clusterEnd = -1
clusterHeight = 0
readCount = 0
tagCount = 0
clusterChrom = "x"
clusterStrand = "x"
clusterCounter = 0


class read:
    def __init__(self, chrom=None, start=None, end=None, block=-1, height=None, id=None, strand=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.block = block
        self.height = height
        self.id = id
        self.strand = strand


def writeHeader():
    print("# blockbuster result file generated %(t)s\n# query file: %(filename)s\n# scale: %(sizescale).1f, minblockheight: %(minblockheight)i, mergeDistance: %(merge)i\n# block_number\tchromosome\tstart_of_block\tend_of_block\tstrand\treadIDs\n"
          % {'t': time.ctime(None), 'filename': args.file, 'sizescale': args.sizescale, 'minblockheight': args.minblockheight, 'merge': args.merge})


# CALCULATE THE STANDARD DEVIATION
def stddev(readMeans, readHeights):
    s = 0
    counter = 0
    for i in range(tagCount):
        if readMeans[i] != -1:
            s += (np.int(readHeights[i]) * readMeans[i])
            counter += np.int(readHeights[i])

    if counter == 0:
        return 0
    mean = s / counter

    s = 0
    for i in range(tagCount):
        if (readMeans[i] != -1):
            s += (np.int(readHeights[i]) * np.power((readMeans[i] - mean), 2))

    return np.sqrt(s / counter)


# CALCULATE THE GAUSSIAN DISTRIBUTIONS OF ALL READS AND SUM THEM UP
def writeSuperGaussian(anchor, distrib, clusterSize, clusterStart):
    for r in anchor:
        if r.block == -1:
            mean = ((r.start + r.end) / 2) - clusterStart
            variance = args.sizescale * (abs(r.end - r.start) / 2)

            x = [i + mean for i in range(int(2 * variance) + 1)]
            y = r.height * norm.pdf(x, mean, variance)
            for i in range(int(2 * variance) + 1):
                if (int(mean + i) < clusterSize):
                    distrib[int(mean + i)] += y[i]
                if (int(mean - i) > 0):
                    distrib[int(mean - i)] += y[i]


# ASSIGN READS TO A BLOCK
def assignReads(anchor, highestPeak, clusterSize, blockCount):
    readMeans = blist([-1])
    readHeights = blist([-1])
    readMeans *= tagCount
    readHeights *= tagCount
    meanCounter = 0
    blockHeight = 0
    block = list()

    counterNew = 0
    counterOld = -1

    while counterOld != counterNew:
        dev = stddev(readMeans, readHeights)
        counterOld = counterNew
        for read in anchor[:]:
            mean = ((read.start + read.end) / 2) - clusterStart
            variance = args.sizescale * (abs(read.end - read.start) / 2)

            if (((mean - variance - dev) <= highestPeak and (mean + variance + dev) >= highestPeak) or (mean >= (highestPeak - args.merge) and mean <= (highestPeak + args.merge))):
                readMeans[meanCounter] = mean
                readHeights[meanCounter] = read.height
                meanCounter += 1
                blockHeight += read.height
                read.block = blockCount
                block.append(read)
                anchor.remove(read)
                counterNew += 1

    if blockHeight >= args.minblockheight:
        return block
    else:
        return []


def assignReadsToBlocks(anchor):
    blockCount = 1
    cluster = dict()
    global readCount
    global clusterEnd
    global clusterStart
    global tagCount

    # create an array with clusterSize entries for the superGaussian distribution
    clusterSize = (clusterEnd - clusterStart)
    distrib = np.zeros(clusterSize, dtype=np.double)

    old = 1
    new = 0

    # run through sorted peaks
    while(old != new):
        old = len(anchor)

        # clean distribution array
        distrib = np.zeros(clusterSize, dtype=np.dtype('d'))

        # write distribution
        writeSuperGaussian(anchor, distrib, clusterSize, clusterStart)
        highestPeakIndex = np.argmax(distrib)
        distrib[highestPeakIndex] = 0

        # assign reads to the highest peak
        block = assignReads(anchor, highestPeakIndex, readCount, blockCount)
        if len(block) != 0:
            cluster[blockCount] = block
            blockCount += 1
        new = len(anchor)

    writeBlocks(cluster)


# WRITE THE READS THAT ARE ASSIGNED TO A BLOCK TO STDOUT
def writeBlocks(cluster):
    global clusterCounter
    global clusterChrom
    global clusterStart
    global clusterEnd
    global clusterStrand
    absTagCount = 0
    absClusterHeight = 0

    # get cluster information
    for block in cluster:
        for read in cluster[block]:
            absClusterHeight += read.height
            absTagCount += 1

    if len(cluster) > 0:
        clusterCounter += 1
        # print header
        print(">cluster_%(cC)i\t%(cCh)s\t%(cS)i\t%(cE)i\t%(cStr)s\t%(absCH).2f\t%(absTC)i\t%(blockNb)i"
              % {'cC': clusterCounter, 'cCh': clusterChrom, 'cS': clusterStart, 'cE': clusterEnd,
                 'cStr': clusterStrand, 'absCH': absClusterHeight, 'absTC': absTagCount, 'blockNb': len(cluster)})

        # print blocks
        if args.printout == 1:
            for block in cluster:
                thisBlockHeight = 0
                thisBlockTags = 0
                thisBlockStart = -1
                thisBlockEnd = -1
                for read in cluster[block]:
                    if thisBlockStart == -1:
                        thisBlockStart = read.start
                        thisBlockEnd = read.end
                    if read.start < thisBlockStart:
                        thisBlockStart = read.start
                    if read.end > thisBlockEnd:
                        thisBlockEnd = read.end
                    thisBlockHeight += read.height
                    thisBlockTags += 1
                print("%(wB)i\t%(cCh)s\t%(tBS)i\t%(tBE)i\t%(cS)s\t%(tBH).2f\t%(tBT)i"
                      % {'wB': block, 'cCh': clusterChrom, 'tBS': thisBlockStart, 'tBE': thisBlockEnd,
                         'cS': clusterStrand, 'tBH': thisBlockHeight, 'tBT': thisBlockTags})

        # print tags
        if args.printout == 2:
            for block in cluster:
                thisBlockHeight = 0
                for read in cluster[block]:
                    print("%(sCh)s\t%(ss)d\t%(se)d\t%(si)s\t%(sh)lf\t%(sst)s\t%(wB)i"
                          % {'sCh': read.chrom, 'ss': read.start, 'se': read.end, 'si': read.id,
                             'sh': read.height, 'sst': read.strand, 'wB': block})


def read_bed_file(filename):
    global clusterHeight
    global clusterStart
    global clusterEnd
    global tagCount
    global clusterStrand
    global clusterChrom
    thisCluster = blist([])
    lastEnd = -1
    lastChrom = "x"
    lastStrand = "x"
    try:
        f = open(filename, "r")
    except:
        print("cannot open %(filename)s\n" % {'filename': filename})
    else:
        header = 0
        for line in iter(f.readline, ''):
            # if "#" at the beginning of line -> header
            if line[0] == '#':
                header = 1
                writeHeader()

            # parse line
            if header != 1:
                chrom, id, strand, info = None, None, None, None
                start, end, height, freq = None, None, None, -1

                if (args.type) == 1:
                    # run through line and split at separator
                    linelist = blist(line.split())
                    try:
                        chrom = linelist[0]
                        start = int(linelist[1])
                        end = int(linelist[2])
                        id = linelist[3]
                        height = float(linelist[4])
                        strand = linelist[5]
                    except:
                        print('\nwrong file format\n')
                        parser.print_help()
                        exit(0)

                elif (args.type) == 2:
                    # run through line and split at separator
                    linelist = blist(line.split())
                    try:
                        info = linelist[1]
                        strand = int(linelist[9])
                        assert strand is not None
                        start = int(linelist[10])
                        end = int(linelist[11])
                        chrom = linelist[12]
                        freq = float(linelist[14])
                    except:
                        print('\nwrong file format\n')
                        parser.print_help()
                        exit(0)

                    # split id|freq
                    q = blist(info.split('|'))
                    for i in range(len(q)):
                        if i == 0:
                            id = q[i]
                        if i == 1:
                            height = float(q[i])

                    if height != -1:
                        height /= freq
                    else:
                        height = 1 / freq

                if (height >= (args.tagFilter)):
                    if (chrom != lastChrom) or (strand != lastStrand) or ((start - lastEnd) > (args.distance)):
                        if ((clusterHeight) > (args.minClusterHeight)):
                            # Analyze Cluster
                            assignReadsToBlocks(thisCluster)

                        thisCluster.clear()

                        # reset cluster dimensions
                        clusterStart = start
                        clusterEnd = end
                        clusterHeight = height
                        tagCount = 1

                    else:
                        # update cluster dimensions
                        if clusterStart > start:
                            clusterStart = start
                        if clusterEnd < end:
                            clusterEnd = end
                        clusterHeight += height
                        tagCount += 1

                    thisRead = read(chrom=chrom, start=start, end=end, height=height, id=id, strand=strand)
                    thisCluster.append(thisRead)

                    clusterChrom = chrom
                    clusterStrand = strand
                    lastChrom = chrom
                    lastStrand = strand
                    lastEnd = end

            elif header == 1:
                header = 0

        assignReadsToBlocks(thisCluster)
        f.close()

try:
    read_bed_file(args.file)
except:
    print('\nError while reading the file.')
    exit(0)
