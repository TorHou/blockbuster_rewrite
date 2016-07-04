import time
import argparse
import numpy as np
# from scipy import stats
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

class cluster:
    def __init__(self):
        self.head = None

    def add_read(self, chrom=None, start=None, end=None, block=-1, height=None, id=None, strand=None):
        new_read = read(chrom=None, start=None, end=None, block=-1, height=None, id=None, strand=None)
        new_read.next = self.head
        self.head = new_read

    def free(self):
        self.head = None

    def assignReadsToBlocks(self):
        blockCount = 1
        global readCount
        global clusterEnd
        global clusterStart

        # create an array with clusterSize entries for the superGaussian distribution
        clusterSize = (clusterEnd - clusterStart)
        distrib = np.zeros(clusterSize, dtype=np.double)

        old = 1
        new = 0

        # run through sorted peaks
        while(old != new):
            old = self.getRest()

            # clean distribution array
            distrib = np.zeros(clusterSize, dtype=np.dtype('d'))

            # write distribution
            self.writeSuperGaussian(distrib, clusterSize)
            highestPeakIndex = np.argmax(distrib)
            distrib[highestPeakIndex] = 0

            # assign reads to the highest peak
            sum = self.assignReads(highestPeakIndex, readCount, blockCount)
            if sum != 0:
                blockCount += 1
            new = self.getRest()

    # GET THE AMOUNT OF READS NOT ASSIGNED TO BLOCKS
    def getRest(self):
        sum = 0
        while self.head != None:
            if self.head.block == -1:
                sum += 1
            self.head = self.head.next
        return sum

    # CALCULATE THE GAUSSIAN DISTRIBUTIONS OF ALL READS AND SUM THEM UP
    def writeSuperGaussian(self, distrib, clusterSize):
        global clusterStart

        while self.head != None:
            if self.head.block == -1:
                mean = ((r.start + r.end) / 2) - clusterStart
                variance = args.sizescale * (abs(r.end - r.start) / 2)

                x = [i + mean for i in range(int(2 * variance) + 1)]
                y = r.height * norm.pdf(x, mean, variance)
                for i in range(int(2 * variance) + 1):
                    x = mean + i
                    if (int(x) < clusterSize):
                        distrib[int(x)] += y[i]
                    if (int(mean - i) > 0):
                        distrib[int(mean - i)] += y[i]
            self.head = self.head.next

    # ASSIGN READS TO A BLOCK
    def assignReads(self, highestPeak, clusterSize, blockCount):
        global tagCount
        global clusterStart
        readMeans = -1 * np.ones(tagCount, dtype=np.double)
        readHeights = -1 * np.ones(tagCount, dtype=np.double)
        meanCounter = 0

        counterNew = 0
        counterOld = -1

        while counterOld != counterNew:
            dev = stddev(readMeans, readHeights)
            counterOld = counterNew
            while self.head != None:
                if self.head.block == -1:
                    mean = ((self.head.start + self.head.end) / 2) - clusterStart
                    variance = args.sizescale * (abs(self.head.end - self.head.start) / 2)

                    if (((mean - variance - dev) <= highestPeak and (mean + variance + dev) >= highestPeak) or (mean >= (highestPeak - args.merge) and mean <= (highestPeak + args.merge))):
                        readMeans[meanCounter] = mean
                        readHeights[meanCounter] = self.head.height
                        meanCounter += 1
                        self.head.block = blockCount
                        counterNew += 1
                self.head = self.head.next

        return counterNew

    # WRITE THE READS THAT ARE ASSIGNED TO A BLOCK TO STDOUT
    def writeBlocks(self):
        global clusterCounter
        global clusterChrom
        global clusterStart
        global clusterEnd
        global clusterStrand
        thisBlock = 0
        blockHeight = 0
        blockNb = 0
        thisTagCount = 0
        thisClusterHeight = 0
        absTagCount = 0
        absClusterHeight = 0
        size = 1

        # get cluster information
        while size > 0:
            thisBlock += 1

            # reset variables
            size = 0
            blockHeight = 0
            thisClusterHeight = 0
            thisTagCount = 0

            # run through linked list of reads
            while self.head != None:
                # if current read is in thisBlock
                if self.head.block == thisBlock:
                    size += 1
                    blockHeight += self.head.height
                    thisClusterHeight += self.head.height
                    thisTagCount += 1
                self.head = self.head.next

            # check if block is high enough
            if (blockHeight >= args.minblockheight) and (size > 0):
                blockNb += 1
                absClusterHeight += thisClusterHeight
                absTagCount += thisTagCount

        if blockNb > 0:
            clusterCounter += 1
            # print header
            print(">cluster_%(cC)i\t%(cCh)s\t%(cS)i\t%(cE)i\t%(cStr)s\t%(absCH).2f\t%(absTC)i\t%(blockNb)i"
                  % {'cC': clusterCounter, 'cCh': clusterChrom, 'cS': clusterStart, 'cE': clusterEnd,
                     'cStr': clusterStrand, 'absCH': absClusterHeight, 'absTC': absTagCount, 'blockNb': blockNb})

            # print blocks
            if args.printout == 1:
                thisBlock = 0
                size = 1
                writeBlock = 0
                while size > 0:
                    thisBlock += 1
                    size = 0
                    thisBlockHeight = 0
                    thisBlockTags = 0
                    thisBlockStart = -1
                    thisBlockEnd = -1
                    while self.head != None:
                        if self.head.block == thisBlock:
                            if thisBlockStart == -1:
                                thisBlockStart = self.head.start
                                thisBlockEnd = self.head.end
                            if self.head.start < thisBlockStart:
                                thisBlockStart = self.head.start
                            if self.head.end > thisBlockEnd:
                                thisBlockEnd = self.head.end
                            thisBlockHeight += self.head.height
                            thisBlockTags += 1
                            size += 1
                        self.head = self.head.next
                    if (thisBlockHeight >= args.minblockheight) and (size > 0):
                        writeBlock += 1
                        print("%(wB)i\t%(cCh)s\t%(tBS)i\t%(tBE)i\t%(cS)s\t%(tBH).2f\t%(tBT)i"
                              % {'wB': writeBlock, 'cCh': clusterChrom, 'tBS': thisBlockStart, 'tBE': thisBlockEnd,
                                 'cS': clusterStrand, 'tBH': thisBlockHeight, 'tBT': thisBlockTags})

            # print tags
            if args.printout == 2:
                thisBlock = 0
                size = 1
                writeBlock = 0
                while size > 0:
                    thisBlockHeight = 0
                    thisBlock += 1
                    size = 0
                    while self.head != None:
                        if self.head.block == thisBlock:
                            thisBlockHeight += self.head.height
                            size += 1
                        self.head = self.head.next
                    if (thisBlockHeight >= args.minblockheight) and (size > 0):
                        writeBlock += 1
                        while self.head != None:
                            if self.head.block == thisBlock:
                                print("%(sCh)s\t%(ss)d\t%(se)d\t%(si)s\t%(sh)lf\t%(sst)s\t%(wB)i"
                                      % {'sCh': self.head.chrom, 'ss': self.head.start, 'se': self.head.end, 'si': self.head.id,
                                         'sh': self.head.height, 'sst': self.head.strand, 'wB': writeBlock})
                            self.head = self.head.next


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


def read_bed_file(filename):
    global clusterHeight
    global clusterStart
    global clusterEnd
    global tagCount
    global clusterStrand
    global clusterChrom
    thisCluster = cluster()
    lastEnd = -1
    lastChrom = "x"
    lastStrand = "x"
    try:
        f = open(filename, "r")
    except:
        print("cannot open %(filename)s\n" % {'filename': filename})
    else:
        header = 0
        while True:
            line = f.readline()
            if not line:
                break
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
                    linelist = line.split()
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
                    linelist = line.split()
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
                    q = info.split('|')
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
                            thisCluster.assignReadsToBlocks()
                            thisCluster.writeBlocks()

                        thisCluster.free()

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

                    thisCluster.add_read(chrom=chrom, start=start, end=end, height=height, id=id, strand=strand)

                    clusterChrom = chrom
                    clusterStrand = strand
                    lastChrom = chrom
                    lastStrand = strand
                    lastEnd = end

            elif header == 1:
                header = 0

        thisCluster.assignReadsToBlocks()
        thisCluster.writeBlocks()
        f.close()

try:
    read_bed_file(args.file)
except:
    print('\nError while reading the file.')
    exit(0)
