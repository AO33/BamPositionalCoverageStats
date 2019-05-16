### Imports ###
import time
import argparse
import pysam
import multiprocessing as mp
############################

### Code ###
class Position:
	''' This is a class that is used to calculate alignment statistics over a specified position in a bam file.
		It is suppose to be initiated by reading a bedFile where the reference and alternate alleles are alread specified.
		There are two methods assosicated with this class that will get the positional info over the specified position,
		and another that will format a string for output to a tsv file.'''
	def __init__(self,chromosome,position,ref,alt):
		self.chromosome = chromosome
		self.position = position
		self.ref = ref
		self.alt = alt
	##################################
	def getInfo(self,samFile,mapQ=60):
		self.mappingQuality = mapQ
		self.nucDict = {"A":0,"C":0,"T":0,"G":0,"N":0}
		self.coverage,self.indelReads = 0.0,0.0
		self.properPairs,self.mapQ = 0.0,0.0
		pysamAdjustedPos = self.position - 1
		for pileupcolumn in samFile.pileup(reference=self.chromosome,max_depth=100000000,start=pysamAdjustedPos,end=self.position):
			if pileupcolumn.pos == (pysamAdjustedPos):
				for pileUpRead in pileupcolumn.pileups:
					### Check if the read has a "deletion over this position" ###
					if pileUpRead.is_del == 1:
						self.indelReads += 1.0
					else:
						### Deal with reads that actually have coverage over this position and not a deletion ###
						self.nucDict[pileUpRead.alignment.query_sequence[pileUpRead.query_position]] += 1.0
						self.coverage += 1
					if pileUpRead.alignment.is_proper_pair:
						self.properPairs += 1
					### For bwa mem, this means a read is of the highest mapping quality and uniquely mapped 
					if pileUpRead.alignment.mapping_quality >= mapQ:
						self.mapQ += 1.0
	########################
	def getAnnotation(self):
		self.ann = self.chromosome+'\t'+str(self.position)+'\t'+self.ref+'\t'+"Coverage:"+str(self.coverage)+'\t'
		self.ann = self.ann+'A'+":"+str(self.nucDict["A"])+'\t'+'C'+":"+str(self.nucDict["C"])+'\t'+'T'+":"+str(self.nucDict["T"])+'\t'+'G'+":"+str(self.nucDict["G"])+'\t'
		self.ann = self.ann+'N'+":"+str(self.nucDict["N"])+'\t'"INDELReads:"+str(self.indelReads)+'\t'+"MAPQ>="+str(self.mappingQuality)+":"+str(self.mapQ)+'\t'+"PROPER_PAIRED:"+str(self.properPairs)

def reportPositionStats(ls,samMan,mapQ=60):
	'''This function takes in a list of bed file entries and will calculate positional statistics over each entry.
		It inititally opens up a bamFile via pysam and then loops through the bedFile entries calling the Position class.object and getInfo method from the object.
		It then returns a list of Position objects where each one has had the getInfo method called on it'''
	samFile = pysam.AlignmentFile(samMan,"rb")
	returnArray = []
	for line in ls:
		cols = line.strip('\r').strip('\n').split('\t')
		chrom,start,ref,alt = cols[0],int(cols[1]),cols[3],cols[4]
		posMan = Position(chrom,start,ref,alt)
		posMan.getInfo(samFile,mapQ=mapQ)
		returnArray.append(posMan)
	###############
	samFile.close()
	return returnArray

def readBedFileToMem(bedFile):
	'''Pretty straight foreward. Reads through a bed file and appends each line to a list'''
	bedFile = open(bedFile,'r')
	bedLines = [ line for line in bedFile ]
	bedFile.close()
	return bedLines

def divideWorkLoad(items,proc=1):
	'''Clever little script that I paritally stole from the internet, that divides up a list into equal parts.
	   Or as equal as possible which is ideal for "parallel processing" because each node will get an equal work load.
	   Other methods Ive found tend to have 1 node that gets an unequal distribution of work so it just sits there.'''
	if proc == 1:
		return items
	######################################
	baskets = [ [] for i in xrange(proc) ]
	for i,item in enumerate(items):
		baskets[ i % proc ].append(item)
	##############################
	baskets = filter(None,baskets)
	return baskets

def grabCoverageSNPs(bedFile,samFileLocation,proc=1):
	''' Master function that will report position statistics for each entry in a bed file.
		Uses multiprocessing library to distribute the work load evenly among cores.
		Runtime is essentially linear with respect to the number of processors which is ideal.
		Will write out the results to a file.
		bedFile = fullFilePathToBedFile. In format of chromosome'\t'position'\t'position'\t'referenceAllele'\t'alternateAllele
		samFileLocation = fullFilePathToBamFile
		I refer to the input bam file as samFileLocation because pysam opens the bam file as a sam essentially.'''
	startTime = time.time()
	if proc != 1:
		bedLineLists = divideWorkLoad(readBedFileToMem(bedFile),proc=proc)
		output = mp.Queue()
		pool = mp.Pool(processes=proc)
		print "Workload divided"
		print "Calculating position stats across "+str(proc)+" cores..."
		results = [ pool.apply_async(reportPositionStats,args=(bedLines,samFileLocation)) for bedLines in bedLineLists ]
		output = [ p.get() for p in results ]
	else:
		output = [ [pObj] for pObj in reportPositionStats(readBedFileToMem(bedFile),samFileLocation) ]
	###################################
	print "Positional stats calculated"
	print "Ordering output..."
	### Now we want to order our positionObjects by chromosome and position (smalles position first) ###
    ### This just makes the output easier to understand/follow instead of it just being essentially random ###
	chrDict = {}
	for posObjs in output:
		for posObj in posObjs:
			if posObj.chromosome not in chrDict:
				chrDict.update({ posObj.chromosome:[posObj] })
			else:
				chrDict[posObj.chromosome].append(posObj)
    #################
	for c in chrDict:
		chrDict[c] = sorted(chrDict[c],key=lambda posObj: posObj.position)
	######################
	stopTime = time.time()
	print "RunTime:"+'\t'+str(stopTime-startTime)
	return chrDict

def writePositionsToFile(outFile,chrDict,chromosomeOrderList="NA"):
	'''Also pretty self explanatory. Writes the results of the positional stats objects to a file.
		outFile = fullFilePathToOutputFile. Should include the name you want to name the outFile
		chromosomeOrderList = [ chrName, chrName, ... ] if you want to order the output file in a chromosome speficic manner you can'''
	if chromosomeOrderList == "NA":
		chromosomeOrderList = [ c for c in chrDict ]
	###########################
	outFile = open(outFile,'w')
	for c in chromosomeOrderList:
		for posObj in chrDict[c]:
			posObj.getAnnotation()
			outFile.write(posObj.ann+'\n')
	###############
	outFile.close()
###################


if __name__ == "__main__":
	###############
	## ARG PARSE ##
	###############
	def mm():
		parser = argparse.ArgumentParser(description='Calculates coverage stats from an alignment file (bam format) for specified positions')
		parser.add_argument("-bedFile",help="Standard bed file format with the fourth and fiftch columns being referenceAllele and alternateAllele respectively",required=True,type=str)
		parser.add_argument("-outFile",help="OutPut file to write results to", required=True, type=str)
		parser.add_argument("-bamFile",help="Location of alignment file in bam format. Note, it needs to be indexed first with something like samtools index.",required=True,type=str)
		parser.add_argument("-proc",help="Number of cores to use. Runs ~linearly with the number of cores give. Drastically speeds up compute time for large data sets",required=True,type=int)
		parser.add_argument("-chrList",help="Comma separated list of chromosome names. NO spaces contained with this argument. Example: chr1,chr2,chr3",required=False,type=str)
		return parser.parse_args()
	###########
	args = mm()
	if args.chrList:
		chromList = args.chrList.split(',')
	else:
		chromList = "NA"
	###########################################################################################################################
	writePositionsToFile(args.outFile,grabCoverageSNPs(args.bedFile,args.bamFile,proc=args.proc),chromosomeOrderList=chromList)


### Code is pretty flexible to add in more arguments like minimum mapping quality etc...
### This is just the basic version now tailored for a particular project ive been working on
