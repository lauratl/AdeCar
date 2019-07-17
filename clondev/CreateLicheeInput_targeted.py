from sys import argv


def parseHeader(linea):
    
    samples=[]
    
    linea = linea.split(" ")
    
    for field in linea:
        if field[:9]=="REF_COUNT":
            samples.append("_".join(field.split("_")[2:5]))
            
    return samples
    

def getMutsInfo(linea):
    global muts
        
    linea = linea.split(" ")
        
    mutId = len(muts)
    muts[mutId]={}
    muts[mutId]["chr"]= linea[0].split(":")[0]
    muts[mutId]["pos"]= linea[0].split(":")[1]
    muts[mutId]["samplesInfo"]= {}
    for sampleIdx in range(len(samples)):
        sample = samples[sampleIdx]
        muts[mutId]["samplesInfo"][sample] = {}
        muts[mutId]["samplesInfo"][sample]["ref_count"]=int(linea[4*sampleIdx+1])
        muts[mutId]["samplesInfo"][sample]["alt_count"]=int(linea[4*sampleIdx+2])
        muts[mutId]["samplesInfo"][sample]["ref_nuc"]=linea[4*sampleIdx+3]
        muts[mutId]["samplesInfo"][sample]["alt_nuc"]=linea[4*sampleIdx+4]
            
def keepMultiallelic():
    majorAllele = getMajotAllele()


def getTotalDepth(mut):
    totalDepth = 0
    for sampleInfo in muts[mut]['samplesInfo'].values():
        totalDepth += int(sampleInfo['alt_count']) + int(sampleInfo['ref_count'])
    return totalDepth
        
def getPerSampleDepth(mut):
    sampleDepths = {}
    for sample in muts[mut]['samplesInfo'].keys():
        sampleInfo = muts[mut]['samplesInfo'][sample]
        sampleDepths[sample] = sampleInfo['alt_count'] + sampleInfo['ref_count']
    return sampleDepths

def getPerSampleVaf(mut):
    sampleVafs = {}
    for sample in muts[mut]['samplesInfo'].keys():
        sampleInfo = muts[mut]['samplesInfo'][sample]
        sampleVafs[sample] = float(sampleInfo['alt_count']) /  (sampleInfo['ref_count'] + sampleInfo['alt_count'])
    return sampleVafs

    
def filterLowDepth(mut):

    sampleDepths = getPerSampleDepth(mut)
    for depth in sampleDepths.values():
        if depth < minDepth:
            return True
    return False


def removeMut(mut):
    global muts
    del muts[mut]
    
def getAltAlleles(mut):
    altAlleles = []
    for sampleInfo in muts[mut]['samplesInfo'].values():
        altAlleles.append(sampleInfo['alt_nuc'])
        
    altAlleles = [x for x in altAlleles if x != "N"]
    return altAlleles

def getAltAllelesCount(mut):
    altAlleles = getAltAlleles(mut)
    altAllelesCount = {}
    for sampleInfo in muts[mut]['samplesInfo'].values():
        altAllele = sampleInfo['alt_nuc']
        if altAllele in altAlleles:
            if altAllele in altAllelesCount.keys():
                altAllelesCount[altAllele]+= sampleInfo['alt_count']
            else:
                altAllelesCount[altAllele]= sampleInfo['alt_count']
    return altAllelesCount

def checkMultiallelic(mut):
    
    altAlleles = getAltAlleles(mut)

    if len(set(altAlleles)) > 1:
        return True
    #elif len(set(altAlleles)) > 2:
    #    print "ERROR: too many alternative alleles"
    #    print mut
    #    print muts[mut]['pos']
    return False

def getMajorAllele(altAllelesCount):
    
    majorAlleleCount = 0
    majorAllele=""
    
    for alt in altAllelesCount.keys():
        if altAllelesCount[alt]>majorAlleleCount:
            majorAllele =alt
            majorAlleleCount= altAllelesCount[alt]
    return majorAllele

def getMinorAllele(altAllelesCount):
    
    minorAlleleCount = float('Inf')
    minorAllele=""
    
    for alt in altAllelesCount.keys():
        if altAllelesCount[alt]< minorAlleleCount:
            minorAllele =alt
            minorAlleleCount= altAllelesCount[alt]
    return minorAllele

def getMinorAlleles(majorAllele):
    altAlleles = getAltAlleles(mut)
    minorAlleles = [x for x in altAlleles if x != majorAllele]
    return minorAlleles
    
    
    
def keepMultiallelic(mut):
    
    altAlleles = getAltAlleles(mut)
    altAllelesCount = getAltAllelesCount(mut)
    majorAllele = getMajorAllele(altAllelesCount)
    
    if len(set(altAlleles)) == 2:
        minorAllele = getMinorAllele(altAllelesCount)
        alleleRatio = float(altAllelesCount[minorAllele]) / altAllelesCount[majorAllele]
        
    elif len(set(altAlleles)) > 2:
        minorAlleles = getMinorAlleles(majorAllele)
        
        minorAllelesCount = []
        for minorAllele in minorAlleles:
            minorAllelesCount.append(altAllelesCount[minorAllele])

        alleleRatio = float(max(minorAllelesCount)) / altAllelesCount[majorAllele]
    
    
    if alleleRatio < 0.2:
        return True
    else:
        return False
            
    
def removeMinorAllele(mut):
    
    global muts
    
    altAllelesCount = getAltAllelesCount(mut)
    majorAllele = getMajorAllele(altAllelesCount)
    
    for sample in muts[mut]['samplesInfo'].keys():
        if muts[mut]['samplesInfo'][sample]['alt_nuc'] != majorAllele:
            muts[mut]['samplesInfo'][sample]['alt_count'] = 0
            muts[mut]['samplesInfo'][sample]['alt_nuc'] = "N"
    
def filterLowMaxVaf(mut):
    sampleVafs = getPerSampleVaf(mut)

    for vaf in sampleVafs.values():
        if vaf > minVaf:
            return False
    return True


def filterGermline(mut):
    sampleVafs = getPerSampleVaf(mut)
    if sampleVafs[healthy] > germlineVaf:
        return True
    return False


def withHealthy():
    if healthy=="":
        return False
    return True

def filterHighFreq(mut):
    
    vafs = getPerSampleVaf(mut)
    for vaf in vafs.values():
        if vaf > maxVafIfNotHealthy:
            return True
    return False


def isPanelError(mut):
    
    panelErrors = { 'chr19' : ["15285135"],
                    'chr21' : ["326385", 
                               "32638549", 
                               "32638550"],
                    'chr8'  : ["13356818"] }
    
    chrom = muts[mut]['chr']
    pos = muts[mut]['pos']
    
    if chrom in panelErrors.keys() and pos in panelErrors[chrom]:
        return True
    return False
        
    
    


def printLicheeHeader():
    
    licheeOutputFile.write("#chr\tposition\tdescription")
    if withHealthy():
        licheeOutputFile.write("\t"+healthy)
    else:
        licheeOutputFile.write("\tNormal")
    for sample in samples:
        if sample != healthy:
            licheeOutputFile.write("\t" +sample)
    
    licheeOutputFile.write("\n")
    
    
def printLicheeMut(mut):
    chrom = muts[mut]['chr']
    pos = muts[mut]['pos']
    refAllele = muts[mut]['samplesInfo'][samples[0]]['ref_nuc']
    altAllele = getAltAlleles(mut)[0]
    
    
    licheeOutputFile.write( chrom + "\t" + pos + "\t" + refAllele + "/" + altAllele )
    
    vafs = getPerSampleVaf(mut)
    
    
    if withHealthy():       
        licheeOutputFile.write( "\t" +  str(vafs[healthy]) )
    else:
        licheeOutputFile.write( "\t" +  "0.0" )
        
        
    for sample in samples:
        if sample != healthy:
            licheeOutputFile.write( "\t" +  str(vafs[sample]) )
    
        
   
    licheeOutputFile.write("\n")
    
    

def printLicheeOutput():
    printLicheeHeader()
    for mut in muts.keys():
        printLicheeMut(mut)



counts = open(argv[1], 'r')
licheeOutputFile = open(argv[2], 'w')

if len(argv) == 4:
    healthy=argv[3]
else:
    healthy=""


minDepth = 20
minVaf = 0.04
germlineVaf = 0.1
maxVafIfNotHealthy = 0.9


muts = {}

for linea in counts:
    linea = linea[:-1]
    
    if linea[:3]=="CON":
        samples = parseHeader(linea)
    else:
        getMutsInfo(linea)


for mut in muts.keys():
    
    if checkMultiallelic(mut):
        if keepMultiallelic(mut):
            removeMinorAllele(mut)
        else:
            removeMut(mut)
            continue
    
    if filterLowDepth(mut):
        removeMut(mut)
        continue
        
    
    if filterLowMaxVaf(mut):
        removeMut(mut)
        continue
    
    if withHealthy() and filterGermline(mut):
        removeMut(mut)
        continue
        
    if (not withHealthy()) and filterHighFreq(mut):
        removeMut(mut)
        continue
        
    if isPanelError(mut):
        print mut
        removeMut(mut)
        continue
    
    
printLicheeOutput()
licheeOutputFile.close()
