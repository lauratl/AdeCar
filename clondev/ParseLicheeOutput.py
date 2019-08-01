# Python 2.7

from sys import argv

if len(argv)!= 3:
    print "Usage "+argv[0]+" <Input: path/to/Lichee/output.Lichee> <Output: /path/to/output.fasta"
    
lichee = open(argv[1], 'r')
output = open(argv[2], 'w')
nodes=False
tree=False
snvinfosection=False
clonesToClusters=False
snvsInNodes = {}
snvsInfo = {}

for linea in lichee:
    if linea[:-1] =="":
        nodes=False
    elif linea[:12]=="Error score:":
        tree=False 
	
        
    elif(nodes):
	
        linea=linea[:-1]
        linea = linea.split("]")
        nodeNumber = linea[0].split("\t")[0]
        snvIds = linea[1].split("\t")[1:]
	snvNumbers = []
	for snvId in snvIds:
		snvNumbers.append(int(snvId[3:]))
      				
        snvsInNodes[nodeNumber] = snvNumbers
        
    elif(tree):
        linea = linea[:-1].split(" -> ")
        #if linea[0]=="0":
        #    continue
        topology_pairs.append([linea[0],linea[1]])
    
    elif(snvinfosection):
        linea= linea[:-1].split(" ")
        snvId = int(linea[0][:-1][3:])
        snvsInfo[snvId] = linea[1:]
    
    elif(clonesToClusters):
	
	new_origins = ["0"]
        while True:
            origins = new_origins
	    new_origins = []
	    for origin in origins:

	        for topology_pair in topology_pairs:
 	            if topology_pair[0]==origin:
			
		
		        if origin != "0":
			    snvsInNodes[topology_pair[1]]+=snvsInNodes[topology_pair[0]]
                        new_origins.append(topology_pair[1])
	    if len(new_origins) == 0:
                break

        


	clonesToClusters = False 
        
    if linea[:5]=="Nodes":
        nodes=True


    elif linea[:12]=="Error score:":
	clonesToClusters=True

    elif linea[:14]=="****Tree 0****":
        tree=True
        nodes=False

	topology = {}
	for node in snvsInNodes.keys():
		topology[node] = [node]

	topology_pairs = []
    elif linea[:9]=="SNV info:":
        snvinfosection = True



    
output.write(">healthy\n")

for snv in snvsInfo.keys():
    output.write(snvsInfo[snv][2].split("/")[0])
output.write("\n")


for clon in sorted(snvsInNodes.keys()):
    output.write(">clon_"+str(clon)+"\n")
    for snv in sorted(snvsInfo.keys()):
	if snv in snvsInNodes[clon]:
            output.write(snvsInfo[snv][2].split("/")[1])
        else:
            output.write(snvsInfo[snv][2].split("/")[0])
    output.write("\n")
            
