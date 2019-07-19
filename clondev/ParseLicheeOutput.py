# Python 2.7

from sys import argv

if len(argv)!= 3:
    print "Usage "+argv[0]+" <Input: path/to/Lichee/output.Lichee> <Output: /path/to/output.fasta"
    
lichee = open(argv[1], 'r')
output = open(argv[2], 'w')
nodes=False
tree=False
snvinfosection=False

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
        snvsInNodes[nodeNumber] = snvIds
        
    elif(tree):
        linea = linea[:-1].split(" -> ")
        if linea[0]=="0":
            continue
        snvsInNodes[linea[1]]+=snvsInNodes[linea[0]]
    
    elif(snvinfosection):
        linea= linea[:-1].split(" ")
        snvId = linea[0][:-1]
        snvsInfo[snvId] = linea[1:]
        
        
        
    if linea[:5]=="Nodes":
        nodes=True
    elif linea[:14]=="****Tree 0****":
        tree=True
        nodes=False
    elif linea[:9]=="SNV info:":
        snvinfosection = True
        
    
output.write(">healthy\n")

for snv in snvsInfo.keys():
    output.write(snvsInfo[snv][2].split("/")[0])
output.write("\n")


for clon in sorted(snvsInNodes.keys()):
    output.write(">clon_"+str(clon)+"\n")
    for snv in snvsInfo.keys():
        if snv in snvsInNodes[clon]:
            output.write(snvsInfo[snv][2].split("/")[1])
        else:
            output.write(snvsInfo[snv][2].split("/")[0])
    output.write("\n")
            
