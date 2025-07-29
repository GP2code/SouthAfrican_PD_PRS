path = "full_path/msp_files"
# file name for msp with chromosome number replaced by "CHROM"
fileName = "query_results_CHROM.msp"

outputFile = open(f"{path}/proportionANCperIndividual.tsv", 'w')
outputFile.write("IID")

dictInds = {}
# number parental populations in file
numberParentals = 5
# setting total length of windows to 0
sizeAllWindows = 0
headerPrint = False
# iterate through chromsomes
for i in range(1,3):
    print(f"Reading chromosome {i}")
    mspFile = f"{path}/{fileName.replace('CHROM',str(i))}"
    with open(mspFile,'r') as msp:
        header = True
        for line in msp:
            if header:
                if headerPrint == False:
                    # extracting which subpopulation goes with which number
                    if line.startswith("#Subpopulation order/codes:"):
                        split = line.strip().split('\t')
                        for j in range(0, numberParentals):
                            # write subpopulations in numerical order
                            if split[j].endswith(str(j)):
                                if j == 0:
                                    # replacing beginning text for first subpopulation
                                    outputFile.write(f"\t{split[j].replace('#Subpopulation order/codes: ','').split('=')[0]}")
                                    print(split[j].replace('#Subpopulation order/codes: ',''))0
                                else:
                                    outputFile.write(f"\t{split[j].split('=')[0]}")
                                    print(split[j])
                            else:
                                # printing error message if the subpopulations are out of order
                                print(f"Problem with {split[j]}, may be out of numerical order")
                        outputFile.write("\n")
                        headerPrint = True
                # saving header values
                if line.startswith("#chm"):
                    headerLine = line
                    split = line.strip().split('\t')
                    for j in range(6,len(split)):
                        # combining haplotypes into one ID
                        indID = split[j].replace(".0","").replace(".1","")
                        if indID not in dictInds:
                            dictInds[indID] = {}
                            # setting length of cM covered by ancestry to 0
                            for anc in range(0, numberParentals):
                                dictInds[indID][anc] = 0.0
                    header = False
                    print("Done reading header")
            else:
                split = line.strip().split('\t')
                # calculating exact size of each window
                sizeWindow = float(split[4])-float(split[3])
                # adding length of window to total ancestral window length for individuals
                for j in range(6,len(headerLine.split('\t'))):
                    indID = headerLine.strip().split('\t')[j].replace(".0","").replace(".1","")
                    anc = int(split[j])
                    dictInds[indID][anc] = dictInds[indID][anc]+sizeWindow
                # adding length of window to total genomic window length
                sizeAllWindows += sizeWindow
    print("Done extracting windows")

for ind in dictInds:
    outputFile.write(ind)
    # calulating average contribution of each ancestry to each individual and adding to output file
    for anc in range(0, numberParentals):
        avg = (dictInds[ind][anc]/sizeAllWindows)/2
        outputFile.write(f"\t{str(avg)}")
    outputFile.write("\n")
print("Done writing file")

# closing output file
outputFile.close()