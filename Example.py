############## ALIGN KNOWN WT SEQUENCES ##########

alba='ATGCCGATGCGCGAAATGTTGAGCAGCGAACATGGAAATGCCCTTGCGAGTGTACCAGTGCAGATGAGTCTCTAGCTCGA'
trem='ATGCCGATGAGCGAAATGTTGAG   CGAACATGGAAATGCCCTTGCGAGAGTACCAGTCCAGATGAGTCTCTAGCTCGA'
# SNPs         ^                                        ^        ^
# gRNA                             ACATGGAAATGCCCTTG

# there is a 3 nucleotide indel between the two genes (and a few SNPs)
trem=trem.replace(' ','') #remove the gap. it was for visualization.




########## MERGE THESE ALLELES INTO A SINGLE DEGENERATE SEQUENCE ##########
degen_matrix = [["-","A","C","T","G"],
                ["A","A","M","W","R"],
                ["C","M","C","Y","S"],
                ["T","W","Y","T","K"],
                ["G","R","S","K","G"]]

def merge(seq1,seq2):
    if len(seq2)>len(seq1): #make the first seqence the longer one
        seq1,seq2=seq2,seq1

    temp_seq = ""
    for i in range(len(seq2)): #for the overlap length
        for x in range(len(degen_matrix[0])):
            for y in range(len(degen_matrix)):
                if seq1[i].upper() == degen_matrix[0][x]:
                    if seq2[i].upper() == degen_matrix[y][0]:
                        temp_seq += degen_matrix[y][x]

    temp_seq+=seq1[len(seq2):] #append the overhang
    return temp_seq

degenerate_WT=merge(alba,trem)




######## DEFINE FUNCTIONS FOR SEARCHING ##############

s2_table = {"A":"ARWMN",
            "C":"CYSMN",
            "T":"TYWKN",
            "G":"GRSKN",
            "R":"RN",
            "Y":"YN",
            "M":"MN",
            "K":"KN",
            "W":"WN",
            "S":"SN"}

def equiv_dict(q,dat):
    #can 'q' be found in the degenerate nucleotides in s2_table
    return dat in s2_table[q]

def downsample(q,dat):
    #downsample using just the first nucleotide.
    temp_possible=[]
    for x in range(len(dat)):
        if equiv_dict(q,dat[x]):
            temp_possible.append(x)
    return temp_possible #returns a list of positions

def search(query,data,loud=False):
    possible=downsample(query[0],data[:-len(query)+1])
    if loud:
        print("length of query:  ",len(query))
        print("length of data:   ",len(data))
    matches=list()
    for l in range(len(possible)):              #for each potential site
        for q_index in range(len(query)):       #for the length of the query
            if equiv_dict(query[q_index],data[possible[l]+q_index]): #if query can degenerately match data (index by index)
                if q_index+1 == len(query):
                    if loud:
                        print("found at : ", possible[l])
                        print("which is : ", data[possible[l]:possible[l]+q_index+1])
                        print("against q: ", query)
                    matches.append(possible[l])
            else:
                break
    if loud:
        print(len(possible), " sites attempted. (as downsampled by first nucleotide)")
    return matches


############### SELECT QUERY SEQUENCES ##############

#alba='ATGCCGATGCGCGAAATGTTGAGCAGCGAACATGGAAATGCCCTTGCGAGTGTACCAGTGCAGATGAGTCTCTAGCTCGA'
#trem='ATGCCGATGAGCGAAATGTTGAG   CGAACATGGAAATGCCCTTGCGAGAGTACCAGTCCAGATGAGTCTCTAGCTCGA'
# SNPs          ^                                        ^        ^
# gRNA                              ACATGGAAATGCCCTTG
# pre_grna                       CGAACA
# post_grna                                                            TGAGTCTCTAGCTCGA



############## GET WT REFERENCE ############

pre_grna='CGAACA'           # Query before gRNA cut site
post_grna='TGAGTCTCTAGCTCGA' # Query after gRNA cut site
# These were selected to avoid SNPs and to be far enough away from the cut site that they won't be affected.

# Find locations of these query sequences in the degenerate_WT sequence. (Ideally you'll get two instances each).
pre_grna_locs=search(pre_grna,degenerate_WT)
post_grna_locs=search(post_grna,degenerate_WT)

#The distance between the pre and post query sequences in the WT
WT_spacing=list()
WT_spacing.append(post_grna_locs[0]-pre_grna_locs[0])
WT_spacing.append(post_grna_locs[1]-pre_grna_locs[1])

print(WT_spacing) #Should be [35,35]




############ GET SAMPLE_1 SPACING ############

# This is from sequencing an unknown read near the indel/gRNA
simulated_indel='GCATGCCGATGMGCGAAATGTTGAGCRRMSAWSRWRRWRMYSYYSKYGAKMGWRYSWRYCMRKRYRRRTSWSTMKCTMGMTCGA'

#Using the same Query sequences as above.
#pre_grna='CGAACA'
#post_grna='TGAGTCTCTAGCTCGA'
pre_grna_locs=search(pre_grna,simulated_indel)
post_grna_locs=search(post_grna,simulated_indel)

# Get the distance between these Query sequences
sample1_spacing=list()
sample1_spacing.append(post_grna_locs[0]-pre_grna_locs[0])
sample1_spacing.append(post_grna_locs[1]-pre_grna_locs[1])




########### COMPARE SIMULATED SPACING TO WT #########

print('WT:',WT_spacing)
print('S1:',sample1_spacing)
# shows there is a +1/+2 insertion from the cas9 cut.
# since both will induce a frameshift, this sample should be looked at more closely.


