#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Removes reciprocal hits (A-B = B-A). Also adds together multiple hits while considering overlaps")

# Add the arguments to the parser
parser.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input file. This assumes a file with the following columns: 'qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen'.")
parser.add_argument("-o", "--output", dest="file_out", required=True,
                    help="Output file. Returns the filtered file to the specified location.")
parser.add_argument("-m","--mode",default="size", dest = "mode",
                    help = "mode of priority for counting hits, are size of hit or percent identity prioritized when calculating overlapping hits? Options: 'size' and 'identity'")
args = parser.parse_args()

# This function simply looks for an overlap between two ranges
def overlap(start1, end1, start2, end2):
    y = range(start1,end1)
    x = range(start2, end2) 

    """how much does the range (start1, end1) overlap with (start2, end2)"""
    r = range(max(x[0],y[0]),min(x[-1],y[-1]+1))
    return len(r)



print("")
print("  Setting hit names")


dirty = {}
clean = {}
removed = 0
accepted = 0
selfHit = 0

# This loops through the input file
for line in open(args.file_in):
    line = line.strip().split()
    seq1 = line[0]
    seq2 = line[1]
    sstart = min(int(line[8]),int(line[9]))
    send = max(int(line[8]),int(line[9]))
    qstart = min(int(line[5]),int(line[6]))
    qend = max(int(line[5]),int(line[6]))

    #Removes selfhits from the final output
    if seq1 == seq2:
        selfHit = selfHit + 1

    else:
        hit = str(seq1) + "\t" + str(seq2)
        hitr = str(seq2) + "\t" + str(seq1)

        #removes reciprocal hit
        if hitr in dirty:
            removed = removed + 1 

        # If the hit is part of multiple fragmented hits they are added in the list of lists of the hit dictionary key
        elif hit in dirty:

            #Removes small hits
            if max(int(line[5]),int(line[6]))-min(int(line[5]),int(line[6])) < 40 or max(int(line[8]),int(line[9]))-min(int(line[8]),int(line[9])) < 40:
                removed = removed +1
            else:
                ssize = max(int(line[9]),int(line[8]))-min(int(line[8]),int(line[9]))
                qsize = max(int(line[6]), int(line[5]))-min(int(line[5]),int(line[6]))
                
                values = [line[2],float(line[3]),line[4],int(line[5]),int(line[6]),int(line[7]),int(line[8]),int(line[9]),int(line[10]), int((float(line[3])/100)*ssize),qsize,ssize]

                dirty[hit].append(values)
        
        #If this is the first time a hit between these two sequences are encountered, then it creates a key of the dictionary dirty and adds it as a list within a list there
        else:
            if max(int(line[5]),int(line[6]))-min(int(line[5]),int(line[6])) < 40 or max(int(line[8]),int(line[9]))-min(int(line[8]),int(line[9])) < 40:
                removed = removed +1
            else:   
                
                ssize = max(int(line[9]),int(line[8]))-min(int(line[8]),int(line[9]))
                qsize = max(int(line[6]), int(line[5]))-min(int(line[5]),int(line[6]))
                
                values = [line[2],float(line[3]),line[4],int(line[5]),int(line[6]),int(line[7]),int(line[8]),int(line[9]),int(line[10]), int((float(line[3])/100)*ssize),qsize,ssize]

                dirty[hit] = [values]

# This big loop goes through the hits in the dirty dictionary, cleaning them up before outputting them
for hit in dirty:

    ## If there is only one hit between sequence A and B, it is immediately put in the clean dictionary
    if len(dirty[hit]) == 1:

        clean[hit] = dirty[hit][0]
        accepted=accepted+1
        
    
    else:

        #Starts to loop through the different hits for each pair of sequences
        for i, list in enumerate(dirty[hit]):
            
            
            #Calculates the number of matching bps between the hits based on the percent identity and length
            list_percBP = (float(list[1])/100)*list[11]

            #Loops through all other hits....
            for j, list2 in enumerate(dirty[hit]):
                
                #.. Except for the current hit
                if list == list2:
                    pass
                
                else:
                    
                    #Checks if there is overlap between the hits in the subject sequence and overlap between query
                    if overlap(min(list[6],list[7]), max(list[6],list[7]),min(list2[6],list2[7]), max(list2[6],list2[7])) != 0 and overlap(min(list[3],list[4]), max(list[3],list[4]),min(list2[3],list2[4]), max(list2[3],list2[4])) != 0: 
                        
                        #calculates the overlap between the two hits in the subject
                        cur_sover = overlap(min(list[6],list[7]), max(list[6],list[7]),min(list2[6],list2[7]), max(list2[6],list2[7]))
                        
                        # This mode considers if a smaller hit is better and then removes the overlap from the first hit instead
                        if args.mode == "identity":

                            #Length of the second hit, and minus the overlap
                            list2_hit = max(list2[6],list2[7])-min(list2[6],list2[7])
                            
                            #removes the overlap from list1 instead when recounting the perc_identity
                            list1_hit_nooverlap = dirty[hit][i][11]-cur_sover
                            list1_hit_perc_BP = (float(list[1])/100)*list1_hit_nooverlap

                            # Nr of BPs that match between query and subject in hit 2, calculated by dividing per_identity with length of hit. 
                            # This also takes into account overlap with first hit to not double count.
                            list2_percBP = (float(list2[1])/100)*list2_hit

                            # This adds together the hits and recalculates the perc_identity of the combined hits and updates hit 1 with the new perc_identity and total hit length
                            new_perc_id = (list1_hit_perc_BP+list2_percBP)/(list1_hit_nooverlap+list2_hit) 

                            dirty[hit][i][1] = round(new_perc_id*100,3)
                            dirty[hit][i][11] += list2_hit

                        #This mode will always remove the overlap from the second hit
                        elif args.mode == "size":
                            
                            #Length of the second hit, and minus the overlap
                            list2_hit = max(list2[6],list2[7])-min(list2[6],list2[7])
                            list2_hit_nooverlap = list2_hit-cur_sover

                            # Nr of BPs that match between query and subject in hit 2, calculated by dividing per_identity with length of hit. 
                            # This also takes into account overlap with first hit to not double count.
                            list2_percBP = (float(list2[1])/100)*list2_hit_nooverlap

                            # This adds together the hits and recalculates the perc_identity of the combined hits and updates hit 1 with the new perc_identity and total hit length
                            new_perc_id = (list_percBP+list2_percBP)/(list[11]+list2_hit_nooverlap) 

                            dirty[hit][i][1] = round(new_perc_id*100,3)
                            dirty[hit][i][11] += list2_hit_nooverlap

                        # Updates the range of the hit to be able to compare the new range with further hits
                        # Range is not updated for non-overlapping reads so this is not reflecting the size of all combined reads, it is only useful for checking overlaps
                        dirty[hit][i][6] = min((list[6]),(list[7]),(list2[6]),(list2[7])) 
                        dirty[hit][i][7] = max((list[6]),(list[7]),(list2[6]),(list2[7]))

                        #Doing the same for query overlap, this is only important for coverage of hit in the query sequence
                        cur_qover  = overlap(min(list[3],list[4]), max(list[3],list[4]),min(list2[3],list2[4]), max(list2[3],list2[4]))
                            
                        list2_qhit = max(list2[3],list2[4])-min(list2[3],list2[4])
                        list2_qhit_nooverlap = list2_qhit-cur_qover

                        dirty[hit][i][10] += list2_qhit_nooverlap


                        removed = removed +1
                        dirty[hit].remove(list2)


                    #Checks if there is overlap only between subject hits
                    elif overlap(min(list[6],list[7]), max(list[6],list[7]),min(list2[6],list2[7]), max(list2[6],list2[7])) != 0 and overlap(min(list[3],list[4]), max(list[3],list[4]),min(list2[3],list2[4]), max(list2[3],list2[4])) == 0: 
                        
                        #calculates the overlap between the two hits in the subject
                        cur_sover = overlap(min(list[6],list[7]), max(list[6],list[7]),min(list2[6],list2[7]), max(list2[6],list2[7]))
                        
                        
                        if args.mode == "identity":

                            #Length of the second hit, and minus the overlap
                            list2_hit = max(list2[6],list2[7])-min(list2[6],list2[7])
                            
                            #removes the overlap from list1 instead when recounting the perc_identity
                            list1_hit_nooverlap = dirty[hit][i][11]-cur_sover
                            list1_hit_perc_BP = (float(list[1])/100)*list1_hit_nooverlap

                            # Nr of BPs that match between query and subject in hit 2, calculated by dividing per_identity with length of hit. 
                            # This also takes into account overlap with first hit to not double count.
                            list2_percBP = (float(list2[1])/100)*list2_hit

                            # This adds together the hits and recalculates the perc_identity of the combined hits and updates hit 1 with the new perc_identity and total hit length
                            new_perc_id = (list1_hit_perc_BP+list2_percBP)/(list1_hit_nooverlap+list2_hit) 

                            dirty[hit][i][1] = round(new_perc_id*100,3)
                            dirty[hit][i][11] += list2_hit

                        elif args.mode == "size":
                            
                            #Length of the second hit, and minus the overlap
                            list2_hit = max(list2[6],list2[7])-min(list2[6],list2[7])
                            list2_hit_nooverlap = list2_hit-cur_sover

                            # Nr of BPs that match between query and subject in hit 2, calculated by dividing per_identity with length of hit. 
                            # This also takes into account overlap with first hit to not double count.
                            list2_percBP = (float(list2[1])/100)*list2_hit_nooverlap

                            # This adds together the hits and recalculates the perc_identity of the combined hits and updates hit 1 with the new perc_identity and total hit length
                            new_perc_id = (list_percBP+list2_percBP)/(list[11]+list2_hit_nooverlap) 

                            dirty[hit][i][1] = round(new_perc_id*100,3)
                            dirty[hit][i][11] += list2_hit_nooverlap



                        # Updates the range of the hit to be able to compare the new range with further hits
                        # Range is not updated for non-overlapping reads so this is not reflecting the size of all combined reads, it is only useful for checking overlaps
                        dirty[hit][i][6] = min((list[6]),(list[7]),(list2[6]),(list2[7])) 
                        dirty[hit][i][7] = max((list[6]),(list[7]),(list2[6]),(list2[7]))

                        dirty[hit][i][10] += (max(list2[3],list2[4])-min(list2[3],list2[4]))
                        
                        removed = removed +1
                        dirty[hit].remove(list2)
                    
                    # Removes tandem repeat expansions / duplications where one sequence in the query has multiple hits in the subject
                    elif min(list[6],list[7]) == min(list[6],list[7]) and max(list[6],list[7]) == max(list2[6],list2[7]):
                        
                        removed = removed+1
                        dirty[hit].remove(list2)


        # Goes through each non-overlapping hit, adds the non-overlapping hit to the first one in the same way as above        
        n=0
        while len(dirty[hit]) > 1: 
            if n==0:
                n+=1
            else:
                
                #Removes potential tandem repeats from the counting of hit (If the query hit is the same but not the subject)
                if max(dirty[hit][0][3], dirty[hit][0][4]) == max(dirty[hit][n][3], dirty[hit][n][4]) and min(dirty[hit][0][3], dirty[hit][0][4]) == min(dirty[hit][n][3], dirty[hit][n][4]):
                    removed = removed+1
                    dirty[hit].pop(n)

                else:
                    cur_hit_len = max(dirty[hit][n][7],dirty[hit][n][6])-min(dirty[hit][n][7],dirty[hit][n][6])
                    pre_hit_percBP = (float(dirty[hit][0][1])/100)*dirty[hit][0][11]
                    cur_hit_percBP = (float(dirty[hit][n][1])/100)*dirty[hit][n][11]

                    qover = overlap(min(dirty[hit][0][3], dirty[hit][n][4]),max(dirty[hit][0][3], dirty[hit][0][4]),min(dirty[hit][n][3], dirty[hit][n][4]),max(dirty[hit][n][3], dirty[hit][n][4]))
            

                    cur_qlen_nooverlap = (max(dirty[hit][n][3], dirty[hit][n][4])-min(dirty[hit][n][3], dirty[hit][n][4]))-qover
                    perc_id = (pre_hit_percBP+cur_hit_percBP)/(dirty[hit][0][11]+dirty[hit][n][11])

                    dirty[hit][0][1] = round(perc_id*100,3)
                    dirty[hit][0][10] += cur_qlen_nooverlap
                    dirty[hit][0][11]+=cur_hit_len

                    removed = removed+1
                    dirty[hit].pop(n)
        
        
        clean[hit] = dirty[hit][0]
        accepted = accepted+1


#Outputs the clean dictionary, which should now only have one combined hit per sequence pair
print("  Writing cleaned file with '", accepted, "' hits", sep="")

with open(args.file_out, "w") as outfile:
    for hit in clean.keys():
        print(hit, *clean[hit], sep='\t', file=outfile)
        

print("    Removed hits:\t", removed)
print("    Self matching:\t", selfHit)
print("Done")
print("")
