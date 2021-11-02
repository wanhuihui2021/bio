import os
import sys
white_list=["ALK","ROS1","RET","NTRK1","NTRK3","NTRK2","NRG1", "FGFR2","FGFR3"]
head="Gene 1\tGene 2\tBreakpoint 1\tBreakpoint 2\tTranscript 1\tTranscript 2\tTotal reads\tUnique reads"
out=open(sys.argv[1][0:-4]+"_selected.xls","w")
select_dict={}
with open(sys.argv[1],"r")as input_file:
    for line in input_file.readlines():
        if line.startswith("#Fusion"):
            line=line.strip().split("_")
            gene1=line[0].split(":")[1].strip()
            gene2=line[4]
            transcript1=line[1]
            transcript2=line[-1].strip().split()[0]
            if "|" in line[1]:
                location1=line[1].split("|")[1]
            else:
                location1=line[1].split(":")[-2]+":"+line[1].split(":")[-1]
            if "|" in line[-1].split()[0]:
                location2=line[-1].split()[0].split("|")[1]
            else:
                location2=line[-1].split()[0].split(":")[-2] +":" +\
                           line[-1].split()[0].split(":")[-1]
            total_reads=line[-1].split()[-2].split(",")[0]
            uniq_reads=line[-1].split()[-1].split(":")[-1].split(")")[0]
            if (gene1 in white_list) or (gene2 in white_list):
                info_list=[gene1,gene2,location1[1:],location2[1:],transcript1,\
                           transcript2,total_reads,uniq_reads]
                if total_reads not in select_dict:
                    select_dict[total_reads]=[]
                    select_dict[total_reads].append(info_list)
                else:
                    select_dict[total_reads].append(info_list)
    if select_dict !={}:
        out.write(head+"\n")
        final_key=max(map(eval,list(select_dict.keys())))
        final_key=str(final_key)
        for i in range(len(select_dict[final_key])):
            string="\t".join(select_dict[final_key][i])
            out.write(string+"\n")
    else:
        out.write(head+"\n")
out.close()
input_file.close()
