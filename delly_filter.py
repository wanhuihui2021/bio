import os
import sys

sample_dir=sys.argv[1]
ra_no=sys.argv[2]
bed_fusion_related=sys.argv[3]
bed_probe=sys.argv[4]

def juge_probe(chr1,chr2,site1,site2,probe_dict):
    flag = False
    if chr1 in probe_dict:
        for tup in probe_dict[chr1]:
            if tup[0] <= int(site1) <= tup[1]:
                flag = True
                break
    if flag == False:
        if chr2 in probe_dict:
            for tup in probe_dict[chr2]:
                if tup[0] <= int(site2) <= tup[1]:
                    flag = True
                    break
    return flag


def judge_promoter(gene1,gene2,site1,site2):
    flag=False
    if gene1=="BCL2" and int(site1) >=60986613:
        flag=True
    elif gene1=="BCL6" and int(site1) >=187463256:
        flag=True
    elif gene1=="MYC" and int(site1)<=128747680:
        flag=True
    elif gene2=="BCL2" and int(site2) >=60986613:
        flag=True
    elif gene2=="BCL6" and int(site2) >=187463256:
        flag=True
    elif gene2=="MYC" and int(site2)<=128747680:
        flag=True
    else:
        flag=False
    return flag


print("--------Extracting bed--------")
with open(sample_dir+"/"+ra_no+".bed","w") as bed_output:
    with open(sample_dir+"/"+ra_no+"_translocation.vcf") as f1:
        for line in f1.readlines():
            line=line.strip().split("\t")
            #记录易位信息点1
            bed_output.write(line[0]+"\t"+str(int(line[1])-1)+"\t"+line[1]+"\n")
            #记录易位点2
            if "]" in line[4]:
                info=line[4].strip().split("]")[1].split(":")
                string =info[0]+"\t"+str(int(info[1])-1)+"\t"+\
                         str(info[1]) +"\n"
                bed_output.write(string)
            else:
                info=line[4].strip().split("[")[1].split(":")
                string =info[0]+"\t"+str(int(info[1])-1)+"\t"+\
                         str(info[1]) +"\n"
                bed_output.write(string)
    bed_output.close()

print("--------bedtools annotating--------")
##注释BCL2/BCL6/MYC/IGH/IHL/IGK等基因位置
cmd1="bedtools intersect -a " +sample_dir+"/"+ra_no+".bed" + \
     " -b "+bed_fusion_related+" -wa -wb |cut -f 1,3,7 "+ \
    '>'+sample_dir+"/"+ra_no+"_anno.bed"
os.system(cmd1)

print("--------collecting annotated sites--------")
site_dict={}
white_list=["IGH", "IGL", "IGK", "MYC", "BCL2", "BCL6","CCND1","CCND2","CCND3","IRF4","DUSP22",\
            "TP63","MALT1","ALK"]
with open(sample_dir+"/"+ra_no+"_anno.bed","r") as f3:
    for line in f3.readlines():
        line=line.strip().split("\t")
        site=line[0]+":"+line[1]
        gene_name=line[-1]
        if site in site_dict:
            if site_dict[site] not in white_list:
                site_dict[site] = gene_name
        else:
            site_dict[site] = gene_name
f3.close()

print("--------extracting vcf-----------")
#解析vcf文件提取有用的信息
anno_out=open(sample_dir+"/"+ra_no+"_anno.xls","w")
head="Gene1\tGene2\tBreakpoint1\tBreakpoint2\tReference pairs\t"+\
      "Variant pairs\tReference junction reads\tVariant junction reads\tQua\n"
anno_out.write(head)
with open(sample_dir+"/"+ra_no+"_translocation.vcf","r") as f2:
        for line in f2.readlines():
            line=line.strip().split("\t")
            con=line[7].strip().split(";")[0]
            seq=" "
            if con=="PRECISE":
                seq=line[7].strip().split(";")[-2].split("=")[1]
            qua=con +" | " +seq+ " | "
            counts="\t".join(line[-1].split(":")[-4:])
            if "]" in line[4]:
                site2=line[4].strip().split("]")[1]
                site1=line[0]+":"+line[1]
                if site1 in site_dict:
                    gene1=site_dict[site1]
                else:
                    gene1="Unknown"
                if site2 in site_dict:
                    gene2=site_dict[site2]
                else:
                    gene2="Unknown"
                if judge_promoter(gene1,gene2,site1.split(":")[1],site2.split(":")[1]):
                    pro="Promoter"
                else:
                    pro="UN"
                anno_out.write(gene1+"\t"+gene2+"\t"+\
                          site1+"\t"+site2+"\t"+counts+"\t"+pro+"_"+qua+"\n")
            else:
                site2=line[4].strip().split("[")[1]
                site1=line[0]+":"+line[1]
                if site1 in site_dict:
                    gene1=site_dict[site1]
                else:
                    gene1="Unknown"
                if site2 in site_dict:
                    gene2=site_dict[site2]
                else:
                    gene2="Unknown"
                if judge_promoter(gene1, gene2, site1.split(":")[1], site2.split(":")[1]):
                    pro = "Promoter"
                else:
                    pro = "UN"
                anno_out.write(gene1 + "\t" + gene2 + "\t" + \
                               site1 + "\t" + site2 + "\t" + counts + "\t" + pro + "_" + qua + "\n")
anno_out.close()
f2.close()

print("---------Selecting target genes---------------")
###挑选出一端在IGH/IGK/IGL,另一端在BCL2/BCL6/MYC上的序列
cmd2= "cat "+sample_dir+"/"+ra_no+"_anno.xls | "+ "grep \"MYC\|BCL2\|BCL6\" | grep \"IG\" >"+ \
    sample_dir+"/"+ra_no+"_filtered.xls"
os.system(cmd2)

cmd4="echo -ne \"Gene1\tGene2\tBreakpoint1\tBreakpoint2\tReference pairs\tVariant pairs\tReference junction reads\tVariant junction reads\tQua\n\" >"+sample_dir+"/"+ra_no+"_all.xls"
os.system(cmd4)

cmd3= "cat "+sample_dir+"/"+ra_no+"_anno.xls | "+ "grep \"MYC\|BCL2\|BCL6\|CCND1\|IRF4\|MALT1\|ALK\|DUSP22\|TP63"+\
    "\|CCND2\|CCND3\" |grep -v \"IMPRECISE\|Unknown\" >>"+sample_dir+"/"+ra_no+"_all.xls"
os.system(cmd3)



print("--------Reading probes--------")
def read_probe():
    probe_dict={}
    with open(bed_probe,"r") as f5:
        for line in f5.readlines():
            line=line.strip().split("\t")
            if line[0] not in probe_dict:
                probe_dict[line[0]]=[]
                tup=(int(line[1]),int(line[2]))
                probe_dict[line[0]].append((min(tup),max(tup)))
            else:
                tup=(int(line[1]),int(line[2]))
                probe_dict[line[0]].append((min(tup),max(tup)))
    f5.close()
    return probe_dict



print("--------filtering probes--------")
##过滤探针位点
common_list=["IGH-BCL2","BCL6-IGH","MYC-IGH","BCL2-IGH","IGH-BCL6","IGH-MYC"]
final_out=open(sample_dir+"/"+ra_no+"_final_selected.xls","w")
with open(sample_dir+"/"+ra_no+"_filtered.xls","r") as f4:
    final_out.write(head)
    for line in f4.readlines():
        rawline=line
        line=line.strip().split("\t")
        chr1=line[2].split(":")[0]
        site1=line[2].split(":")[1]
        chr2=line[3].split(":")[0]
        site2=line[3].split(":")[1]
        dict_probe=read_probe()
        com=line[0]+"-"+line[1]
        if juge_probe(chr1,chr2,site1,site2,dict_probe):
            if com in common_list:
                final_out.write(rawline)
f4.close()
final_out.close()
print("---------Delly filter done--------")

