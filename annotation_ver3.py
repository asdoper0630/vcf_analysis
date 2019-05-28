vcf_fn =input('vcf filename:')

with open(vcf_fn,'r') as vcf:
    vcf_lines = vcf.readlines()
with open('Homo_sapiens.GRCh38.96.chr.gff3') as gff:
    gff_lines = gff.readlines()

chr_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'X', 'Y'] # for문 돌리는 기준

gff_chposid_list = []

for gff_line in gff_lines:
    glsp = gff_line.split()

    if len(glsp) < 6:
        pass
    elif glsp[2] == 'mRNA' and glsp[6] == '+' :
            gff_chposid_list.append((glsp[0],int(glsp[3]),int(glsp[4]),glsp[8].split(';')[2].split('-')[0][5:]))

gff_dict = {}

for gff_chposid in gff_chposid_list:
    if gff_chposid[0] in gff_dict:
        gff_dict[gff_chposid[0]].append((gff_chposid[1],gff_chposid[2],gff_chposid[3]))
    elif gff_chposid[0] not in gff_dict:
        gff_dict[gff_chposid[0]] = [(gff_chposid[1],gff_chposid[2],gff_chposid[3])]
    else:
        print('error')

vcf_chpos_dic = {}

for vcf_line in vcf_lines[:-1]:
    vlsp = vcf_line.split()

    if vlsp[0][0] != '#':
        if vlsp[0] not in vcf_chpos_dic:
            vcf_chpos_dic[vlsp[0]] = [vcf_line.strip()]
        elif vlsp[0] in vcf_chpos_dic:
            vcf_chpos_dic[vlsp[0]].append(vcf_line.strip())


result = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTERINFO\tGENE_ID\n'

found = False

for chrom in chrom_list: # 1~23~MTXY
    if chrom not in vcf_chpos_dic:
        continue
    for vcf in vcf_chpos_dic[chrom]:
        vcf_num = int(vcf.split()[1])
        for gff in gff_dict[chrom]:
            stt_num = gff[0]
            end_num = gff[1]
            gene_ID = gff[2]
            if vcf_num >= stt_num and vcf_num <= end_num:
                result += vcf + '\t' + gene_ID + '\n'
                found = True
        if found == False:
            result += vcf + '\t.\n'
with open(vcf_fn[:-4]+'_GeneID_added'+vcf_fn[-4:]+'.txt','w') as sv:
    sv.write(result)
