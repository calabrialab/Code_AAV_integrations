import sys

import pysam
import pandas as pd
import argparse


def getListFromCigar2(cigar):
    pos = 0
    cigar_listf = []
    for i in range(len(cigar)):
        if cigar[i].isalpha():
            cigar_listf.append((int(cigar[pos:i]), cigar[i]))
            pos = i + 1
    return cigar_listf

def getListFromCigar(cigar):
    pos = 0
    cigar_listf = []
    for i in range(len(cigar)):
        if cigar[i].isalpha():
            if cigar[i]!='I' and cigar[i]!='D' :
                cigar_listf.append((int(cigar[pos:i]), cigar[i]))
            elif cigar[i]=='D':
                cigar_listf.append((int(cigar[pos:i]),'M'))
            pos = i + 1
    cleaned_cigar_list = []
    cleaned_cigar_list.append(cigar_listf[0])
    j=0
    for i in range(1,len(cigar_listf)):
        if cigar_listf[i][1]==cigar_listf[i-1][1]:
            matches = cleaned_cigar_list[j][0] + cigar_listf[i][0]
            cleaned_cigar_list[j] = (matches, cigar_listf[i][1])
        else:
            cleaned_cigar_list.append(cigar_listf[i])
            j += 1

    return cleaned_cigar_list

def getStartingMismatch(cigar_listf):
    no_match = 0
    for i in range(len(cigar_listf)):
        bases, cigar_type = cigar_listf[i]
        if cigar_type != 'M':
            no_match += bases
        else:
            return (no_match, i)


def getClosingMismatch(cigar_listf):
    no_match, pos = getStartingMismatch(cigar_listf[::-1])
    return no_match, len(cigar_listf) - pos


def checkChrPos(dict_list, chr_to_check, pos_to_check):
    for i_chr, i_pos in dict_list:
        if i_chr == chr_to_check and -3000 < pos_to_check - i_pos < 3000:
            return True
    return False

def getMatchesFromCigarString(cigar_str):
    cigar_list = getListFromCigar(cigar_str)
    matches = [int(x[0]) for x in cigar_list if x[1]=='M']
    return sum(matches)
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input bam file')
parser.add_argument('-o', '--output', help='path/basename for outputs')
parser.add_argument('-F', '--input_F4', help='input F4 bam file')
parser.add_argument('-l', '--lower_gap', help='maximum negative gap accepted')
parser.add_argument('-g', '--higher_gap', help='maximum positive gap accepted')
# parser.add_arxgument('-R', '--input_F320', help='input F320 bam file')
args = parser.parse_args()
file = args.input
file_F4 = args.input_F4
# file_F320 = args.input_F320
output = args.output
max_gap = int(args.higher_gap)
min_gap = -int(args.lower_gap)
aav_gap_threshold = 50

dict_r2 = {}
r1_step3 = {}
r2_step3 = {}
r2_samfile = pysam.AlignmentFile(file_F4, "rb")
for read in r2_samfile.fetch():
    if read.is_read2:
        if read.query_name in dict_r2:
            dict_r2[read.query_name].append([read.reference_name, read.reference_start])
            r2_step3[read.query_name].append(read)
        else:
            dict_r2[read.query_name] = [[read.reference_name, read.reference_start]]
            r2_step3[read.query_name] = [read]
    else:
        if read.query_name in r1_step3:
            r1_step3[read.query_name].append(read)
        else:
            r1_step3[read.query_name] = [read]

samfile = pysam.AlignmentFile(file, "rb")
list_results = []
list_only_rear = []
list_no_integration = [[], []]
r1_counter = 0
r1_noIS_counter = 0
r1_onlyaav_counter = 0
r1_except_counter = 0
# pairedreads = pysam.AlignmentFile("allpaired.bam", "wb", template=samfile)
for read in samfile.fetch():
    if read.is_read1:
        r1_counter += 1
        try:
            if read.has_tag("SA"):
                sa = read.get_tag("SA")
                all_alignments = [nl.split(",") for nl in sa.split(";")[:-1]]
                strand = "-" if read.is_reverse else "+"
                all_alignments.append(
                    [read.reference_name, read.reference_start, strand, read.cigarstring, read.mapq, 0])
                chrV_alignments = []
                others_alignments = []
                for i in range(len(all_alignments)):
                    if all_alignments[i][0] == "chrV":
                        chrV_alignments.append(all_alignments[i])
                    else:
                        others_alignments.append(all_alignments[i])
                results = {}
                for x_oth in others_alignments:
                    x_oth.append(getMatchesFromCigarString(x_oth[3]))
                starts_chrV = []
                chrV_cigars_list = []
                for i in range(len(chrV_alignments)):
                    cigar_list = getListFromCigar(chrV_alignments[i][3])
                    chrV_cigars_list.append(cigar_list)
                    c = getStartingMismatch(cigar_list) if chrV_alignments[i][2] == '+' else getClosingMismatch(
                        cigar_list)
                    if chrV_alignments[i][2] == '-':
                        c = (c[0], c[1] - 1) # adjusting the pos if reverse strand as explained above
                    starts_chrV.append((c, i))

                starts_chrV.sort(key=lambda tup: tup[0][0])
                start, index = starts_chrV[0][0]
                elem = starts_chrV[0][1]
                end = start + chrV_cigars_list[elem][index][0]
                other_chrV_aln = len(chrV_alignments) - 1
                results['name'] = read.query_name
                results['n_aav_aln'] = 1
                results['start_aav'] = start
                results['aav_matches'] = end-start
                aav_string = ",".join([str(x) for x in chrV_alignments[elem]])
                aav_string += ";"
                last_chrV_alignment = chrV_alignments[elem]
                last_matches = chrV_cigars_list[elem][index][0]
                new_aav_string =",".join(str(x) for x in ["chrV",chrV_alignments[elem][1],int(chrV_alignments[elem][1])+int(last_matches),chrV_alignments[elem][2]])

                for i in range(1, other_chrV_aln + 1):
                    start, index = starts_chrV[i][0]
                    elem = starts_chrV[i][1]
                    aav_gap = start - end
                    if abs(aav_gap) <= aav_gap_threshold:
                        end = start + chrV_cigars_list[elem][index][0]
                        results['n_aav_aln'] += 1
                        aav_string += ",".join([str(x) for x in chrV_alignments[elem]])
                        aav_string += ";"

                        last_chrV_alignment = chrV_alignments[elem]
                        last_matches = chrV_cigars_list[elem][index][0]
                        new_aav_string +=";"+ ",".join(str(x) for x in ["chrV", chrV_alignments[elem][1], int(chrV_alignments[elem][1]) + int(last_matches),
                             chrV_alignments[elem][2]])

                        results['aav_matches'] += chrV_cigars_list[elem][index][0]
                        if aav_gap<0:
                            results['aav_matches'] += aav_gap

                results['input_file']=file_F4
                results['aav_alignments'] = aav_string
                results['aav_alignments_start_end'] = new_aav_string
                results['end_aav'] = end
                results['aav_last_start'] = int(last_chrV_alignment[1])
                results['aav_last_end'] = int(last_chrV_alignment[1])+int(last_matches)
                results['aav_last_strand'] = last_chrV_alignment[2]
                if last_chrV_alignment[2]=='+':
                    results['junction_locus'] = results['aav_last_end']
                else:
                    results['junction_locus'] = results['aav_last_start']
                starts_others = []
                others_cigars_list = []
                results['gap'] = None
                results['other'] = None
                others_alignments.sort(key=lambda x: int(x[-1]), reverse=True)
                for i in range(len(others_alignments)):
                    cigar_list = getListFromCigar(others_alignments[i][3])
                    others_cigars_list.append(cigar_list)
                    start, index = getStartingMismatch(cigar_list) if others_alignments[i][
                                                                          2] == '+' else getClosingMismatch(cigar_list)
                    others_alignments[i].append(start - end)

                #others_alignments.sort(key=lambda x: abs(int(x[-1])), reverse=False)

                dict_r2_chr_list = [c[0] for c in dict_r2[read.query_name]]
                for i in range(len(others_alignments)):
                    # cigar_list = getListFromCigar(others_alignments[i][3])
                    # others_cigars_list.append(cigar_list)
                    # start, index = getStartingMismatch(cigar_list) if others_alignments[i][2] == '+' else getClosingMismatch(cigar_list)
                    gap = others_alignments[i][-1]
                    if max_gap >= gap >= min_gap and others_alignments[i][0] in dict_r2_chr_list and checkChrPos(
                            dict_r2[read.query_name], others_alignments[i][0], int(others_alignments[i][1])):
                        results['gap'] = gap
                        results['other'] = others_alignments[i]
                        results['target_chr'] = others_alignments[i][0]
                        results['target_start'] = others_alignments[i][1]
                        cigar_list = getListFromCigar(others_alignments[i][3])
                        start, index = getStartingMismatch(cigar_list) if others_alignments[i][
                                                                              2] == '+' else getClosingMismatch(
                            cigar_list)
                        if others_alignments[i][2]=='-':
                            index -= 1
                        results['target_end'] = int(others_alignments[i][1]) + int(cigar_list[index][0])
                        results['target_strand'] = others_alignments[i][2]
                        results['target_cigar'] = others_alignments[i][3]
                        results['target_matches'] = int(cigar_list[index][0])
                        results['total_matches'] = results['aav_matches'] + results['target_matches']
                        results['from_junc_to_plusN'] = None
                        results['from_junc_to_plus20'] = None
                        results['from_minusN_to_IS'] = None
                        original_sequence = read.get_forward_sequence()
                        if others_alignments[i][2]=='+':
                            results['integration_locus'] = results['target_start']
                        else:
                            results['integration_locus'] = results['target_end']
                        if results['gap']>=0:
                            results['seq_gap'] = original_sequence[results['end_aav']:results['end_aav'] + results['gap']]
                            results['from_junc_to_plusN'] = original_sequence[results['end_aav']:results['end_aav'] + 30]
                            results['from_junc_to_plus20'] = original_sequence[results['end_aav']:results['end_aav'] + 20]
                            loc_x = results['end_aav']+results['gap']
                            results['from_minusN_to_IS'] = original_sequence[loc_x-12:loc_x]
                        else:
                            results['total_matches'] += results['gap']
                            results['seq_gap'] = original_sequence[results['end_aav'] + results['gap']:results['end_aav']]
                            loc_x = results['end_aav'] + results['gap']
                            results['from_junc_to_plusN'] = original_sequence[loc_x:loc_x + 30]
                            loc_x = results['end_aav'] + results['gap']
                            results['from_minusN_to_IS'] = original_sequence[results['end_aav'] - 12:results['end_aav']]
                        results['integration'] = ":".join(str(x) for x in others_alignments[i][0:3])
                        break
                if results['gap'] is not None:
                    r1_v2 = results['name']
                    list_reads_v2 = r2_step3[r1_v2]
                    r1_integration_v2 = results['integration'].split(":")
                    r1_integration_chr_v2 = r1_integration_v2[0]
                    r1_integration_pos_v2 = int(r1_integration_v2[1])
                    r1_integration_strand_v2 = r1_integration_v2[2]

                    for read_v2 in list_reads_v2:
                        strand_v2 = "-" if read_v2.is_reverse else "+"
                        if read_v2.reference_name == r1_integration_chr_v2 and (
                                -3000 < read_v2.reference_start - r1_integration_pos_v2 < 3000):
                            target_alignment = [read_v2.reference_name, read_v2.pos, strand_v2, read_v2.cigarstring, read_v2.mapq, 0]
                            results['cigar_r2'] = read_v2.cigarstring
                            results['pos_r2'] = read_v2.reference_start
                            results['strand_r2'] = strand_v2
                            results['target_matches_r2'] = getMatchesFromCigarString(read_v2.cigarstring)
                            break
                    list_results.append(results)

                elif results['n_aav_aln'] > 1:
                    list_only_rear.append(results)
                    r1_onlyaav_counter += 1
                    if len(others_alignments) == 0:
                        r1_noIS_counter += 1
                else:
                    r1_onlyaav_counter += 1
                    if len(others_alignments) == 0:
                        r1_noIS_counter += 1
                    list_no_integration[0].append(read.query_name)
                    list_no_integration[1].append(sa)

            else:
                if read.reference_name=="chrV":
                    r1_onlyaav_counter += 1
                    r1_noIS_counter += 1
                list_no_integration.append(read.query_name)
        except Exception as ex:
            r1_except_counter += 1
            print(read.query_name)
            print(output, "EXCEPTION")
            print(ex)

df_r1 = pd.DataFrame(list_results)
try:
    df_r1['GCperc'] = df_r1.seq_gap.apply(lambda x: (x.count('G') + x.count('C')) / len(x) if len(x) > 3 else None)
except:
    pass
df_r1.to_csv(output + ".R1.tsv", sep="\t", index=False)
print(output, " FINISHED")
