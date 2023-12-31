#!/usr/bin/env python2.7

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--infoFiles", default="${ref_infos}", help="")
parser.add_argument("--datasets", default="${datasets}", help="")
parser.add_argument("--out_prefix", default="${out_prefix}", help="")
parser.add_argument("--infoCutoff", default="${impute_info_cutoff}", help="")

args = parser.parse_args()


def filter_info(infoFiles, datasets, infoCutoff, out_prefix):
    """
    Return:
        well_imputed: certainy >= 1
        SNP_concordance: concord_type0 != -1
    """
    well_imputed = {}
    SNP_concordance = {}
    count = 0
    infoFiles = infoFiles.split(',')
    datasets = datasets.split(',')
    outWell_imputed_out = open(out_prefix + "_well_imputed.tsv", 'w')
    outWell_imputed_snp_out = open(out_prefix + "_well_imputed_snp.tsv", 'w')
    outSNP_accuracy_out = open(out_prefix + "_accuracy.tsv", 'w')
    for infoFile in infoFiles:
        dataset = datasets[infoFiles.index(infoFile)]
        well_imputed[dataset] = []
        SNP_concordance[dataset] = []
        header = []
        # print infoFile
        for line in open(infoFile):
            data = line.strip().split()
            if "SNP" in line and "Rsq" in line:
                if len(header) == 0:
                    header = data
                    info_idx = header.index("Rsq")
                    if 'EmpRsq' in header:
                        conc_idx = header.index("EmpRsq")
                    outWell_imputed_out.writelines('\\t'.join(["GROUPS"] + data) + '\\n')
                    outWell_imputed_snp_out.writelines(data[1] + '\\n')
                    outSNP_accuracy_out.writelines('\\t'.join(["GROUPS"] + data) + '\\n')
            else:
                # print info_idx, data
                if data[info_idx] != '-' and data[info_idx] != '.' and float(data[info_idx]) >= float(infoCutoff):
                    outWell_imputed_out.writelines('\\t'.join([dataset] + data) + '\\n')
                    outWell_imputed_snp_out.writelines(data[1] + '\\n')
                if 'EmpRsq' in header:
                    if data[conc_idx] != '-':
                        outSNP_accuracy_out.writelines('\\t'.join([dataset] + data) + '\\n')
                count += 1
    outWell_imputed_out.close()
    outWell_imputed_snp_out.close()
    outSNP_accuracy_out.close()

if args.infoFiles and args.infoCutoff:
    filter_info(args.infoFiles, args.datasets, args.infoCutoff, args.out_prefix)

