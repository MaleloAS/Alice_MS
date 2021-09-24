#!/usr/bin/env python3

import pysam
import argparse

parser = argparse.ArgumentParser(description= 'Add tags to VCF')
parser.add_argument('invcf', default= '-', nargs= '?', help= 'Input VCF [%(default)s]')
args = parser.parse_args()

invcf = pysam.VariantFile(args.invcf, "r")

invcf.header.formats.add("ALT_AF", "A", "Float", "Alternate allele frequencies")
invcf.header.formats.add("SUM_ALT", "1", "Integer", "Sum of observations from all alternate alleles")
invcf.header.formats.add("SUM_ALT_AF", "1", "Float", "Sum of alternate allele frequencies")

outvcf = pysam.VariantFile('-', 'w', header= invcf.header)

for variant in invcf:
    for sample in variant.samples:
        ad = variant.samples[sample]['AD']
        dp = float(sum(ad))
        alt_af = [x/dp for x in ad[1:]]
        sumalt = sum(ad) - ad[0]
        variant.samples[sample]['ALT_AF'] = alt_af
        variant.samples[sample]['SUM_ALT'] = sumalt
        variant.samples[sample]['SUM_ALT_AF'] = sum(alt_af)
    outvcf.write(variant)
invcf.close()
