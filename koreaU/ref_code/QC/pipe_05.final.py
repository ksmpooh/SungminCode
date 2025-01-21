import os
import time
import hail as hl
from pprint import pprint

from hail.plot import show
from gnomad_qc.v2.resources.sample_qc import *
from gnomad_qc.v2.resources import get_gnomad_data

import argparse
import sys
import logging
import csv
import pandas as pd
import numpy as np
import glob
import time
from tqdm import tqdm


from gnomad.utils.filtering import filter_to_autosomes
from pprint import pprint

from bokeh.io import output_notebook, show
from bokeh.layouts import gridplot
from bokeh.models import Span
from bokeh.transform import factor_cmap, factor_mark
from bokeh.plotting import figure, output_file, save

os.chdir("/data1/mycho/WGS_AD_2011/")

vcf_dir = "/data2/Oneomics/VCF_hg38.1-6/batch1-7/"
fn_vcf = "SMC_WGS.batch1-7.vcf.gz"

#dir_output = "/data2/Oneomics/VCF_hg38.1-6/batch1-7/data"
dir_output = "/data1/mycho/WGS_AD_2011/4.VCF.QCed.set1-7/"
prefix = "SMC_batch1-7_hg38"
fn_LCR="/data1/mycho/Source/annotation/LCR-hs38.auto_XY.bed"
fn_meta="/data1/mycho/WGS_AD_2011/1.data/phenotype/WGS_1824.pheno_update.230726"
fn_sample_list="/data1/mycho/WGS_AD_2011/4.VCF.QCed.set1-7/data/sample.list"

config = { 'spark.driver.memory':'200g', 'spark.local.dir' : '/data1/mycho/WGS_AD_2011/temp'}
hl.init(master='local[40]', spark_conf=config,  default_reference='GRCh38')
###01.1 load matrix file
mt = hl.read_matrix_table(dir_output +'/data/' + prefix +'.variantQC2.final.mt')
print(dir_output +'/data/' + prefix +'.variantQC2.final.mt')

samples_to_remove = {"WGS_0349","WGS_0507","WGS_0522","WGS_0527","WGS_0528","WGS_0465","WGS_0505","WGS_0506","WGS_0509","WGS_0510","WGS_0511","WGS_0524","WGS_0526","WGS_0532","WGS_0533","WGS_0539","WGS_0547","WGS_0548","WGS_0568","WGS_0590","WGS_0595"}

set_to_remove = hl.literal(samples_to_remove)
mt = mt.filter_cols(~set_to_remove.contains(mt['s']))
print('After Sample QC, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))

#hail variantQC & sampleQC
mt=hl.variant_qc(mt)
mt=hl.sample_qc(mt)

##not working
#old_sample_names = mt.s.collect()
#mapping = {"WGS_0285":"WGS_0930", "WGS_0693":"WGS_0929", "WGS_1639":"WGS_1637"}
#mt = mt.annotate_cols(new_s = mapping.get(mt.s, mt.s))
#mt = mt.key_cols_by(s=mt.new_s)
#mt = mt.drop('new_s')

#https://github.com/dnanexus/OpenBio/blob/master/hail_tutorial/replace_id.ipynb
mapping_table = hl.import_table('/data1/mycho/WGS_AD_2011/1.data/phenotype/WGS_1824.pheno_update.230726_IDconversion.txt', key='old_sample_id')
mt = mt.annotate_cols(**mapping_table[mt.s])
mt = mt.key_cols_by(s = mt.new_sample_id).drop("new_sample_id")

print('After Sample QC, %d samples and %d variants remain.' % (mt.count_cols(), mt.count_rows()))
mt.write(dir_output +'/data/' + prefix +'.final.230808.mt', overwrite=True)
hl.export_vcf(mt, dir_output +'/data/' + prefix +'.QCed.final_230808.vcf.bgz')
hl.export_plink(mt, dir_output +'/data/' + prefix +'.QCed.final_230808', ind_id = mt.s, fam_id= mt.s)


#mt = hl.read_matrix_table(dir_output +'/data/' + prefix +'.final.230808.mt')

#for chromosome in range(1,23):
#    chromosome_mt = mt.filter_rows(mt.locus.contig == "chr" + str(chromosome))
#    chromosome_mt.write(dir_output +'/data/' + prefix +'.chr'+str(chromosome) +'.final.230808.mt', overwrite=True)
#    hl.export_vcf(chromosome_mt, dir_output +'/data/' + prefix + '.chr'+str(chromosome) +'.QCed.final_230808.vcf.bgz')


