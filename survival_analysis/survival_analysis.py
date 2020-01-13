import pandas as pd
import numpy as np
import seaborn as sns
import gzip
import tarfile
import shutil
import os
import pdb
import glob
from pathlib import Path
from multiprocessing.dummy import Pool as ThreadPool
from sklearn import preprocessing
from lifelines import KaplanMeierFitter
from scipy.stats import ranksums, fisher_exact, ttest_ind
from subprocess import PIPE, run
# Plotting helpers
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 10, 5
sns.set()

verheek = pd.read_table("GSE6891_series_matrix.txt", comment="!")
verheek_surv = pd.read_csv("OS_EFS_GSE6891.csv")
verheek_info = pd.read_table("GSE6891_series_matrix.txt", skiprows=43,
                             nrows=25)

probes = pd.read_table('GPL570-55999.txt', skiprows=16)
probes.index = probes.ID.values

verheek.index = verheek.ID_REF.values
verheek.drop('ID_REF', axis=1, inplace=True)
verheek = verheek.T
verheek_info.index = verheek_info["!Sample_title"]
verheek_info.drop("!Sample_title", axis=1, inplace=True)
verheek_info = verheek_info.T
verheek_info.index = verheek_info.index.str.replace('AML ', "")
verheek_surv.index = verheek_surv.volgnummer.values
verheek_surv.drop('volgnummer', axis=1, inplace=True)
verheek_surv.index = verheek_surv.index.map(str)
verheek_s_i = pd.merge(verheek_info, verheek_surv,
                       left_index=True, right_index=True)
verheek_s_i.index = verheek_s_i["!Sample_geo_accession"].values
verheek_merged = pd.merge(verheek_s_i, verheek,
                          left_index=True, right_index=True)
verheek_merged.dropna(inplace=True)
# fix broken columns and name them better
columns = ['disease state', 'cell type', 'idh1', 'idh2', 'gender',
           'age', 'score', 'risk', 'cebpa', 'karyotype', 'npm1',
           'flt3 itd mutation', 'flt3 tkd mutation', 'n-ras mutation',
           'k-ras mutation', 'evi1 expression', 'cebpa mutation']
verheek_merged.columns.values[8:25] = columns
verheek_merged = verheek_merged[verheek_merged['score'] != 'score: FAB M3']
# make binary column for deaths
lut = dict(zip(verheek_merged.osi.unique(), [0, 1]))
verheek_merged['osi_bin'] = verheek_merged.osi.map(lut)
min_max_scaler = preprocessing.MinMaxScaler()
# probe 206589_at = GFI1
verheek_mm = verheek_merged.copy()
verheek_mm.loc[:, '1007_s_at':'AFFX-TrpnX-M_at'] = min_max_scaler.fit_transform(
    verheek_mm.loc[:, '1007_s_at':'AFFX-TrpnX-M_at'])
verheek_mm_NN = verheek_merged.copy()
verheek_mm_NN = verheek_mm_NN.query("karyotype == 'karyotype: NN'")
verheek_mm_NN.loc[:, '1007_s_at':'AFFX-TrpnX-M_at'] = min_max_scaler.fit_transform(
    verheek_mm_NN.loc[:, '1007_s_at':'AFFX-TrpnX-M_at'])
verheek_mm_noNN = verheek_merged.copy()
verheek_mm_noNN = verheek_mm_noNN.query("karyotype != 'karyotype: NN'")
verheek_mm_noNN.loc[:, '1007_s_at':'AFFX-TrpnX-M_at'] = min_max_scaler.fit_transform(
    verheek_mm_noNN.loc[:, '1007_s_at':'AFFX-TrpnX-M_at'])

## loop to look for most significant genes
x = {}
y = 0
for i in verheek_mm.loc[:, '1007_s_at':'AFFX-TrpnX-M_at'].columns:
    y += 1
    gene_high = verheek_mm.loc[verheek_mm.loc[:, i] > 0.7]
    gene_low = verheek_mm.loc[verheek_mm.loc[:, i] < 0.3]
    if gene_high.shape[0] > 20:
        if gene_low.shape[0] < 400:
            z = ranksums(gene_high.os, gene_low.os)[1]
            if z < 0.05:
                print(i, y)
                x[i] = [gene_high.shape[0], gene_low.shape[0], z]

# change to dataframe
sig = pd.DataFrame.from_dict(x, orient='index')
sig.columns = ['high_no', 'low_no', 'sig']
sig = sig.join(probes['Gene Symbol'])
sig.to_csv('sig.csv')
# resume from here no need to test all sig again
sig = pd.read_csv('sig.csv', index_col=0)

gene = '217975_at' # WBP5
# gene = sig.sort_values('sig').index[4]
symbol = probes.loc[gene]['Gene Symbol']
file_name = symbol+"_"+gene

# gene = sig.sort_values('sig').head(10).index[1]

# function to select high, mid and low expressing patients for later analysis
def high_mid_low(verheek_mm):
    gene_high = verheek_mm.loc[verheek_mm.loc[:, gene] > 0.7]
    gene_mid = verheek_mm.loc[(verheek_mm.loc[:, gene] <= 0.7) & (verheek_mm.loc[:, gene] >= 0.3)]
    gene_low = verheek_mm.loc[verheek_mm.loc[:, gene] < 0.3]
    return(gene_high,gene_mid,gene_low)

gene_high,gene_mid,gene_low = high_mid_low(verheek_mm)
gene_highNN,gene_midNN,gene_lowNN = high_mid_low(verheek_mm_NN)
gene_highnoNN,gene_midnoNN,gene_lownoNN = high_mid_low(verheek_mm_noNN)

# function to randomly sample the patients until they are balanced
# BEWARE! It is not always possible to balance genes due to the 
# nature of the mutations or chromosomal abnormalities
def rand_samp(checks,file_name,
        gene_high,gene_mid,gene_low,
        verheek_merged, verheek_mm):
    table = pd.DataFrame()
    for i in checks:
        x = pd.DataFrame(dict(high=gene_high[i].value_counts(),
                              low=gene_low[i].value_counts()))
        table = table.append(x)
    table.fillna(0, inplace=True)
    p = []
    for i, j in table.values:
        x = [gene_low.shape[0]-j, j], [gene_high.shape[0]-i, i]
        z = fisher_exact(x)[1]
        p.append(z)
    table['p'] = p
    table.to_csv(file_name + "_pre_balanced.csv")
    # random sampling until fisher is not significant or 10000 iterations
    c = 0
    sig = 0
    sig2 = 1
    while c < 10000 and sig2 > 0.05:
        c = c + 1
        print('Iteration number = {}'.format(str(c)))
        table = pd.DataFrame()
        if gene_high.shape[0] > 20:
            gene_high_samp = gene_high.sample(
                np.random.randint(20, gene_high.shape[0]))
        else:
            gene_high_samp = gene_high.copy()
        if gene_mid.shape[0] > 20:
            gene_mid_samp = gene_mid.sample(
                np.random.randint(20, gene_mid.shape[0]))
        else:
            gene_mid_samp = gene_mid.copy()
        if gene_mid.shape[0] > 20:
            gene_low_samp = gene_low.sample(
                np.random.randint(20, gene_low.shape[0]))
        else:
            gene_low_samp = gene_low.copy()
        for i in checks:
            x = pd.DataFrame(dict(high=gene_high_samp[i].value_counts(),
                                  low=gene_low_samp[i].value_counts()))
            table = table.append(x)
        table.fillna(0, inplace=True)
        p = [fisher_exact(([gene_low_samp.shape[0]-j, j],
            [gene_high_samp.shape[0]-i, i]))[1]
             for i, j in table.values]
        # p = []
        # for i,j in table.values:
        #     x = [gene_low.shape[0]-j,j],[gene_high.shape[0]-i,i]
        #     z = fisher_exact(x)[1]
        #     p.append(z)
        table['p'] = p
        p.sort()
        sig = p[0]
        print("Lowest fisher significance = {:.8f}".format(sig))
        if sig > 0.05:
            sig2 = ranksums(gene_high_samp.os, gene_low_samp.os)[1]
            print("Survival significance = {:.3f}".format(sig2))
    a = gene_high_samp.index.values.tolist() \
            + gene_mid_samp.index.values.tolist() + \
        gene_low_samp.index.values.tolist()
    print("exporting balanced unscaled to csv...")
    verheek_merged.loc[a].to_csv(file_name+"_balanced_unscaled.csv")
    print("gzipping balanced unscaled csv...")
    with open(file_name+"_balanced_unscaled.csv", 'rb') as f_in:
        with gzip.open(file_name+"_balanced_unscaled.csv.gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file_name+"_balanced_unscaled.csv")
    print("exporting balanced scaled to csv...")
    verheek_mm.loc[a].to_csv(file_name+"_balanced.csv")
    print("gzipping balanced scaled csv...")
    with open(file_name+"_balanced.csv", 'rb') as f_in:
        with gzip.open(file_name+"_balanced.csv.gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file_name+"_balanced.csv")

checks = ['idh1', 'idh2', 'gender', 'score', 'cebpa', 'karyotype',
          'npm1', 'flt3 itd mutation', 'flt3 tkd mutation',
          'n-ras mutation', 'k-ras mutation', 'evi1 expression', ]



checks2 = ['idh1', 'idh2', 'gender', 'score', 'cebpa',
          'npm1', 'flt3 itd mutation', 'flt3 tkd mutation',
          'n-ras mutation', 'k-ras mutation', 'evi1 expression', ]
file_name_NN = file_name + "_NN"
file_name_noNN = file_name + "noNN"

rand_samp(checks,file_name,gene_high,gene_mid,gene_low,
        verheek_merged, verheek_mm)
rand_samp(checks2,file_name_NN,gene_highNN,gene_midNN,gene_lowNN,
        verheek_merged, verheek_mm)
rand_samp(checks2,file_name_noNN,gene_highnoNN,gene_midnoNN,gene_lownoNN ,
        verheek_merged, verheek_mm)
# resume from here no need to run all again 

# random sampling can take a long time, this function allows resuming after
# the datasets are balanced
def resume(file_name):
    my_file = Path(file_name+"tar.gz")
    if my_file.is_file():
        tar = tarfile.open(file_name+".tar.gz")
        bal = tar.extractfile(file_name+"_balanced.csv.gz")
        with gzip.open(bal, 'rb') as f:
            balanced = pd.read_csv(f, index_col=0)
        tar.close()
    else:
        with gzip.open(file_name+"_balanced.csv.gz", 'rb') as f:
            balanced = pd.read_csv(f, index_col=0)
    return(balanced)

balanced = resume(file_name)
balancedNN = resume(file_name_NN)
balancednoNN = resume(file_name_noNN)


def survival(df, probe, surv_type):
    '''function to do stats and plot survival for a specific probe
    from a dataframe. type = os or efs'''
    gene_high = df.loc[df.loc[:, probe] > 0.7]
    gene_mid = df.loc[(df.loc[:, probe] <= 0.7) & (df.loc[:, probe] >= 0.3)]
    gene_low = df.loc[df.loc[:, probe] < 0.3]
    kmf = KaplanMeierFitter()
    ax = plt.subplot(111)
    if surv_type == 'os':
        kmf.fit(gene_low.os, gene_low.osi_bin, label=['low'])
        kmf.survival_function_.plot(ax=ax)
        #kmf.fit(gene_mid.os, gene_mid.osi_bin, label=['mid'])
        #kmf.survival_function_.plot(ax=ax)
        kmf.fit(gene_high.os, gene_high.osi_bin, label=['high'])
        kmf.survival_function_.plot(ax=ax)
    if surv_type == 'ef':
        kmf.fit(gene_low.efs, gene_low.osi_bin, label=['low'])
        kmf.survival_function_.plot(ax=ax)
        #kmf.fit(gene_mid.efs, gene_mid.osi_bin, label=['mid'])
        #kmf.survival_function_.plot(ax=ax)
        kmf.fit(gene_high.efs, gene_high.osi_bin, label=['high'])
        kmf.survival_function_.plot(ax=ax)
    plt.title(file_name+"_"+surv_type)
    plt.ylabel(surv_type+' (days)')
    plt.ylim(ymin=0)
    plt.gcf().subplots_adjust(bottom=0.25)
    if surv_type == 'os':
        plt.figtext(0.5,0,
        'high (70-100%) expression n = {}, mean = {:.2f}\n\
        low (0-30%) expression = {}, mean = {:.2f}\n\
        p-value is {:.5f}'.format(gene_high.shape[0],gene_high.os.mean(),
                                gene_low.shape[0],gene_low.os.mean(),
                            ranksums(gene_high.os, gene_low.os)[1]),
                ha='center',
                wrap=True)
    if surv_type == 'ef':
        plt.figtext(0.5,0,
        'high (70-100%) expression n = {}, mean = {:.2f}\n\
        low (0-30%) expression = {}, mean = {:.2f}\n\
        p-value is {:.5f}'.format(gene_high.shape[0],gene_high.efs.mean(),
                                gene_low.shape[0],gene_low.efs.mean(),
                            ranksums(gene_high.efs, gene_low.efs)[1]),
                ha='center',
                wrap=True)
    kmf2 = plt.gcf()
    # test significance
    return(kmf2)
def surv_plots(full,balanced,file_name):
    x = survival(verheek_mm, gene, 'os')
    x.savefig(file_name+"_pre-balance_overall_survival.svg", dpi=300)
    plt.close()
    x = survival(balanced, gene, 'os')
    x.savefig(file_name+"_post-balance_overall_survival.svg", dpi=300)
    plt.close()
    x = survival(verheek_mm, gene, 'ef')
    x.savefig(file_name+"_pre-balance_ef_survival.svg", dpi=300)
    plt.close()
    x = survival(balanced, gene, 'ef')
    x.savefig(file_name+"_post-balance_ef_survival.svg", dpi=300)
    plt.close()

surv_plots(verheek_mm,balanced,file_name)
surv_plots(verheek_mm_NN,balancedNN,file_name_NN)
surv_plots(verheek_mm_noNN,balancednoNN,file_name_noNN)

# extract patient information for downstream analysis in R
def high_low_to_R(balanced, file_name, GSE):
    gene_high_bal = balanced.loc[balanced.loc[:, gene] > 0.7]
    gene_low_bal = balanced.loc[balanced.loc[:, gene] < 0.3]
    high_and_low = gene_high_bal.index.values.tolist() + \
        gene_low_bal.index.values.tolist()
    if GSE == "GSE6891":
        mid = np.logical_not(verheek.index.isin(high_and_low))
        mid = verheek.index[mid]
    if GSE == "GSE15434":
        mid = np.logical_not(kohl.index.isin(high_and_low))
        mid = kohl.index[mid]
    if GSE == "GSE1159":
        mid = np.logical_not(valk.index.isin(high_and_low))
        mid = valk.index[mid]
    a = '0'
    b = []
    for i in range(len(gene_high.index)):
        b.append(a)
    a = '1'
    for i in range(len(gene_low.index)):
        b.append(a)
    a = 'X'
    for i in range(len(mid)):
        b.append(a)
    a = gene_high.index.values.tolist() + gene_low.index.values.tolist() + \
        mid.values.tolist()
    mapping = dict(zip(a, b))
    if GSE == "GSE6891":
        mapped = verheek.index.map(mapping)
    if GSE == "GSE15434":
        mapped = kohl.index.map(mapping)
    if GSE == "GSE1159":
        mapped = valk.index.map(mapping)
    str1 = ''.join(str(e) for e in mapped)
    with open(file_name+"high_low_"+GSE+"_mapping.txt", "w") as f:
        f.write(str1)
    highDf = pd.DataFrame(gene_high_bal.index.tolist())
    highDf['highvslow'] = 'high'
    lowDf = pd.DataFrame(gene_low_bal.index.tolist())
    lowDf['highvslow'] = 'low'
    highvslowDf = pd.concat([highDf,lowDf])
    highvslowDf.columns = ["patient", "highvlow"]
    highvslowDf.to_csv(file_name+"_"+GSE+"_highvslow.csv", index=False)

high_low_to_R(balanced, file_name, "GSE6891")
high_low_to_R(balancedNN,file_name_NN,"GSE6891")
high_low_to_R(balancednoNN,file_name_noNN,"GSE6891")

# check gene in other datasets then check balancing
kohl = pd.read_table("GSE15434_series_matrix.txt", comment="!")
kohl_info = pd.read_table("GSE15434_series_matrix.txt", 
        skiprows=39, nrows=40)
kohl.index = kohl.ID_REF.values
kohl.drop('ID_REF', axis= 1, inplace=True)
kohl = kohl.T
kohl_info.index = kohl_info["!Sample_geo_accession"].values
kohl_info.drop("!Sample_geo_accession", axis=1, inplace=True)
kohl_info = kohl_info.T
kohl_m = pd.merge(kohl_info, kohl, left_index=True, right_index=True)
columns = kohl_m.iloc[1,7:16].str.replace(":.*","").values.tolist()
kohl_m.columns.values[7:16] = columns
kohl_mm = kohl_m.copy()
kohl_mm.loc[:, '1007_s_at':'AFFX-TrpnX-M_at'] = min_max_scaler.fit_transform(
    kohl_mm.loc[:, '1007_s_at':'AFFX-TrpnX-M_at'])
gene_high = kohl_mm.loc[kohl_mm.loc[:, gene] > 0.7]
gene_mid = kohl_mm.loc[(kohl_mm.loc[:, gene] <= 0.7) & (kohl_mm.loc[:, gene] >= 0.3)]
gene_low = kohl_mm.loc[kohl_mm.loc[:, gene] < 0.3]
checks = ["npm1", "flt3-itd", "cebpa"]
table = pd.DataFrame()
for i in checks:
    x = pd.DataFrame(dict(high=gene_high[i].value_counts(),
                          low=gene_low[i].value_counts()))
    table = table.append(x)
table.fillna(0, inplace=True)
p = []
for i, j in table.values:
    x = [gene_low.shape[0]-j, j], [gene_high.shape[0]-i, i]
    z = fisher_exact(x)[1]
    p.append(z)
table['p'] = p
table.to_csv(file_name + "_kohl_pre_balanced.csv")

# random sampling until fisher is not significant or 10000 iterations

c = 0
sig = 0
sig2 = 1
while c < 10000 and sig2 > 0.05:
    c = c + 1
    print('Iteration number = {}'.format(str(c)))
    table = pd.DataFrame()
    if gene_high.shape[0] > 20:
        gene_high_samp = gene_high.sample(
            np.random.randint(20, gene_high.shape[0]))
    else:
        gene_high_samp = gene_high.copy()
    if gene_mid.shape[0] > 20:
        gene_mid_samp = gene_mid.sample(
            np.random.randint(20, gene_mid.shape[0]))
    else:
        gene_mid_samp = gene_mid.copy()
    if gene_mid.shape[0] > 20:
        gene_low_samp = gene_low.sample(
            np.random.randint(20, gene_low.shape[0]))
    else:
        gene_low_samp = gene_low.copy()
    for i in checks:
        x = pd.DataFrame(dict(high=gene_high_samp[i].value_counts(),
                              low=gene_low_samp[i].value_counts()))
        table = table.append(x)
    table.fillna(0, inplace=True)
    p = [fisher_exact(([gene_low_samp.shape[0]-j, j], [gene_high_samp.shape[0]-i, i]))[1]
         for i, j in table.values]
    # p = []
    # for i,j in table.values:
    #     x = [gene_low.shape[0]-j,j],[gene_high.shape[0]-i,i]
    #     z = fisher_exact(x)[1]
    #     p.append(z)
    table['p'] = p
    p.sort()
    sig = p[0]
    print("Lowest fisher significance = {:.8f}".format(sig))
    if sig > 0.05:
        break
a = gene_high_samp.index.values.tolist() + \
        gene_mid_samp.index.values.tolist() + \
    gene_low_samp.index.values.tolist()
print("exporting balanced unscaled to csv...")

kohl_m.loc[a].to_csv(file_name+"_balanced_kohl_unscaled.csv")
kohl_bal_unscaled = kohl.loc[a]
kohl_bal_scaled = kohl_mm.loc[a]
high_low_to_R(kohl_bal_scaled,file_name,"GSE15434")

list_to_check = ['HOXA','HOXB1','MEIS','PBX']
probe_gene = list()
for i in list_to_check:
    probe_gene.append(probes.loc[probes['Gene Symbol'].str.contains(i, na=False)].ID.values.tolist())
flat_list = [item for sublist in probe_gene for item in sublist]
flat_list.append(gene)
columns_to_keep = kohl_m.loc[a].columns[:40].tolist()
columns_to_keep = columns_to_keep +  flat_list
kohl_m.loc[a,columns_to_keep].to_csv(file_name + "_balanced_kohl_selected_genes.csv")
hafer = pd.read_csv("Hafer_WBP5_AML.csv")
hafer.index = hafer.Probeset_id.values
genes = hafer.Gene_id
hafer.drop(["Probeset_id","Gene_id"], axis=1,inplace=True)
hafer = hafer.T
hafer = hafer.rename(columns={"Probeset_id":"Status"})
hafer = hafer.loc[hafer["Status"].isin(['T(8;21)', 'INV(16)', 'MLL', 'NN', 'COMPLEX'])]
hafer_mm = hafer.copy()
hafer_mm.loc[:, '217975_at':'217975_at'] = min_max_scaler.fit_transform(
    hafer_mm.loc[:, '217975_at':'217975_at'])
lb = preprocessing.LabelBinarizer()
status = pd.DataFrame(lb.fit_transform(hafer.Status))
status.index = hafer.index
status.columns = lb.classes_.tolist() 
hafer_mm = pd.merge(status, hafer_mm, left_index=True, right_index=True)
hafer_mm.drop("Status", axis=1, inplace=True)
gene_high = hafer_mm.loc[hafer_mm.loc[:, gene] > 0.7]
gene_mid = hafer_mm.loc[(hafer_mm.loc[:, gene] <= 0.7) & (hafer_mm.loc[:, gene] >= 0.3)]
gene_low = hafer_mm.loc[hafer_mm.loc[:, gene] < 0.3]
checks = lb.classes_.tolist()
table = pd.DataFrame()
for i in checks:
    x = pd.DataFrame(dict(high=gene_high[i].value_counts(),
                          low=gene_low[i].value_counts()))
    table = table.append(x)
table.fillna(0, inplace=True)
p = []
for i, j in table.values:
    x = [gene_low.shape[0]-j, j], [gene_high.shape[0]-i, i]
    z = fisher_exact(x)[1]
    p.append(z)
table['p'] = p
table.to_csv(file_name + "_hafer_pre_balanced.csv")
# random sampling until fisher is not significant or 10000 iterations
c = 0
sig = 0
sig2 = 1
while c < 10000 and sig2 > 0.05:
    c = c + 1
    print('Iteration number = {}'.format(str(c)))
    table = pd.DataFrame()
    if gene_high.shape[0] > 20:
        gene_high_samp = gene_high.sample(
            np.random.randint(20, gene_high.shape[0]))
    else:
        gene_high_samp = gene_high.copy()
    if gene_mid.shape[0] > 20:
        gene_mid_samp = gene_mid.sample(
            np.random.randint(20, gene_mid.shape[0]))
    else:
        gene_mid_samp = gene_mid.copy()
    if gene_mid.shape[0] > 20:
        gene_low_samp = gene_low.sample(
            np.random.randint(20, gene_low.shape[0]))
    else:
        gene_low_samp = gene_low.copy()
    for i in checks:
        x = pd.DataFrame(dict(high=gene_high_samp[i].value_counts(),
                              low=gene_low_samp[i].value_counts()))
        table = table.append(x)
    table.fillna(0, inplace=True)
    p = [fisher_exact(([gene_low_samp.shape[0]-j, j], [gene_high_samp.shape[0]-i, i]))[1]
         for i, j in table.values]
    # p = []
    # for i,j in table.values:
    #     x = [gene_low.shape[0]-j,j],[gene_high.shape[0]-i,i]
    #     z = fisher_exact(x)[1]
    #     p.append(z)
    table['p'] = p
    p.sort()
    sig = p[0]
    print("Lowest fisher significance = {:.8f}".format(sig))
    if sig > 0.05:
        break
a = gene_high_samp.index.values.tolist() + gene_mid_samp.index.values.tolist() + \
    gene_low_samp.index.values.tolist()

print("exporting balanced unscaled to csv...")
hafer.loc[a].to_csv(file_name+"_balanced_hafer_unscaled.csv")
list_to_check = ['HOXA','HOXB1','MEIS','PBX']
probe_gene = list()
for i in list_to_check:
    probe_gene.append(probes.loc[probes['Gene Symbol'].str.contains(i, na=False)].ID.values.tolist())
flat_list = [item for sublist in probe_gene for item in sublist]
flat_list.append(gene)
columns_to_keep = hafer.loc[a].columns[:40].tolist()
columns_to_keep = columns_to_keep +  flat_list
hafer.loc[a,columns_to_keep].to_csv(file_name + "_balanced_hafer_selected_genes.csv")

# check gene in other datasets then check balancing
valk = pd.read_table("GSE1159_series_matrix.txt", comment="!")
valk_info = pd.read_table("GSE1159_series_matrix.txt", skiprows=41, nrows=24)
valk_more_info = pd.read_csv("GSE1159_moreinfo.csv")
valk_more_info.index = valk_more_info['Array Samples'].values
valk_more_info.drop('Array Samples', axis=1, inplace=True)
valk.index = valk.ID_REF.values
valk.drop('ID_REF', axis= 1, inplace=True)
valk = valk.T
valk_info.index = valk_info["!Sample_geo_accession"].values
valk_info.drop("!Sample_geo_accession", axis=1, inplace=True)
valk_info = valk_info.T
valk_m = pd.merge(valk_info, valk, left_index=True, right_index=True)
columns = valk_m.iloc[1,7:16].str.replace(":.*","").values.tolist()
valk_m.columns.values[7:16] = columns
valk_m = pd.merge(valk_more_info, valk_m, left_index=True,right_index=True)
valk_m = valk_m[valk_m['FAB'] != 'M3']
valk_mm = valk_m.copy()
valk_mm.loc[:, '1007_s_at':'AFFX-TrpnX-M_at'] = min_max_scaler.fit_transform(
    valk_mm.loc[:, '1007_s_at':'AFFX-TrpnX-M_at'])
gene_high = valk_mm.loc[valk_mm.loc[:, gene] > 0.7]
gene_mid = valk_mm.loc[(valk_mm.loc[:, gene] <= 0.7) & (valk_mm.loc[:, gene] >= 0.3)]
gene_low = valk_mm.loc[valk_mm.loc[:, gene] < 0.3]
checks = ["FAB","Cytogenetics", "FLT3 ITD", "FLT3 TKD", "NRAS","KRAS","EVI1",
        "CEBPA"]
table = pd.DataFrame()
for i in checks:
    x = pd.DataFrame(dict(high=gene_high[i].value_counts(),
                          low=gene_low[i].value_counts()))
    table = table.append(x)
table.fillna(0, inplace=True)
p = []
for i, j in table.values:
    x = [gene_low.shape[0]-j, j], [gene_high.shape[0]-i, i]
    z = fisher_exact(x)[1]
    p.append(z)
table['p'] = p
table.to_csv(file_name + "_valk_pre_balanced.csv")
# random sampling until fisher is not significant or 10000 iterations
c = 0
sig = 0
sig2 = 1
while c < 10000 and sig2 > 0.05:
    c = c + 1
    print('Iteration number = {}'.format(str(c)))
    table = pd.DataFrame()
    if gene_high.shape[0] > 20:
        gene_high_samp = gene_high.sample(
            np.random.randint(20, gene_high.shape[0]))
    else:
        gene_high_samp = gene_high.copy()
    if gene_mid.shape[0] > 20:
        gene_mid_samp = gene_mid.sample(
            np.random.randint(20, gene_mid.shape[0]))
    else:
        gene_mid_samp = gene_mid.copy()
    if gene_mid.shape[0] > 20:
        gene_low_samp = gene_low.sample(
            np.random.randint(20, gene_low.shape[0]))
    else:
        gene_low_samp = gene_low.copy()
    for i in checks:
        x = pd.DataFrame(dict(high=gene_high_samp[i].value_counts(),
                              low=gene_low_samp[i].value_counts()))
        table = table.append(x)
    table.fillna(0, inplace=True)
    p = [fisher_exact(([gene_low_samp.shape[0]-j, j],
        [gene_high_samp.shape[0]-i, i]))[1]
         for i, j in table.values]
    # p = []
    # for i,j in table.values:
    #     x = [gene_low.shape[0]-j,j],[gene_high.shape[0]-i,i]
    #     z = fisher_exact(x)[1]
    #     p.append(z)
    table['p'] = p
    p.sort()
    sig = p[0]
    print("Lowest fisher significance = {:.8f}".format(sig))
    if sig > 0.05:
        break
a = gene_high_samp.index.values.tolist() + \
        gene_mid_samp.index.values.tolist() + \
    gene_low_samp.index.values.tolist()
print("exporting balanced unscaled to csv...")
valk_m.loc[a].to_csv(file_name+"_balanced_valk_unscaled.csv")
valk_bal_unscaled = valk.loc[a]
valk_bal_scaled = valk_mm.loc[a]

high_low_to_R(valk_bal_scaled,file_name,"GSE1159")

file_name = symbol+"_"+gene+"_gse37642_"

gse37642 = pd.read_table("GSE37642-GPL96_series_matrix.txt.gz",
        skiprows=80, sep="\t")
gse37642_info = pd.read_table("GSE37642-GPL96_series_matrix.txt.gz",
        skiprows=42,nrows=25)
gse37642_surv = pd.read_table("GSE37642_Survival_data.txt.gz",
        header=None, sep="\t", skiprows=1)
gse37642_surv.columns = ["GSM","patient","os","life_status"]
gse37642.index = gse37642.ID_REF.values
gse37642.drop('ID_REF', axis=1, inplace=True)
gse37642 = gse37642.T
index_col = "!Sample_title"
gse37642_info.index = gse37642_info[index_col].values
gse37642_info.drop(index_col, axis=1, inplace=True)
gse37642_info = gse37642_info.T
id_column = "!Sample_geo_accession"
gse37642_info.index = gse37642_info[id_column].values
gse37642_info.drop(id_column, axis=1, inplace=True)
gse37642_surv.index = gse37642_surv.GSM.values
gse37642_surv.drop('GSM', axis=1, inplace=True)
gse37642_surv.index = gse37642_surv.index.map(str)
gse37642_s_i = pd.merge(gse37642_info, gse37642_surv,
                       left_index=True, right_index=True)
gse37642_merged = pd.merge(gse37642_s_i, gse37642,
                          left_index=True, right_index=True)
gse37642_merged.columns.values[8:9] = "fab"
gse37642_merged.columns.values[10:12] = ["runx1_fusion","runx1_mutation"]
gse37642_merged = gse37642_merged[gse37642_merged['fab'] != 'fab: 3']
gse37642_merged = gse37642_merged[gse37642_merged['fab'] != 'fab: 3v']
gse37642_merged = gse37642_merged[gse37642_merged['life_status'].isin(['alive','dead'])]
# make binary column for deaths
lut = dict(zip(gse37642_merged.life_status.unique(), [1, 0]))
gse37642_merged['osi_bin'] = gse37642_merged.life_status.map(lut)
min_max_scaler = preprocessing.MinMaxScaler()
gse37642_mm = gse37642_merged.copy()
gse37642_mm.loc[:, '1007_s_at':'AFFX-TrpnX-M_at'] = min_max_scaler.fit_transform(
    gse37642_mm.loc[:, '1007_s_at':'AFFX-TrpnX-M_at'])
x = survival(gse37642_mm, gene, 'os')
x.savefig(file_name+"_pre-balance_overall_survival.svg", dpi=300)
plt.close()

gene_high,gene_mid,gene_low = high_mid_low(gse37642_mm)
checks = ["fab","runx1_fusion","runx1_mutation"]

rand_samp(checks,file_name,gene_high,gene_mid,gene_low,
        gse37642_merged, gse37642_mm)

balanced = resume(file_name)
x = survival(balanced, gene, 'os')

x.savefig(file_name+"_post-balance_overall_survival.svg", dpi=300)
plt.close()

# GSE12417
file_name = symbol+"_"+gene+"_gse12417_"
gse12417 = pd.read_table("GSE12417-GPL570_series_matrix.txt.gz",
        skiprows=71, sep="\t")
gse12417_info = pd.read_table("GSE12417-GPL570_series_matrix.txt.gz",
        skiprows=36,nrows=25)

def merge_dfs(x, x_info):
    x.index = x.ID_REF.values
    x.drop('ID_REF', axis=1, inplace=True)
    x = x.T
    index_col = "!Sample_title"
    x_info.index = x_info[index_col].values
    x_info.drop(index_col, axis=1, inplace=True)
    x_info = x_info.T
    id_column = "!Sample_geo_accession"
    x_info.index = x_info[id_column].values
    x_info.drop(id_column, axis=1, inplace=True)
    x_merged = pd.merge(x_info, x,
                              left_index=True, right_index=True)
    return(x_merged)

gse12417_a = merge_dfs(gse12417, gse12417_info)
table = "GSE12417-GPL96_series_matrix.txt.gz"
gse12417 = pd.read_table(table,skiprows=71, sep="\t")
gse12417_info = pd.read_table(table,skiprows=36,nrows=25)
gse12417_b = merge_dfs(gse12417, gse12417_info)
# WBP5 not in this array
#table = "GSE12417/GSE12417-GPL97_series_matrix.txt.gz"
#gse12417 = pd.read_table(table,skiprows=70, sep="\t")
#gse12417_info = pd.read_table(table,skiprows=36,nrows=25)
#gse12417_c = merge_dfs(gse12417, gse12417_info)
gse12417_a = gse12417_a[['!Sample_characteristics_ch1',gene]]
gse12417_b = gse12417_b[['!Sample_characteristics_ch1',gene]]
gse12417 = pd.concat([gse12417_a, gse12417_b])
new = gse12417["!Sample_characteristics_ch1"].str.split(";", expand=True)
gse12417['disease'] = new[0]
gse12417['os'] = new[2]
gse12417['osi'] = new[3]
gse12417.drop(columns=["!Sample_characteristics_ch1"], inplace=True)
gse12417['os'] = gse12417['os'].str.replace('OS = ','')
gse12417['os'] = gse12417['os'].str.replace(' days','')
gse12417['os'] = gse12417['os'].astype(int)
gse12417['osi'] = gse12417['osi'].str.replace('0','alive') 
gse12417['osi'] = gse12417['osi'].str.replace('1','dead') 
lut = dict(zip(gse12417.osi.unique(), [1, 0]))
gse12417['osi_bin'] = gse12417.osi.map(lut)
gse12417['disease'] = gse12417['disease'].str.replace("\(.*\)", "") 
gse12417['disease'] = gse12417['disease'].str.replace(",", "") 
gse12417 = gse12417.loc[gse12417['disease'] != 'MDS normal karyotype  MDS RAEB']
gse12417 = gse12417.loc[gse12417['disease'] != 'MDS RAEB normal karyotype  MDS RAEB']
min_max_scaler = preprocessing.MinMaxScaler()
gse12417_mm = gse12417.copy()
gse12417_mm.loc[:, '217975_at':'217975_at'] = min_max_scaler.fit_transform(
        gse12417_mm.loc[:, '217975_at':'217975_at'])
x = survival(gse12417_mm, gene, 'os')
x.savefig(file_name+"_pre-balance_overall_survival.svg", dpi=300)
plt.close()

gene_high,gene_mid,gene_low = high_mid_low(gse12417_mm)
checks = ["disease"]
rand_samp(checks,file_name,gene_high,gene_mid,gene_low,
        gse12417, gse12417_mm)
balanced = resume(file_name)
x = survival(balanced, gene, 'os')
x.savefig(file_name+"_post-balance_overall_survival.svg", dpi=300)
plt.close()


