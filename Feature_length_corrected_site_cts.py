import pandas as pd
site_tbl=pd.read_csv("Annotated_Results_CpG_intersection.csv",sep='\t',low_memory=False)
out_f=open("CpG_sitefound_stats_072620.tsv",'w')
out_f.write("Annotation\tCount\n")
for anno,counts in site_tbl.groupby("Annotation"):
    out_f.write(anno+"\t"+str(counts['start'].count())+"\n")
out_f.close()

site_tbl=pd.read_csv("Annotated_Results_CHG_intersection.csv",sep='\t',low_memory=False)
out_f=open("CHG_sitefound_stats_072620.tsv",'w')
out_f.write("Annotation\tCount\n")
for anno,counts in site_tbl.groupby("Annotation"):
    out_f.write(anno+"\t"+str(counts['start'].count())+"\n")
out_f.close()

site_tbl=pd.read_csv("Annotated_Results_CHH_intersection.csv",sep='\t',low_memory=False)
out_f=open("CHH_sitefound_stats_072620.tsv",'w')
out_f.write("Annotation\tCount\n")
for anno,counts in site_tbl.groupby("Annotation"):
    out_f.write(anno+"\t"+str(counts['start'].count())+"\n")
out_f.close()
del site_tbl
called_cpg=pd.read_csv("CpG_sitefound_stats_072620.tsv",sep="\t").rename(columns={"Count":"CpGs"}).fillna(0.0)
called_chg=pd.read_csv("CHG_sitefound_stats_072620.tsv",sep="\t").rename(columns={"Count":"CHGs"}).fillna(0.0)
called_chh=pd.read_csv("CHH_sitefound_stats_072620.tsv",sep="\t").rename(columns={"Count":"CHHs"}).fillna(0.0)
len_df=pd.read_csv("Pico_feature_lengths_2.tsv",sep="\t",low_memory=False)
len_df['Annotation']=len_df["Annotation"].replace(",","|",regex=True)
norm_counts=pd.merge(len_df,called_cpg,on="Annotation",how='left')
norm_counts=pd.merge(norm_counts,called_chg,on="Annotation",how='left')
norm_counts=pd.merge(norm_counts,called_chh,on="Annotation",how='left')
norm_counts['CPGpkb']=norm_counts['CpGs']*1000/norm_counts["Length"]
norm_counts['CHGpkb']=norm_counts['CHGs']*1000/norm_counts["Length"]
norm_counts['CHHpkb']=norm_counts['CHHs']*1000/norm_counts["Length"]
norm_counts=norm_counts.fillna(0.0)
norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)].to_csv("Norm_counts_called_Genebody_072620.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==True)].to_csv("Norm_counts_called_Promoter_072620.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts["Annotation"].str.contains("_T")==True)].to_csv("Norm_counts_called_Terminator_072620.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts["Annotation"].str.contains("_IGR_")==True)].to_csv("Norm_counts_called_Intergenic_072620.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts['CPGpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)].to_csv("Norm_counts_called_Genebody_072620_CpGnonzero.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts['CPGpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==True)].to_csv("Norm_counts_called_Promoter_072620_CpGnonzero.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts['CPGpkb']>0) & (norm_counts["Annotation"].str.contains("_T")==True)].to_csv("Norm_counts_called_Terminator_072620_CpGnonzero.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts['CPGpkb']>0) & (norm_counts["Annotation"].str.contains("_IGR_")==True)].to_csv("Norm_counts_called_Intergenic_072620_CpGnonzero.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts['CHGpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)].to_csv("Norm_counts_called_Genebody_072620_CHGnonzero.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts['CHGpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==True)].to_csv("Norm_counts_called_Promoter_072620_CHGnonzero.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts['CHGpkb']>0) & (norm_counts["Annotation"].str.contains("_T")==True)].to_csv("Norm_counts_called_Terminator_072620_CHGnonzero.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts['CHGpkb']>0) & (norm_counts["Annotation"].str.contains("_IGR_")==True)].to_csv("Norm_counts_called_Intergenic_072620_CHGnonzero.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts['CHHpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)].to_csv("Norm_counts_called_Genebody_072620_CHHnonzero.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts['CHHpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==True)].to_csv("Norm_counts_called_Promoter_072620_CHHnonzero.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts['CHHpkb']>0) & (norm_counts["Annotation"].str.contains("_T")==True)].to_csv("Norm_counts_called_Terminator_072620_CHHnonzero.tsv",sep='\t',index=False)
norm_counts.loc[(norm_counts['CHHpkb']>0) & (norm_counts["Annotation"].str.contains("_IGR_")==True)].to_csv("Norm_counts_called_Intergenic_072620_CHHnonzero.tsv",sep='\t',index=False)

cpg_gb_mean=norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)]['CPGpkb'].mean()
cpg_gb_std=norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)]['CPGpkb'].std()
cpg_p_mean=norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==True)]['CPGpkb'].mean()
cpg_p_std=norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==True)]['CPGpkb'].std()
cpg_t_mean=norm_counts.loc[(norm_counts["Annotation"].str.contains("_T")==True)]['CPGpkb'].mean()
cpg_t_std=norm_counts.loc[(norm_counts["Annotation"].str.contains("_T")==True)]['CPGpkb'].std()
cpg_i_mean=norm_counts.loc[(norm_counts["Annotation"].str.contains("_IGR_")==True)]['CPGpkb'].mean()
cpg_i_std=norm_counts.loc[(norm_counts["Annotation"].str.contains("_IGR_")==True)]['CPGpkb'].std()
chg_gb_mean=norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)]['CHGpkb'].mean()
chg_gb_std=norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)]['CHGpkb'].std()
chg_p_mean=norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==True)]['CHGpkb'].mean()
chg_p_std=norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==True)]['CHGpkb'].std()
chg_t_mean=norm_counts.loc[(norm_counts["Annotation"].str.contains("_T")==True)]['CHGpkb'].mean()
chg_t_std=norm_counts.loc[(norm_counts["Annotation"].str.contains("_T")==True)]['CHGpkb'].std()
chg_i_mean=norm_counts.loc[(norm_counts["Annotation"].str.contains("_IGR_")==True)]['CHGpkb'].mean()
chg_i_std=norm_counts.loc[(norm_counts["Annotation"].str.contains("_IGR_")==True)]['CHGpkb'].std()
chh_gb_mean=norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)]['CHHpkb'].mean()
chh_gb_std=norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)]['CHHpkb'].std()
chh_p_mean=norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==True)]['CHHpkb'].mean()
chh_p_std=norm_counts.loc[(norm_counts["Annotation"].str.contains("_P")==True)]['CHHpkb'].std()
chh_t_mean=norm_counts.loc[(norm_counts["Annotation"].str.contains("_T")==True)]['CHHpkb'].mean()
chh_t_std=norm_counts.loc[(norm_counts["Annotation"].str.contains("_T")==True)]['CHHpkb'].std()
chh_i_mean=norm_counts.loc[(norm_counts["Annotation"].str.contains("_IGR_")==True)]['CHHpkb'].mean()
chh_i_std=norm_counts.loc[(norm_counts["Annotation"].str.contains("_IGR_")==True)]['CHHpkb'].std()
cpg_nz_gb_mean=norm_counts.loc[(norm_counts['CPGpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)]['CPGpkb'].mean()
cpg_nz_gb_std=norm_counts.loc[(norm_counts['CPGpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)]['CPGpkb'].std()
cpg_nz_p_mean=norm_counts.loc[(norm_counts['CPGpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==True)]['CPGpkb'].mean()
cpg_nz_p_std=norm_counts.loc[(norm_counts['CPGpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==True)]['CPGpkb'].std()
cpg_nz_t_mean=norm_counts.loc[(norm_counts['CPGpkb']>0) & (norm_counts["Annotation"].str.contains("_T")==True)]['CPGpkb'].mean()
cpg_nz_t_std=norm_counts.loc[(norm_counts['CPGpkb']>0) & (norm_counts["Annotation"].str.contains("_T")==True)]['CPGpkb'].std()
cpg_nz_i_mean=norm_counts.loc[(norm_counts['CPGpkb']>0) & (norm_counts["Annotation"].str.contains("_IGR_")==True)]['CPGpkb'].mean()
cpg_nz_i_std=norm_counts.loc[(norm_counts['CPGpkb']>0) & (norm_counts["Annotation"].str.contains("_IGR_")==True)]['CPGpkb'].std()
chg_nz_gb_mean=norm_counts.loc[(norm_counts['CHGpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)]['CHGpkb'].mean()
chg_nz_gb_std=norm_counts.loc[(norm_counts['CHGpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)]['CHGpkb'].std()
chg_nz_p_mean=norm_counts.loc[(norm_counts['CHGpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==True)]['CHGpkb'].mean()
chg_nz_p_std=norm_counts.loc[(norm_counts['CHGpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==True)]['CHGpkb'].std()
chg_nz_t_mean=norm_counts.loc[(norm_counts['CHGpkb']>0) & (norm_counts["Annotation"].str.contains("_T")==True)]['CHGpkb'].mean()
chg_nz_t_std=norm_counts.loc[(norm_counts['CHGpkb']>0) & (norm_counts["Annotation"].str.contains("_T")==True)]['CHGpkb'].std()
chg_nz_i_mean=norm_counts.loc[(norm_counts['CHGpkb']>0) & (norm_counts["Annotation"].str.contains("_IGR_")==True)]['CHGpkb'].mean()
chg_nz_i_std=norm_counts.loc[(norm_counts['CHGpkb']>0) & (norm_counts["Annotation"].str.contains("_IGR_")==True)]['CHGpkb'].std()
chh_nz_gb_mean=norm_counts.loc[(norm_counts['CHHpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)]['CHHpkb'].mean()
chh_nz_gb_std=norm_counts.loc[(norm_counts['CHHpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==False) & (norm_counts["Annotation"].str.contains("_T")==False) & (norm_counts["Annotation"].str.contains("_IGR_")==False)]['CHHpkb'].std()
chh_nz_p_mean=norm_counts.loc[(norm_counts['CHHpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==True)]['CHHpkb'].mean()
chh_nz_p_std=norm_counts.loc[(norm_counts['CHHpkb']>0) & (norm_counts["Annotation"].str.contains("_P")==True)]['CHHpkb'].std()
chh_nz_t_mean=norm_counts.loc[(norm_counts['CHHpkb']>0) & (norm_counts["Annotation"].str.contains("_T")==True)]['CHHpkb'].mean()
chh_nz_t_std=norm_counts.loc[(norm_counts['CHHpkb']>0) & (norm_counts["Annotation"].str.contains("_T")==True)]['CHHpkb'].std()
chh_nz_i_mean=norm_counts.loc[(norm_counts['CHHpkb']>0) & (norm_counts["Annotation"].str.contains("_IGR_")==True)]['CHHpkb'].mean()
chh_nz_i_std=norm_counts.loc[(norm_counts['CHHpkb']>0) & (norm_counts["Annotation"].str.contains("_IGR_")==True)]['CHHpkb'].std()

print("CpGs/kb (average)\t"+"\t".join([str(cpg_gb_mean),str(cpg_p_mean),str(cpg_t_mean),str(cpg_i_mean)])+"\n")
print("CpGs/kb (standard deviation)\t"+"\t".join([str(cpg_gb_std),str(cpg_p_std),str(cpg_t_std),str(cpg_i_std)])+"\n")
print("CHGs/kb (average)\t"+"\t".join([str(chg_gb_mean),str(chg_p_mean),str(chg_t_mean),str(chg_i_mean)])+"\n")
print("CHGs/kb (standard deviation)\t"+"\t".join([str(chg_gb_std),str(chg_p_std),str(chg_t_std),str(chg_i_std)])+"\n")
print("CHHs/kb (average)\t"+"\t".join([str(chh_gb_mean),str(chh_p_mean),str(chh_t_mean),str(chh_i_mean)])+"\n")
print("CHHs/kb (standard deviation)\t"+"\t".join([str(chh_gb_std),str(chh_p_std),str(chh_t_std),str(chh_i_std)])+"\n")

print("Non-zero CpGs/kb (average)\t"+"\t".join([str(cpg_nz_gb_mean),str(cpg_nz_p_mean),str(cpg_nz_t_mean),str(cpg_nz_i_mean)])+"\n")
print("Non-zero CpGs/kb (standard deviation)\t"+"\t".join([str(cpg_nz_gb_std),str(cpg_nz_p_std),str(cpg_nz_t_std),str(cpg_nz_i_std)])+"\n")
print("Non-zero CHGs/kb (average)\t"+"\t".join([str(chg_nz_gb_mean),str(chg_nz_p_mean),str(chg_nz_t_mean),str(chg_nz_i_mean)])+"\n")
print("Non-zero CHGs/kb (standard deviation)\t"+"\t".join([str(chg_nz_gb_std),str(chg_nz_p_std),str(chg_nz_t_std),str(chg_nz_i_std)])+"\n")
print("Non-zero CHHs/kb (average)\t"+"\t".join([str(chh_nz_gb_mean),str(chh_nz_p_mean),str(chh_nz_t_mean),str(chh_nz_i_mean)])+"\n")
print("Non-zero CHHs/kb (standard deviation)\t"+"\t".join([str(chh_nz_gb_std),str(chh_nz_p_std),str(chh_nz_t_std),str(chh_nz_i_std)])+"\n")