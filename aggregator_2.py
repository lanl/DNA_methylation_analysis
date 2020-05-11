# coding: utf-8
for i in ["7-25","7-26","7-27","7-28","7-31"]:
    df1=pd.read_csv("Calls/"+i+"/CpG_result_table.txt_agt20.txt",sep='\t')
    df2=pd.read_csv("Calls/"+i+"/CHG_result_table.txt_agt20.txt",sep='\t')
    df3=pd.read_csv("Calls/"+i+"/CHH_result_table.txt_agt20.txt",sep='\t')
    all_df=pd.concat([df1,df2,df3])
    all_df.to_csv("Calls/"+i+"/all_result_agt20_table.txt",sep='\t')

