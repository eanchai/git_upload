import pickle

with open('/home/ubuntu/NGS_data/pickle/fastqc_result.pickle','rb') as file:
            data = pickle.load(file)
            
for sample_id in data.keys():
    for read_id in data[sample_id].keys():
        update_table = data[sample_id][read_id]['sge']
        update_table.insert(1, 'SEQ_DATE', 'date')
        update_table.insert(2, 'FC_ID', 'fc_id')
        update_table.insert(3, 'SAMPLE_ID', sample_id)
        update_table.insert(5, 'IDX', 0)
        update_table.insert(4, 'FASTQ_TYPE', read_id)
        update_table['CREATE_DATE'] = 'null'
        update_table = update_table.rename(columns={'Base':"BASE"})
        
update_table

for sample_id in data.keys():
    for read_id in data[sample_id].keys():
        update_table = data[sample_id][read_id]['bsq']
        update_table.insert(1, 'SEQ_DATE', 'date')
        update_table.insert(2, 'FC_ID', 'fc_id')
        update_table.insert(3, 'SAMPLE_ID', sample_id)
        update_table.insert(5, 'IDX', 0)
        update_table.insert(4, 'FASTQ_TYPE', read_id)
        update_table['CREATE_DATE'] = 'null'
update_table = update_table.rename(columns={'Base':"BASE"})
update_table = update_table.rename(columns={'10th Percentile':'TENTH_PERCENTILE'})
update_table = update_table.rename(columns={'90th Percentile':'NINETIETH_PERCENTILE'})

        
update_table

for sample_id in data.keys():
    for read_id in data[sample_id].keys():
        update_table = data[sample_id][read_id]['bsc']
        update_table.insert(1, 'SEQ_DATE', 'date')
        update_table.insert(2, 'FC_ID', 'fc_id')
        update_table.insert(3, 'SAMPLE_ID', sample_id)
        update_table.insert(5, 'IDX', 0)
        update_table.insert(4, 'FASTQ_TYPE', read_id)
        update_table['CREATE_DATE'] = 'null'
        update_table = update_table.rename(columns={'Base':"BASE"})
        

 for _, col in update_table.columns.tolist().replace(" ","_")
 
col_list = []
for col in update_table.columns.tolist():
    new_col = col.replace(" ","_")
    col_list.append(new_col)
    
print(col_list)

update_table.columns = col_list