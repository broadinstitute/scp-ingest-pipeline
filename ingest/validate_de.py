import pandas as pd
import itertools
import numpy as np
file_path = "/Users/mvelyuns/scp-ingest-pipeline/ingest/pairwise.csv"
data = pd.read_csv(file_path)

genes = data.iloc[:, 0].tolist()
rest = data[data.columns[1:]]
col = list(rest.columns)

col_values = {}
for i in col:
    col_values[i] = rest[i].tolist()

#note: my initial files had pval, qval, log2fc. David's files have qval, mean, log2fc. For the purposes of this validation I will be using his column values/formatting. If we choose to remove mean again this code can be changed

#question: should user also input type names to not rely solely on parsing?

#below parsing relies on this format: 'type_0'_'type_1'_qval
#necessary to include quotes
clean_val = []

#cleans column header of unnecessary characters and isolates column name, returns clean_val containing only cleaned values
for i in col:
    str = i.split("'")
    col_names = []
    for j in str:
        if (j != '') and (j != "_"):
            col_names.append(j.strip("_"))
    clean_val.append(col_names)

types = []
qual = []

#isolate the types and values in submitted file
#expected format: 
#types: ['type_0', 'type_1', 'type_2', 'type_3']
#qual: ['qval', 'log2fc', 'mean']
for i in clean_val:
    if i[0] not in types:
        types.append(i[0])
    if i[1] not in types:
        types.append(i[1])
    if i[2] not in qual:
        qual.append(i[2])

#verify all pairwise combinations are present and have each metric (qval, log2fc, mean)
#correct pairwise combinations:
correct_pairs = list(itertools.combinations(types,2))

correct_list = []
for i in correct_pairs:
    q1 = [i[0], i[1], qual[0]]
    q2 = [i[0], i[1], qual[1]]
    q3 = [i[0], i[1], qual[2]]
    correct_list.append(q1)
    correct_list.append(q2)
    correct_list.append(q3)

if correct_list != clean_val:
    raise ValueError("incorrect inputted values")


#create individual files
#desired format:
#for ex, if we have 
#'type_0'_'type_1'_qval
#'type_0'_'type_1'_log2fc
#'type_0'_'type_1'_mean
#final format should have type 0 type 1 in the title, and genes, qval, log2fc, mean as columns
names_dict = {}
all_group = []
for i in range(len(types)):
    curr_type = types[i]
    type_group = []

    for j in range(len(clean_val)):
        curr_val = clean_val[j][0]
        comp_val = clean_val[j][1]
        file_naming = "'" + curr_val + "'" + "_" + "'" + comp_val + "'"
        names_dict[file_naming] = []
        real_title = col[j]
        if curr_type == curr_val:
            type_group.append(real_title)

    all_group.append(type_group)

all_group_fin = [ele for ele in all_group if ele != []]
    
grouped_lists = []

for i in all_group_fin:
    for j in range(0, len(i), 3):
        x = j
        grouped_lists.append(i[x:x+3])

for i in names_dict:
    for j in grouped_lists:
        for k in j:
            if i in k:
                names_dict[i].append(k)

#Now we have all the columns grouped in lists by pairwise comparison, with qval, log2fc, mean
#have to pair with corresponding gene for that row
#dictionary format:
# comparison name: [[gene, qval, log2fc, mean] [gene, qval, log2fc, mean] etc...] 
keys = names_dict.keys()
file_d = dict.fromkeys(keys, [])

def check_type(name):
    for i in names_dict:
        if name in names_dict[i]:
            return i

individual_vals = {}
test = []
for i in grouped_lists:
    list_i = []
    for j in i:
        if "qval" in j:
            n = "qval"
        if "log2fc" in j:
            n = "log2fc"
        if "mean" in j:
            n = "mean"
        
        f_name = check_type(j)
        col_v = rest[j].tolist()
        list_i.append(col_v)

    file_d[f_name] = genes, list_i[0], list_i[1], list_i[2]

fin_column_names = ["genes", "qval", "log2fc", "mean"]

for i in file_d: 
    arr = np.array(file_d[i])
    t_arr = arr.transpose()
    inner_df = pd.DataFrame(data = t_arr, columns = fin_column_names)
    inner_df.to_csv("ingest/{}.tsv".format(i), sep ='\t')

#final result: individual files for each comparison