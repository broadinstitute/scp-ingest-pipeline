import pandas as pd
import numpy as np
import sys

#input file name as CLI
#python ingest/validate_de.py /Users/mvelyuns/scp-ingest-pipeline/ingest/pairwise.csv

if len(sys.argv) > 1:
    file_path = sys.argv[1]
else:
    file_path = "/Users/mvelyuns/scp-ingest-pipeline/ingest/pairwise.csv"   

#only necessary for this specific example since we plan to accept all other files in long format to begin with 
def convert_wide_to_long(data):
    data.rename(columns={ data.columns[0]: "genes"}, inplace=True)
    ls = list(data.columns)[1:]
    long = pd.melt(data, id_vars='genes', var_name = "comparisons", value_vars= ls)
    return long


#TYPICAL USE: convert from long format (intended input) to wide format, since parsing is easier with wide format in this case
def convert_long_to_wide(data):
    wide = pd.pivot(data, index="genes", columns="comparisons", values="value")
    wide.reset_index(inplace=True)
    return wide

#left to do:
# 1. make necessary changes to support long format 
# 2. DONE -- remove parts that need all pairwise + error message
# 3. DONE -- file name as CLI
# 4. DONE -- make functions
# 5. random small stuff
# 6. DONE -- add main function
# 7. DONE -- fix file name

#outputs full list of column names
#outputs full list of genes
#outputs "rest" which is columns isolated from genes
#col_values outputs dict with the comparison characteristic and all the numeric data in an array, for ex 'type_0'_'type_1'_qval: [0, 1, 2]

def get_data_by_col(data):  
    genes = data.iloc[:, 0].tolist()
    rest = data[data.columns[1:]]
    col = list(rest.columns)
    col_values = {}
    for i in col:
        col_values[i] = rest[i].tolist()
    return col, genes, rest

#note: my initial files had pval, qval, log2fc. David's files have qval, mean, log2fc. For the purposes of this validation I will be using his column values/formatting. If we choose to remove mean again this code can be changed

#question: should user also input type names to not rely solely on parsing?

#below parsing relies on this format: 'type_0'_'type_1'_qval
#necessary to include quotes

#cleans column header of unnecessary characters and isolates column name, returns clean_val containing only cleaned values

#qual is the properties to measure by, in the example's case it is ['qval', 'log2fc', 'mean']
    #parse from columns rather than assume these are the three we are using
#types is the types being compared, in this case ['type_0', 'type_1', 'type_2', 'type_3']

def get_types_and_properties(col):
    clean_val = []

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
    return types, clean_val, qual

#create individual files
#desired format:
#for ex, if we have 
#'type_0'_'type_1'_qval
#'type_0'_'type_1'_log2fc
#'type_0'_'type_1'_mean
#final format should have type 0 type 1 in the title, and genes, qval, log2fc, mean as columns

def check_type(names_dict, name):
    for i in names_dict:
        if name in names_dict[i]:
            return i

def generate_individual_files(col, genes, rest, types, clean_val, qual):
    names_dict = {}
    all_group = []
    for i in range(len(types)):
        curr_type = types[i]
        type_group = []

        for j in range(len(clean_val)):
            curr_val = clean_val[j][0]
            comp_val = clean_val[j][1]
            file_naming = f"{curr_val}_{comp_val}"
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
                k = k.replace("'", "")
                if i in k:
                    names_dict[i].append(k)
    
#Now we have all the columns grouped in lists by pairwise comparison, with qval, log2fc, mean
#have to pair with corresponding gene for that row
#dictionary format:
# comparison name: [[gene, qval, log2fc, mean] [gene, qval, log2fc, mean] etc...] 
    keys = names_dict.keys()
    file_d = dict.fromkeys(keys, [])

    for i in grouped_lists:
        list_i = []
        for j in i:        
            new_j = j.replace("'", "")
            f_name = check_type(names_dict, new_j)
            col_v = rest[j].tolist()
            list_i.append(col_v)

        file_d[f_name] = genes, list_i[0], list_i[1], list_i[2]

    qual.insert(0, "genes")

    print(qual)
    for i in file_d: 
        arr = np.array(file_d[i])
        t_arr = arr.transpose()
        print(t_arr)
        inner_df = pd.DataFrame(data = t_arr, columns = qual)
        inner_df.to_csv("ingest/{}.tsv".format(i), sep ='\t')

#final result: individual files for each comparison

if __name__ == '__main__':
    data = pd.read_csv(file_path)
    long_format = convert_wide_to_long(data)
    wide_format = convert_long_to_wide(long_format)

    col = get_data_by_col(wide_format)[0]
    genes = get_data_by_col(wide_format)[1]
    rest = get_data_by_col(wide_format)[2]
    types = get_types_and_properties(col)[0]
    clean_val = get_types_and_properties(col)[1]
    qual = get_types_and_properties(col)[2]

    generate_individual_files(col, genes, rest, types, clean_val, qual)