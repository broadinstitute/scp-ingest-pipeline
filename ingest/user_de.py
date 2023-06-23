import pandas as pd
import numpy as np
import sys
import csv

cluster_name = "All_Cells_UMAP"
clean_annotation = "General_Celltype"
annot_scope = "study"
method = "wilcoxon"

# input file name as CLI
# example:
# python ingest/user_de.py /Users/mvelyuns/scp-ingest-pipeline/ingest/pairwise.csv

if len(sys.argv) > 1:
    file_path = sys.argv[1]
else:
    # file_path = "ingest/pairwise.csv"   
    # file_path = "ingest/pairwise_and_rest.csv"   
    # file_path = "ingest/test_one_vs_rest.csv"  
    file_path = "ingest/test_long_format.csv"

# only necessary for this specific example since we plan to accept all other files in long format to begin with 
def convert_wide_to_long(data):
    data.rename(columns={ data.columns[0]: "genes"}, inplace=True)
    ls = list(data.columns)[1:]
    long = pd.melt(data, id_vars='genes', var_name = "comparisons", value_vars= ls)
    return long

def convert_long_to_wide(data):
    """
    TYPICAL USE: convert from long format (intended input) to wide format, since parsing is easier with wide format in this case
    """
    data["combined"] = f"{data['group']}--{data['comparison_group']}"
    frames = []
    metrics = ["qval", "mean", "log2fc"]
    for metric in metrics:
        wide_metric = pd.pivot(data, index="genes", columns="combined", values=metric)
        wide_metric = wide_metric.add_suffix(f"--{metric}")
        frames.append(wide_metric)

    result = pd.concat(frames, axis=1, join='inner')
    result.columns.name = ' ' 
    result = result.reset_index()
    return result

def get_data_by_col(data):  
    """
    outputs full list of column names
    outputs full list of genes
    outputs "rest" which is columns isolated from genes
    col_values outputs dict with the comparison characteristic and all the numeric data in an array, for ex 'type_0'_'type_1'_qval: [0, 1, 2]
    """
    genes = data.iloc[:, 0].tolist()
    rest = data[data.columns[1:]]
    col = list(rest.columns)

    # split col into two lists: pairwise, one_vs_rest
    split_values = {"one_vs_rest": [], "pairwise": []}
    for i in col:
        if "rest" in i:
            split_values["one_vs_rest"].append(i)
        else:
            split_values["pairwise"].append(i)

    return col, genes, rest, split_values
# note: my initial files had pval, qval, log2fc. David's files have qval, mean, log2fc. For the purposes of this validation I will be using his column values/formatting.

def get_groups_and_properties(col):
    """
    takes in column names 

    below parsing relies on this format: 'type_0'_'type_1'_qval - necessary to include quotes

    cleans column header of unnecessary characters and isolates column name, returns clean_val containing only cleaned values

    qual is the properties to measure by, in the example's case it is ['qval', 'log2fc', 'mean']
     -- important note: parse from columns rather than assume these are the three we are using

     groups is the group being compared, in this case ['type_0', 'type_1', 'type_2', 'type_3']
    """
    clean_val = []

    for i in col:
        str = i.split("--")
        col_names = []
        for j in str:
            j = j.replace("'", "")
            if (j != '') and (j != "_"):
                col_names.append(j.strip("_"))
        clean_val.append(col_names)

    groups = []
    qual = []

# isolate the groups and values in submitted file
# expected format: 
# groups: ['group_0', 'group_1', 'group_2', 'group_3']
# qual: ['qval', 'log2fc', 'mean']
    for i in clean_val:
        if i[0] not in groups:
            groups.append(i[0])
        if i[1] not in groups:
            groups.append(i[1])
        if i[2] not in qual:
            qual.append(i[2])

    #TODO: Report this error to Sentry
    if ("log2fc" not in qual) or ("qval" not in qual):
        raise Exception("Comparisons must include at least log2fc and qval to be valid")
    return groups, clean_val, qual

def check_group(names_dict, name):
    """
    helper function to return the comparison based on the comparison with the value

    takes in dictionary that stores all names, and given name

    for example, input type_2_type_3_mean returns type_2_type_3
    """
    for i in names_dict:
        if name in names_dict[i]:
            return i
        

def sort_comparison(ls):
    """
    Naturally sort groups in a pairwise comparison; specially handle one-vs-rest 
    this should take in only a list of two, ex (type1 type2)
    
    https://en.wikipedia.org/wiki/Natural_sort_order

    """
    if any(i.isdigit() for i in ls):
        sorted_arr = sorted(ls, key=lambda x: int("".join([i for i in x if i.isdigit()])))
        return sorted_arr
    elif "rest" == ls[1]:
        return ls
    elif "rest" == ls[0]:
        return [ls[1], ls[0]]
    else:
        return sorted(ls)    

def generate_individual_files(col, genes, rest, groups, clean_val, qual):
    """
    create individual files for each comparison, pairwise or rest, with all the metrics being used (ex qval, log2fc, mean)
    desired format:
    for ex, if we have 
    'type_0'--'type_1'--qval
    type_0'--'type_1'--log2fc
    'type_0'--'type_1'--mean
    final format should have type 0 type 1 in the title, and genes, qval, log2fc, mean as columns
    """
    for i in clean_val:
        val_to_sort = [i[0], i[1]]
        sorted_list = sort_comparison(val_to_sort)
        i[0], i[1] = sorted_list

    names_dict = {}
    all_group = []
    for i in range(len(groups)):
        curr_group = groups[i]
        type_group = []

        for j in range(len(clean_val)):
            curr_val = clean_val[j][0]
            comp_val = clean_val[j][1]
            file_naming = f"{curr_val}--{comp_val}"
            names_dict[file_naming] = []
            real_title = col[j]
            if curr_group == curr_val:
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

# Now we have all the columns grouped in lists by pairwise comparison, with qval, log2fc, mean
# have to pair with corresponding gene for that row
# dictionary format:
# comparison name: [[gene, qval, log2fc, mean] [gene, qval, log2fc, mean] etc...] 
    keys = names_dict.keys()
    file_d = dict.fromkeys(keys, [])

    for i in grouped_lists:
        list_i = []
        for j in i:        
            f_name = check_group(names_dict, j)
            col_v = rest[j].tolist()
            list_i.append(col_v)

        file_d[f_name] = genes, list_i[0], list_i[1], list_i[2]

    qual.insert(0, "genes")

    final_files_to_find = []
    for i in file_d: 
        arr = np.array(file_d[i])
        t_arr = arr.transpose()
        inner_df = pd.DataFrame(data = t_arr, columns = qual)

        if "rest" in i:
            i = i.split("--")[0]

        tsv_name = f'ingest/{cluster_name}--{clean_annotation}--{i}--{annot_scope}--{method}.tsv'

        inner_df.to_csv(tsv_name, sep ='\t')
        #final_files_to_find.append("ingest/{}.tsv".format(i))

    #return final_files_to_find

# final result: individual files for each comparison

def generate_manifest(clean_val, clean_val_p, qual):
    """
    create manifest file of each comparison in the initial data
    if the comparison is with rest, rest is omitted and just the type is written
    """
    file_names_one_vs_rest = []
    for i in clean_val:
        clean_comparison = '--'.join([ x for x in i if x != "rest" and x not in qual])
        if clean_comparison not in file_names_one_vs_rest:
            file_names_one_vs_rest.append(clean_comparison)

    file_names_pairwise = []
    for i in clean_val_p:
        values = [i[0], i[1]]
        if values not in file_names_pairwise:
            file_names_pairwise.append(values)


    with open('ingest/manifest.tsv', 'w', newline='') as f:
        tsv_output = csv.writer(f, delimiter='\t')
        if len(file_names_one_vs_rest) != 0:
            for value in range(len(file_names_one_vs_rest)):
                tsv_output.writerow([file_names_one_vs_rest[value]])
        
        if len(file_names_pairwise) != 0:
            for value in range(len(file_names_pairwise)):
                tsv_output.writerow([file_names_pairwise[value][0], file_names_pairwise[value][1],])

if __name__ == '__main__':
    clean_val = []
    clean_val_p = []
    qual = []
    qual_p = []

    data = pd.read_csv(file_path)
    wide_format = convert_long_to_wide(data)
    print(wide_format)
    data_by_col = get_data_by_col(wide_format)
    col, genes, rest, split_values = data_by_col
    pairwise = split_values["pairwise"]
    one_vs_rest = split_values["one_vs_rest"]

    if len(one_vs_rest) != 0:
        groups_and_props = get_groups_and_properties(one_vs_rest)
        groups, clean_val, qual = groups_and_props
        generate_individual_files(one_vs_rest, genes, rest, groups, clean_val, qual)
            
    if len(pairwise) != 0:
        groups_and_props_p = get_groups_and_properties(pairwise)
        groups_p, clean_val_p, qual = groups_and_props_p
        generate_individual_files(pairwise, genes, rest, groups_p, clean_val_p, qual)

    generate_manifest(clean_val, clean_val_p, qual)

    
    