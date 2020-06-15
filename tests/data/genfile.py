import random
num_cells = 1_000_000
num_genes = 300
with open('test_db_data.txt', 'w') as exp_file:
    header_row = 'GENE'
    for i in range(num_cells):
        header_row += (",cell_foo_" + str(i))
    exp_file.write(header_row)
    exp_file.write('\n')
    for i in range(num_genes):
        gene_row = 'gf' + str(i)
        is_expressed = random.randint(1,3) == 1
        for j in range(num_cells):
            col_value = ',0'
            if is_expressed and random.randint(1,7) == 1:
                col_value = ',' + str(random.uniform(0,6))
            gene_row += col_value

        exp_file.write(gene_row)
        exp_file.write('\n')

