#! /usr/bin/env ruby

# firestore_query.rb
#
# PREREQUISITES
#
# Ruby >= 2.3
# Gems: google-cloud-firestore, parallel
# Files: GCP service account JSON credentials for Firestore access, dense gene expression matrix

require 'optparse'
require 'rubygems'
require 'google/cloud/firestore'
require 'parallel'

credentials_filepath = 'broad-singlecellportal-staging-firestore.json'
matrix_filepath = 'expression_matrix_example.txt'
project_id = 'broad-singlecellportal-staging'
max_lines = nil
num_threads = 10
verbose = false

# options parsing
OptionParser.new do |opts|
    opts.banner = "Usage: ruby firestore_query.rb [options]"

    opts.on("-c", "--credentials CREDENTIALS", String, "Filepath of service account JSON credentials file for Firestore access") do |creds|
        credentials_filepath = creds.strip
    end

    opts.on("-m", "--matrix-file MATRIX_FILE", String, "Filepath of input gene expression matrix to ingest and query") do |matrix|
        matrix_filepath = matrix.strip
    end

    opts.on("-p", "--project-id PROJECT", String, "Google Cloud Platform project id") do |project|
        project_id = project.strip
    end

    opts.on("-l", "--max-lines LINES", String, "Limit of number of lines to ingest from expression matrix") do |lines|
        max_lines = lines.strip.to_i
    end

    opts.on("-t", "--threads THREADS", String, "Number of parallel threads when ingesting/deleting data in Firestore") do |threads|
        num_threads = threads.strip.to_i
    end

    opts.on("-v", "--verbose", "Print detailed ingest/validation messages") do
        verbose = true
    end

    opts.on("-h", "--help", "Prints this help") do
        puts "\n#{opts}\n"
        exit
    end
end.parse!

# create test collection to use
collection_name = "test-genes-#{SecureRandom.uuid}"

puts ""
puts "credentials file: #{credentials_filepath}"
puts "matrix file: #{matrix_filepath}"
puts "Google project ID: #{project_id}"
puts "Maximum line limit: #{max_lines}"
puts "Number of threads: #{num_threads}"
puts "Using test collection name of #{collection_name}"
puts ""

start_time = Time.now

# instantiate Firestore client
if File.exists?(credentials_filepath)
    firestore = Google::Cloud::Firestore.new(
        project_id: project_id,
        credentials: credentials_filepath
    )
else
    puts "Credentials file not found at: #{credentials_filepath}"
    exit(1)
end

# open raw matrix data and store list of cell names
matrix_data = {}
if File.exists?(matrix_filepath)
    matrix_file = File.open(matrix_filepath, 'rb')
    source_file_name = File.basename(matrix_filepath)
else
    puts "Matrix file not found at: #{matrix_filepath}"
    exit(1)
end

cell_names = matrix_file.readline.split(/[\t,]/).map {|cell| cell.strip.to_sym}
if cell_names.first.empty? || cell_names.first.downcase == 'gene'
    cell_names.shift
end

# read matrix and store significant (non-zero) data
print "ingesting #{matrix_filepath}... "
total_genes = 0
while !matrix_file.eof && (max_lines.nil? || matrix_file.lineno <= max_lines)
    # read & transform data into structure to validate with later
    data = matrix_file.readline.split(/[\t,]/).map(&:strip)
    gene_name = data.shift
    scores = data.map(&:to_f)
    significant_scores = {}
    scores.each_with_index do |score, index|
        unless score == 0.0
            significant_scores[cell_names[index]] = score
        end
    end
    matrix_data[gene_name] = significant_scores
    total_genes += 1
    if matrix_file.lineno % 1000 == 0
        puts "read #{matrix_file.lineno} genes " if verbose
    end
end
puts "complete! #{total_genes} genes read"

write_start = Time.now
# using significant data, write to Firestore
print "writing to Firestore in parallel... "
total_batches = total_genes < 100 ? 1 : total_genes / 100
matrix_batches = matrix_data.to_a.each_slice(100)
Parallel.each_with_index(matrix_batches, in_threads: num_threads) do |data_batch, batch_index|
    puts "writing batch #{batch_index + 1} of #{matrix_batches.size}" if verbose
    firestore.batch do |fs_batch|
        data_batch.to_h.each do |gene, scores|
            gene_data = {
                name: gene,
                source_file_name: source_file_name
            }
            scores_data = {
                scores: scores
            }
            fs_batch.set("#{collection_name}/#{gene}", gene_data)
            fs_batch.set("#{collection_name}/#{gene}/expression_scores/scores", scores_data)
        end
    end
end
write_end = Time.now
write_minutes, write_seconds = (write_end - write_start).divmod 60.0
puts "Firestore write complete! ingest runtime: #{write_minutes} minutes, #{write_seconds.round} seconds"
print "validating writes against source data... "

# validate all genes are queryable, and expression data returned is identical
valid_genes = 0
query_start = Time.now
Parallel.each_with_index(matrix_data, in_threads: num_threads) do |(gene, scores), index|
    if verbose
        puts ""
        puts "validating #{gene}"
    end
    query = firestore.col(collection_name).where(:name, :==, gene).get
    if query.nil?
        puts "### #{gene} not found!!" if verbose
    else
        gene_ref = query.first.ref
        gene_data = gene_ref.col('expression_scores').doc('scores').get.data[:scores]
        if scores == gene_data
            valid_genes += 1
            if verbose
                puts "#{gene} is valid!"
                puts "source: #{scores.sort.to_h}"
                puts "query:  #{gene_data.sort.to_h}"
            end
        else
            puts "#{gene} is invalid!" if verbose
            scores.each do |cell, score|
                unless score == gene_data[cell]
                    puts "### #{cell} is incorrect, expected #{score} but found #{gene_data[cell]} ###" if verbose
                end
            end
        end
        puts "" if verbose
    end
    if (index + 1) % 1000 == 0
        puts "processed #{index + 1} genes"
    end
end
query_end = Time.now
puts "validation complete! #{valid_genes} out of #{total_genes} genes validated"
query_minutes, query_seconds = (query_end - query_start).divmod 60.0
puts "validation runtime: #{query_minutes} minutes, #{query_seconds.round} seconds"

# clean up all documents, print statistics
print "cleaning up... "
gene_documents = firestore.col(collection_name).get
delete_start = Time.now
Parallel.map(gene_documents, in_threads: num_threads) do |gene_document|
    ref = gene_document.ref
    ref.col('expression_scores').get do |sub_col|
        sub_ref = sub_col.ref
        sub_ref.delete
    end
    ref.delete
end
delete_end = Time.now
end_time = Time.now
delete_minutes, delete_seconds = (delete_end - delete_start).divmod 60.0
puts "Firestore cleanup complete! deletion runtime: #{delete_minutes} minutes, #{delete_seconds.round} seconds"
total_minutes, total_seconds = (end_time - start_time).divmod 60.0
puts "total runtime: #{total_minutes} minutes, #{total_seconds.round} seconds"
exit_status = total_genes == valid_genes ? 0 : 1
exit(exit_status)