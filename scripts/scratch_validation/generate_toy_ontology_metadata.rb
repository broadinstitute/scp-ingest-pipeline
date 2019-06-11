#!/ usr/bin/env ruby

require 'optparse'
require 'securerandom'

VALID_CELL_TYPES = %w[CL_0002419 CL_0000863 CL_0000875]
VALID_DISEASES = %w[MONDO_0018882 MONDO_0021042 MONDO_0007455]
VALID_ORGANS = %w[UBERON_0001264 UBERON_0002046 UBERON_0000948]
VALID_DATA = []

error_frequency = 10
max_lines = 100
delimiter = "\t"

# options parsing
OptionParser.new do |opts|
    opts.banner = "Usage: ruby generate_toy_ontology_metadata.rb [options]"

    opts.on("-e", "--error-frequency ERROR", Integer, "Frequency of invalid ontology IDs") do |frequency|
        error_frequency = frequency
    end
    
    opts.on("-n", "--num-lines LINES", Integer, "Total number of lines to generate") do |lines|
        max_lines = lines
    end

    opts.on("-d", "--delimiter", String, "Output file delimiter") do |delim|
        delimiter = "#{delim}".gsub(/ /, '')
    end

    opts.on("-h", "--help", "Prints this help") do
        puts "\n#{opts}\n"
        exit
    end
end.parse!

if error_frequency < 10
    puts "Error frequency is too high: #{error_frequency}, defaulting back to 10"
end

newfile = File.new('cell_metadata_toy.txt', 'w+')
headers = ["Cell ID", "Cell Type", "Disease", "Organ"]
newfile.write headers.join(delimiter) + "\n"
error_counter = 0
puts "writing #{max_lines} of synthetic data at an error rate of 1 error every #{error_frequency} values"
1.upto(max_lines) do |line_number|
    cell_name = "CELL_#{line_number}"
    line = "#{cell_name}#{delimiter}"
    cell_type = VALID_CELL_TYPES.sample
    disease = VALID_DISEASES.sample
    organ = VALID_ORGANS.sample
    [cell_type, disease, organ].each do |value|
        error_counter += 1
        if error_frequency > 0 && error_counter == error_frequency
            value = "#{SecureRandom.alphanumeric(4).upcase}_000#{(SecureRandom.rand * 5000 + 1000).floor}"
            error_counter = 0
        end
        line += value + delimiter
    end
    newfile.write line + "\n"
end
puts "complete!"
newfile.close