#!/usr/bin/env ruby
require 'bio'
require 'pp'

BASES = %w{A C G T}

alignment = Bio::Alignment::OriginalAlignment.readfiles(ARGV.shift)
sites = Enumerator.new{|y| alignment.each_site{|site| y << site }}

deripped_consensus = Enumerator.new do |yielder|
  alignment
    .consensus_string(0.1, {:gap_mode => -1})
    .chars
    .zip(sites)
    .each do |consensus, site|
    riplike = FALSE
    counts = Hash[BASES.zip(BASES.map{|base| site.count{|nuc| nuc.upcase == base}})]
    
    if counts["T"] == counts.max && counts["C"] > 0
      yielder << "C"
    elsif counts["A"] == counts.max && counts["G"] > 0
      yielder << "G"
    elsif consensus == "?"
      yielder << "N"
    else
      yielder << consensus
    end
  end
end

puts Bio::Sequence::NA.new(deripped_consensus.to_a.join).upcase.to_fasta(alignment.keys.first, 80)



