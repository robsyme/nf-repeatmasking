#!/usr/bin/env ruby
require 'bio'
require 'pp'

HEADER_MATCH = /^\d+ \d+\.?\d* \d+\.?\d* \d+\.?\d* (?<refSeqid>\S+) (?<refStart>\d+) (?<refStop>\d+) \(\d+\) (?<complement>C)? ?(?<repSeqid>\S+) \(?(?<left>\d+)\)? (?<middle>\d+) \(?(?<right>\d+)\)? (?<matchID>\S+) \d+$/
ALIGNMENT_MATCH = /^[C ] (?<seqid>\S+)\s+(?<start>\d+) (?<bases>\S+) (?<end>\d+)/


class Alignment
  attr_reader :refSeq, :repSeq
  def initialize(reference, repeat)
    refMatch = reference.match(ALIGNMENT_MATCH)
    repMatch = repeat.match(ALIGNMENT_MATCH)

    @repStart =

    @refSeq = refMatch['bases']
    @repSeq = repMatch['bases']
  end

  def +(alignment)
    @refSeq += alignment.refSeq
    @repSeq += alignment.repSeq
    self
  end

  def repBasesWithoutGaps
    Bio::Sequence::NA.new(@refSeq
      .chars
      .zip(@repSeq.chars)
      .find_all{|refBase, repBase| repBase != '-'}
      .map{|refBase, repBase| refBase}
      .join
      .upcase)
  end
end

class AlignMatch
  attr_reader :queryID, :matchID
  def initialize(lines)
    lines = lines.reject{|line| line == ''}
    header = lines.shift
    match = header.match(HEADER_MATCH)
    @hitSeqid = match['refSeqid']
    @hitStart = match['refStart'].to_i
    @hitStop = match['refStop'].to_i
    @complement = match['complement'] == "C"
    @queryStop = match['middle'].to_i
    if @complement
      @queryStart = match['right'].to_i
    else
      @queryStart = match['left'].to_i
    end
    @queryID = match['repSeqid']
    @matchID = match['matchID']

    @alignment = lines
                   .find_all{|line| line =~ ALIGNMENT_MATCH}
                   .each_slice(2)
                   .map{|a| Alignment.new(*a)}
                   .inject(:+)
  end

  def newHit
    hit = Bio::Sequence::NA.new("-" * (@queryStart - 1))
    hit += @complement ? @alignment.repBasesWithoutGaps.complement : @alignment.repBasesWithoutGaps
    hit
  end
end

queries = Hash[Bio::FlatFile.open(ARGV.pop).map{|e| [e.entry_id, e.naseq]}]

alignments = Hash.new{|hash,key| aln = Bio::Alignment.new([]); aln.add_seq(queries[key], key); hash[key] = aln}

File.open(ARGV.shift)
  .map{|line| line.chomp}
  .slice_before(HEADER_MATCH)
  .map{|a| AlignMatch.new(a)}
  .each{|am| alignments[am.queryID].add_seq(am.newHit, am.matchID)}

alignments.each do |repSeqID, alignment|
  File.open("alignment.#{repSeqID.split("#").first}.fasta",'w') do |out|
    out.puts alignment.normalize.output_fasta.lines.map{|line| line =~ /^>/ ? line : line.upcase}
  end
end
