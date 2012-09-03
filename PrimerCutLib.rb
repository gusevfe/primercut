require 'optparse'
require 'bio'
require 'bio-samtools'
require 'ostruct'

Result = Struct.new(:aln, :with_T, :cut, :drop)

class Array
  def Array.load_regions(s)
    r = Array.new

    if File.exist?(s) then
      r = File.readlines(s)
    else
      r = s.split(",")
    end

    r = r.map { |x| Region.new(x) }

    return r
  end
end

class Region
  attr_accessor :chr, :start, :finish

  def initialize(s)
    begin
      raise if s !~ /^\s*[^:]+:\d+-\d+\s*$/
      @chr = s[/([^:]*):/, 1]
      @start = s[/([^:]*):([^-]*)-/, 2].to_i
      @finish = s[/([^:]*):([^-]*)-(.*)$/, 3].to_i
    rescue
      begin
        raise if s !~ /^\s*[^\s]+\s+\d+\s+\d+\s*$/
        s = s.split
        @chr = s[0]
        @start = s[1].to_i
        @finish = s[2].to_i
      rescue
        raise "Failed to parse #{s.inspect} as a region."
      end
    end

  end

  def to_s
    "#{chr}:#{start}-#{finish}"
  end
end

class Bio::DB::Tag
  def ==(x)
    self.tag == x.tag && self.type == x.type && self.value == x.value
  end
end

class Bio::DB::Alignment
  def sam
     data = []
     data << self.qname 
     data << self.flag  
     data << self.rname 
     data << self.pos   
     data << self.mapq  
     data << self.cigar 
     data << self.mrnm  
     data << self.mpos  
     data << self.isize 
     data << self.seq   
     data << self.qual 
     self.tags.each { |tag, x| 
       data << [x.tag, x.type, x.value].join(":")}
     
     data.map { |x| x.to_s }.join("\t")
  end

  def ==(x)
    return false if self.qname != x.qname
    return false if self.flag  != x.flag
    return false if self.rname != x.rname
    return false if self.pos   != x.pos 
    return false if self.mapq  != x.mapq
    return false if self.cigar != x.cigar
    return false if self.mrnm  != x.mrnm
    return false if self.mpos  != x.mpos
    return false if self.isize != x.isize
    return false if self.seq   != x.seq 
    return false if self.qual  != x.qual
    return false if self.tags  != x.tags
    return true
  end

  def cigar_hard_clip(cigar, x)
    hard = cigar[0][1] == "H" ? cigar[0][0] : 0
    hard += x
    left = x
    cigar = cigar.drop_while do |n, op|
      next false if left <= n
      left -= n
      next true
    end

    cigar[0][0] -= left
    cigar.unshift [hard, "H"]

    cigar
  end

  def remove_primer(region)
    result = Result.new

    result.aln = self.dup
    result.with_T = false
    result.cut = false
    result.drop = false

    return result if self.flag & 0x04 > 0 or region.chr != self.rname
    return result if self.pos < region.start - 100 or self.pos > region.start + 100

    # 10M3D5I => [[10, "M"], [3, "D"], [5, "I"]]
    cigar = self.cigar.scan(/(\d+[A-Z])/).map { |x| x[0].split(/(\d+)/)[1..-1]}.map { |a, b| [a.to_i,b]}
    calend = self.pos + cigar.find_all { |n, type| ["N", "M", "D"].include?(type) }.inject(0) { |sum, x| sum += x[0]}

    fixed_region_start = region.start
    fixed_region_finish = region.finish
    
    if self.pos == fixed_region_start - 1 then
      char = self.seq[0]
        
      if char == "T" then
        fixed_region_start -= 1
        result.with_T = true
      end

      if char == "A" then
        fixed_region_start -= 1
        result.with_T = true
      end 
    elsif calend - 1 == fixed_region_finish then
      char = self.seq[-1]
        if char == "A" then
          fixed_region_finish += 1
          result.with_T = true
        end
        if char == "T" then
          fixed_region_finish += 1
          result.with_T = true
        end 
    end

    if self.pos >= fixed_region_start and self.pos < fixed_region_finish then
      excess = fixed_region_finish - self.pos
      if excess >= self.seq.length then
        result.drop = true
        return result
      end
      result.aln.pos = fixed_region_finish
      result.aln.seq = self.seq[excess .. -1]
      result.aln.qual = self.qual[excess .. -1]
      result.aln.cigar = cigar_hard_clip(cigar, excess).map {|n, op| "#{n}#{op}"}.join
      result.cut = true
    end

    if calend > fixed_region_start and calend <= fixed_region_finish then
      excess = calend - fixed_region_start
      if excess >= self.seq.length then
        result.drop = true
        return result
      end
      result.aln.seq = self.seq[0 .. -excess - 1]
      result.aln.qual = self.qual[0 .. -excess - 1]
      result.aln.cigar = cigar_hard_clip(cigar.reverse, excess).reverse.map {|n, op| "#{n}#{op}"}.join
      result.cut = true
    end

    return result
  end

  def remove_primer_from_regions(regions)
    result = Result.new

    result.aln = self.dup
    result.with_T = false
    result.cut = false
    result.drop = false

    regions.each do |region|
      r = result.aln.remove_primer(region)

      result.aln = r.aln
      result.with_T ||= r.with_T
      result.cut ||= r.cut
      result.drop ||= r.drop

      break if result.drop == true
    end

    return result
  end
end
