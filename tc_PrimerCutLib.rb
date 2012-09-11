require_relative './PrimerCutLib.rb'
require 'test/unit'
require 'tempfile'

class TestRegion < Test::Unit::TestCase
  def test_init
    r = Region.new("1:1-100")
    assert_equal("1", r.chr)
    assert_equal(1, r.start)
    assert_equal(100, r.finish)

    r = Region.new("1\t1\t100")
    assert_equal("1", r.chr)
    assert_equal(1, r.start)
    assert_equal(100, r.finish)

    s = "1:100-XXX"
    assert_raise RuntimeError, "Failed to parse #{s.inspect} as a region."  do 
      Region.new(s)
    end

    s = "1 100\t1XXX"
    assert_raise RuntimeError, "Failed to parse #{s.inspect} as a region."  do 
      Region.new(s)
    end
  end
end

class TestRegionArray < Test::Unit::TestCase
  def test_load_regions_from_string
    regions = Array.load_regions("1:100-200,X:400-500")
    assert_equal(2, regions.size) 

    assert_equal("1", regions[0].chr) 
    assert_equal(100, regions[0].start) 
    assert_equal(200, regions[0].finish) 

    assert_equal("X", regions[1].chr) 
    assert_equal(400, regions[1].start) 
    assert_equal(500, regions[1].finish) 
  end

  def test_load_regions_from_file
    tmp = Tempfile.new 'bed'

    tmp.puts "1\t100\t200"
    tmp.puts "X 400 500"

    tmp.close

    regions = Array.load_regions(tmp.path)
    assert_equal(2, regions.size) 

    assert_equal("1", regions[0].chr) 
    assert_equal(100, regions[0].start) 
    assert_equal(200, regions[0].finish) 

    assert_equal("X", regions[1].chr) 
    assert_equal(400, regions[1].start) 
    assert_equal(500, regions[1].finish) 
  end
end

class TestAlignment < Test::Unit::TestCase
   
  def setup
    @aln = Bio::DB::Alignment.new
    @aln.sam = %{B7_610:8:68:570:705     99      seq2    910     99      
                 35M     =       1100    225     
                 GACAGATTTAAAAACATGAACTAACTATATGCTGG     
                 <<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8< 
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")
  end

  def test_not_overlapped_left_intact
    r = @aln.remove_primer(Region.new("1:1-100"))
    assert_equal(@aln, r.aln)
    assert_equal(false, r.drop)

    r = @aln.remove_primer(Region.new("seq1:910-915"))
    assert_equal(@aln, r.aln)
    assert_equal(false, r.drop)

    r = @aln.remove_primer(Region.new("seq2:1520-1525"))
    assert_equal(@aln, r.aln)
    assert_equal(false, r.drop)

    r = @aln.remove_primer(Region.new("seq2:900-910"))
    assert_equal(@aln, r.aln)
    assert_equal(false, r.drop)

    r = @aln.remove_primer(Region.new("seq2:945-950"))
    assert_equal(@aln, r.aln)
    assert_equal(false, r.drop)
  end

  def test_fully_inside_left_intact
    r = @aln.remove_primer(Region.new("seq2:915-920"))
    assert_equal(@aln, r.aln)
    assert_equal(false, r.drop)

    r = @aln.remove_primer(Region.new("seq2:911-920"))
    assert_equal(@aln, r.aln)
    assert_equal(false, r.drop)

    r = @aln.remove_primer(Region.new("seq2:930-944"))
    assert_equal(@aln, r.aln)
    assert_equal(false, r.drop)
  end

  def test_fully_inside_region_dropped
    r = @aln.remove_primer(Region.new("seq2:900-1000"))
    assert_equal(true, r.drop)

    r = @aln.remove_primer(Region.new("seq2:910-945"))
    assert_equal(true, r.drop)
  end

  def test_on_left_end_trimmed
    r = @aln.remove_primer(Region.new("seq2:910-915"))
    assert_not_equal(@aln, r.aln)
    assert_equal(r.with_T, false)
    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal(915, r.aln.pos)
    assert_equal("5H30M", r.aln.cigar)
    assert_equal("ATTTAAAAACATGAACTAACTATATGCTGG", r.aln.seq)
    assert_equal("<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<", r.aln.qual)
  end

  def test_on_left_end_with_A_extened_trimmed
    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    909     99      
                 36M     =       1100    225     
                 AGACAGATTTAAAAACATGAACTAACTATATGCTGG     
                 X<<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8< 
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")

    r = aln.remove_primer(Region.new("seq2:910-915"))
    assert_not_equal(aln, r.aln)
    assert_equal(r.with_T, true)
    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal(915, r.aln.pos)
    assert_equal("6H30M", r.aln.cigar)
    assert_equal("ATTTAAAAACATGAACTAACTATATGCTGG", r.aln.seq)
    assert_equal("<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<", r.aln.qual)
  end

  def test_on_left_end_with_T_extened_trimmed
    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    909     99      
                 36M     =       1100    225     
                 TGACAGATTTAAAAACATGAACTAACTATATGCTGG     
                 X<<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8< 
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")

    r = aln.remove_primer(Region.new("seq2:910-915"))
    assert_not_equal(aln, r.aln)
    assert_equal(r.with_T, true)
    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal(915, r.aln.pos)
    assert_equal("6H30M", r.aln.cigar)
    assert_equal("ATTTAAAAACATGAACTAACTATATGCTGG", r.aln.seq)
    assert_equal("<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<", r.aln.qual)
  end

  def test_on_right_end_trimmed
    r = @aln.remove_primer(Region.new("seq2:930-945"))
    assert_not_equal(@aln, r.aln)
    assert_equal(r.with_T, false)
    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal(910, r.aln.pos)
    assert_equal("20M15H", r.aln.cigar)
    assert_equal("GACAGATTTAAAAACATGAA", r.aln.seq)
    assert_equal("<<<<<<<<<<<<<<<<<<<<", r.aln.qual)
  end

  def test_on_right_end_with_A_extended_trimmed
    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    910     99      
                 36M     =       1100    225     
                 GACAGATTTAAAAACATGAACTAACTATATGCTGGA
                 <<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<X
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")

    r = aln.remove_primer(Region.new("seq2:930-945"))
    assert_not_equal(@aln, r.aln)
    assert_equal(r.with_T, true)
    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal(910, r.aln.pos)
    assert_equal("20M16H", r.aln.cigar)
    assert_equal("GACAGATTTAAAAACATGAA", r.aln.seq)
    assert_equal("<<<<<<<<<<<<<<<<<<<<", r.aln.qual)
  end

  def test_on_right_end_with_T_extended_trimmed
    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    910     99      
                 36M     =       1100    225     
                 GACAGATTTAAAAACATGAACTAACTATATGCTGGT
                 <<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<X
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")

    r = aln.remove_primer(Region.new("seq2:930-945"))
    assert_not_equal(@aln, r.aln)
    assert_equal(r.with_T, true)
    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal(910, r.aln.pos)
    assert_equal("20M16H", r.aln.cigar)
    assert_equal("GACAGATTTAAAAACATGAA", r.aln.seq)
    assert_equal("<<<<<<<<<<<<<<<<<<<<", r.aln.qual)
  end

  def test_unmapped_intact
    aln = Bio::DB::Alignment.new
    aln.sam = %{EAS54_65:7:56:57:985    117     seq2    1519    0     
                 *       =       1519    0       
                 TTCTGTCTTCTCTCCTGTCTTCTTTTCTCTTCTTT     
                 <9'<.<7<<2<<;77<7<<<<7<7<<<<7<<<2<<  
                 MF:i:192
                }.split.join("\t")

    assert_equal(aln, aln.remove_primer(Region.new("1:1-100")).aln)
    assert_equal(aln, aln.remove_primer(Region.new("seq2:1-100")).aln)
    assert_equal(aln, aln.remove_primer(Region.new("seq2:1520-1525")).aln)
  end

  def test_single_region_same_as_regions
    ["1:1-100", "seq1:910-915", "seq2:1520-1525", "seq2:900-910", "seq2:945-950", "seq2:915-920", "seq2:911-920",
     "seq2:930-944", "seq2:900-1000", "seq2:910-945", "seq2:910-915", "seq2:930-945"].each do |r|
      assert_equal(@aln.remove_primer_from_regions([Region.new(r)]), @aln.remove_primer(Region.new(r)), "Region = #{r.inspect}")
    end

    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    910     99      
                 36M     =       1100    225     
                 GACAGATTTAAAAACATGAACTAACTATATGCTGGA
                 <<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<X
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
          }.split.join("\t")

    r = "seq2:930-945" 
    assert_equal(aln.remove_primer_from_regions([Region.new(r)]), aln.remove_primer(Region.new(r)), "Region = #{r.inspect}")

    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    909     99      
                 36M     =       1100    225     
                 TGACAGATTTAAAAACATGAACTAACTATATGCTGG     
                 X<<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8< 
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")

    r = "seq2:910-915" 
    assert_equal(aln.remove_primer_from_regions([Region.new(r)]), aln.remove_primer(Region.new(r)), "Region = #{r.inspect}")
  end
     
  def test_single_region_same_as_regions_unmapped
    aln = Bio::DB::Alignment.new
    aln.sam = %{EAS54_65:7:56:57:985    117     seq2    1519    0     
                 *       =       1519    0       
                 TTCTGTCTTCTCTCCTGTCTTCTTTTCTCTTCTTT     
                 <9'<.<7<<2<<;77<7<<<<7<7<<<<7<<<2<<  
                 MF:i:192
                }.split.join("\t")

    assert_equal(aln.remove_primer_from_regions([Region.new("1:1-100")]), aln.remove_primer(Region.new("1:1-100")))
    assert_equal(aln.remove_primer_from_regions([Region.new("seq2:1-100")]), aln.remove_primer(Region.new("seq2:1-100")))
    assert_equal(aln.remove_primer_from_regions([Region.new("seq2:1520-1525")]), aln.remove_primer(Region.new("seq2:1520-1525")))
  end

  def test_several_overlapped_left_intact
    regions = ["1:1-100", "seq1:910-915", "seq2:100-200", "seq2:900-910", "seq2:945-950"].map { |x| Region.new(x) }
    r = @aln.remove_primer_from_regions(regions)
    assert_equal(@aln, r.aln)
    assert_equal(false, r.drop)
    assert_equal(false, r.with_T)
    assert_equal(false, r.cut)
  end

  def test_several_unmapped_intact
    aln = Bio::DB::Alignment.new
    aln.sam = %{EAS54_65:7:56:57:985    117     seq2    1519    0     
                 *       =       1519    0       
                 TTCTGTCTTCTCTCCTGTCTTCTTTTCTCTTCTTT     
                 <9'<.<7<<2<<;77<7<<<<7<7<<<<7<<<2<<  
                 MF:i:192
                }.split.join("\t")

    regions = ["1:1-100", "seq2:1-100", "seq2:1520-1525"].map { |x| Region.new(x) }
    r = aln.remove_primer_from_regions(regions)
    assert_equal(aln, r.aln)
    assert_equal(false, r.drop)
    assert_equal(false, r.with_T)
    assert_equal(false, r.cut)
  end

  def test_on_left_end_trimmed_regions
    regions = ["1:1-100", "seq2:1-100", "seq2:910-915"].map { |x| Region.new(x) }
    r = @aln.remove_primer_from_regions(regions)
    assert_not_equal(@aln, r.aln)
    assert_equal(r.with_T, false)
    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal(915, r.aln.pos)
    assert_equal("5H30M", r.aln.cigar)
    assert_equal("ATTTAAAAACATGAACTAACTATATGCTGG", r.aln.seq)
    assert_equal("<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<", r.aln.qual)
  end

  def test_on_left_end_with_A_extened_trimmed_regions
    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    909     99      
                 36M     =       1100    225     
                 AGACAGATTTAAAAACATGAACTAACTATATGCTGG     
                 X<<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8< 
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")

    regions = ["1:1-100", "seq2:1-100", "seq2:910-915"].map { |x| Region.new(x) }
    r = aln.remove_primer_from_regions(regions)
    assert_not_equal(aln, r.aln)
    assert_equal(r.with_T, true)
    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal(915, r.aln.pos)
    assert_equal("6H30M", r.aln.cigar)
    assert_equal("ATTTAAAAACATGAACTAACTATATGCTGG", r.aln.seq)
    assert_equal("<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<", r.aln.qual)
  end

  def test_on_left_end_with_G_extened_trimmed_regions
    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    909     99      
                 36M     =       1100    225     
                 GGACAGATTTAAAAACATGAACTAACTATATGCTGG     
                 X<<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8< 
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")

    regions = ["1:1-100", "seq2:1-100", "seq2:910-915"].map { |x| Region.new(x) }

    r = aln.remove_primer_from_regions(regions)
    assert_equal(aln, r.aln)
    assert_equal(r.with_T, false)
    assert_equal(r.drop, false)
    assert_equal(r.cut, false)
  end

  def test_on_left_end_with_GG_extened_trimmed_regions
    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    909     99      
                 37M     =       1100    225     
                 GGGACAGATTTAAAAACATGAACTAACTATATGCTGG     
                 XX<<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8< 
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")

    regions = ["1:1-100", "seq2:1-100", "seq2:910-915"].map { |x| Region.new(x) }

    r = aln.remove_primer_from_regions(regions)
    assert_equal(aln, r.aln)
    assert_equal(r.with_T, false)
    assert_equal(r.drop, false)
    assert_equal(r.cut, false)
  end

  def test_on_left_end_with_T_extened_trimmed_regions
    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    909     99      
                 36M     =       1100    225     
                 TGACAGATTTAAAAACATGAACTAACTATATGCTGG     
                 X<<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8< 
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")

    regions = ["1:1-100", "seq2:1-100", "seq2:910-915"].map { |x| Region.new(x) }
    r = aln.remove_primer_from_regions(regions)
    assert_not_equal(aln, r.aln)
    assert_equal(r.with_T, true)
    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal(915, r.aln.pos)
    assert_equal("6H30M", r.aln.cigar)
    assert_equal("ATTTAAAAACATGAACTAACTATATGCTGG", r.aln.seq)
    assert_equal("<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<", r.aln.qual)
  end

  def test_on_right_end_trimmed_regions
    regions = ["1:1-100", "seq2:1-100", "seq2:930-945"].map { |x| Region.new(x) }
    r = @aln.remove_primer_from_regions(regions)
    assert_not_equal(@aln, r.aln)
    assert_equal(r.with_T, false)
    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal(910, r.aln.pos)
    assert_equal("20M15H", r.aln.cigar)
    assert_equal("GACAGATTTAAAAACATGAA", r.aln.seq)
    assert_equal("<<<<<<<<<<<<<<<<<<<<", r.aln.qual)
  end

  def test_on_right_end_with_A_extended_trimmed_regions
    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    910     99      
                 36M     =       1100    225     
                 GACAGATTTAAAAACATGAACTAACTATATGCTGGA
                 <<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<X
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")

    regions = ["1:1-100", "seq2:1-100", "seq2:930-945"].map { |x| Region.new(x) }
    r = aln.remove_primer_from_regions(regions)
    assert_not_equal(aln, r.aln)
    assert_equal(r.with_T, true)
    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal(910, r.aln.pos)
    assert_equal("20M16H", r.aln.cigar)
    assert_equal("GACAGATTTAAAAACATGAA", r.aln.seq)
    assert_equal("<<<<<<<<<<<<<<<<<<<<", r.aln.qual)
  end

  def test_on_right_end_with_C_extended_trimmed_regions
    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    910     99      
                 36M     =       1100    225     
                 GACAGATTTAAAAACATGAACTAACTATATGCTGGC
                 <<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<X
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")

    regions = ["1:1-100", "seq2:1-100", "seq2:930-945"].map { |x| Region.new(x) }
    r = aln.remove_primer_from_regions(regions)
    assert_equal(aln, r.aln)
    assert_equal(r.with_T, false)
    assert_equal(r.drop, false)
    assert_equal(r.cut, false)
  end

  def test_on_right_end_with_CC_extended_trimmed_regions
    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    910     99      
                 37M     =       1100    225     
                 GACAGATTTAAAAACATGAACTAACTATATGCTGGCC
                 <<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<XX
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")

    regions = ["1:1-100", "seq2:1-100", "seq2:930-945"].map { |x| Region.new(x) }
    r = aln.remove_primer_from_regions(regions)
    assert_equal(aln, r.aln)
    assert_equal(r.with_T, false)
    assert_equal(r.drop, false)
    assert_equal(r.cut, false)
  end

  def test_on_left_end_trimmed_regions_two_from_begin
    regions = ["1:1-100", "seq2:915-920", "seq2:910-915"].map { |x| Region.new(x) }
  end

  def test_on_right_end_trimmed_two
    regions = ["seq1:910-915", "seq2:930-945", "seq2:925-930"].map { |x| Region.new(x) }
  end

  def test_on_left_end_trimmed_regions_two_from_begin_and_end
    regions = ["1:1-100", "seq2:910-915", "seq2:930-945"].map { |x| Region.new(x) }
    r = nil
    assert_nothing_thrown() { r = @aln.remove_primer_from_regions(regions) }
    assert_not_equal(@aln, r.aln)
    assert_equal(r.with_T, false)
    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal(915, r.aln.pos)
    assert_equal("5H15M15H", r.aln.cigar)
    assert_equal("ATTTAAAAACATGAA", r.aln.seq)
    assert_equal("<<<<<<<<<<<<<<<", r.aln.qual)
  end

  def test_on_left_end_trimmed_regions_two_from_begin_and_end_extened
    ["A", "T"].each do |ex|
    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    909     99      
                 36M     =       1100    225     
                 #{ex}GACAGATTTAAAAACATGAACTAACTATATGCTGG     
                 X<<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8< 
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")
    regions = ["1:1-100", "seq2:910-915", "seq2:930-945"].map { |x| Region.new(x) }
    r = nil
    assert_nothing_thrown(ex) { r = aln.remove_primer_from_regions(regions) }
    assert_not_equal(aln, r.aln, ex)
    assert_equal(r.with_T, true, ex)
    assert_equal(r.drop, false, ex)
    assert_equal(r.cut, true, ex)
    assert_equal(915, r.aln.pos, ex)
    assert_equal("6H15M15H", r.aln.cigar, ex)
    assert_equal("ATTTAAAAACATGAA", r.aln.seq, ex)
    assert_equal("<<<<<<<<<<<<<<<", r.aln.qual, ex)
    end
  end

  def test_on_right_end_trimmed_regions_two_from_begin_and_end_extened
    ["A", "T"].each do |ex|
    aln = Bio::DB::Alignment.new
    aln.sam = %{B7_610:8:68:570:705     99      seq2    910     99      
                 36M     =       1100    225     
                 GACAGATTTAAAAACATGAACTAACTATATGCTGG#{ex}
                 <<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<X
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")
    regions = ["1:1-100", "seq2:910-915", "seq2:930-945"].map { |x| Region.new(x) }
    r = nil
    assert_nothing_thrown(ex) { r = aln.remove_primer_from_regions(regions) }
    assert_not_equal(aln, r.aln, ex)
    assert_equal(r.with_T, true, ex)
    assert_equal(r.drop, false, ex)
    assert_equal(r.cut, true, ex)
    assert_equal(915, r.aln.pos, ex)
    assert_equal("5H15M16H", r.aln.cigar, ex)
    assert_equal("ATTTAAAAACATGAA", r.aln.seq, ex)
    assert_equal("<<<<<<<<<<<<<<<", r.aln.qual, ex)
    end
  end
end

class TestAlignmentBWA < Test::Unit::TestCase
  def setup
    @aln = Bio::DB::Alignment.new
    @aln.sam = %{B7_610:8:68:570:705     99      seq2    910     99      
                 3I3D29M     =       1100    225     
                 GACTTTAAAAACATGAACTAACTATATGCTGG     
                 <<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<8< 
                 MF:i:18 Aq:i:30 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
                }.split.join("\t")
  end

  def test_on_left_end_trimmed
    r = @aln.remove_primer(Region.new("seq2:910-915"))
    assert_not_equal(@aln, r.aln)
    assert_equal(r.with_T, false)
    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal(915, r.aln.pos)
    #assert_equal("5H1D29M", r.aln.cigar)
    #assert_equal("ATTTAAAAACATGAACTAACTATATGCTGG", r.aln.seq)
    assert_equal("<<<<<<<<<<<<<<<<<<<<<;<<<<<<8<", r.aln.qual)
  end
end

class TestRealWorld < Test::Unit::TestCase
  def test_one
    aln = Bio::DB::Alignment.new
    aln.sam = %{1031_1902_651	16	1	150238458	5	8H27M15H	
                *	0	0	NTNNNNNTGGAGAAATACAGGGCGAGG	
                *-)&%%'/66DFD@DJJJJJJJJJJJ5	RG:Z:test	NH:i:1	
                CM:i:2	NM:i:0	CQ:Z:1%B674:+@>3<BA)86BB62:/:;5?6'4/4-B/-)8/))'-)/37)0,	
                CS:Z:T20213001022122110223300211330022201222101221223102}.split.join("\t")
    
    r = aln.remove_primer(Region.new("1:150238466-150238485"))

    assert_equal(r.drop, false)
    assert_equal(r.cut, true)
    assert_equal("8H8M34H", r.aln.cigar)
    assert_equal("NTNNNNNT", r.aln.seq)
  end
end
