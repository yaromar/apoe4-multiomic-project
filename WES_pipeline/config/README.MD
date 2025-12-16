These files contain the BED files for the IDT xGen exome kit, along with the YCGA-designed spikein 
to improve RefGene coding region coverage.  Each of the bed files is a sorted, merged bed file (meaning
the regions are sorted in chromosome order, chr1, chr2, ..., chr22, chrX, chrY, chrM, and overlapping
regions are merged together, so that the bed file contains only non-overlapping regions).

   * hg19_idt_target_regions.bed - The target regions of the IDT xGen exome kit
   * hg19_spikein_regions.bed - The target regions for the YCGA-designed spikein
   * hg19_idt_and_spikein_target_regions.bed - The combined bed file, with all regions

The IDT xGen exome kit was designed using hg19/GRCh37, and so no official hg38/GRCh38 regions are
provided.  These files are crossmapped files using the hg19 bed files, and retain all but ~11kb of
the original size.

   * hg38_idt_target_regions.bed - The crossmapped IDT xGen target regions
   * hg38_spikein_regions.bed - The crossmapped YCGA-designed spikein regions
   * hg38_idt_and_spikein_target_regions.bed - The combined bed file
