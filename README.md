# CreateFragmentMatrix
Converts input bam and vcf files to a fragment matrix

### Dependencies
*  [Pysam v.0.13](http://pysam.readthedocs.io/en/latest/#)
### Usage:
```
python3 FragMatrixCreator /path/to/bam /path/to/vcf genomic_region /output/folder
```
### Output:
Fragment matrix /output/folder/genomic_region.frags

### Optional arguments
```
--output_prefix PREFIX      add prefix to output file: /output/folder/prefix_genomic_region.frags
--genotypes                 output additional /output/folder/genomic_region.genotypes file
--se                        set this flag if the reads are single-end
```
### Output format:
#### Fragment matrix:
A fragment matrix always starts with with `>` sign and a genomic region name after it:  
```>scaffold2314|size104419```  
Next lines are fragments, their representation has 5 tab delimited columns:  
```NS500442:8:H194MBGXX:3:13604:2918:8976   352 00-00   F)-F@   55```
1. *Read name*
2. *Start position* (index of the corresponding polymorphic site in input vcf, 1-based)
3. *Fragment sequence*, `0` is REF, `1` is first ALT, 2 is second, etc., whereas `-` means no information available
4. *Base call qualities* in Phred+33 format
  * In case of complex polymorphic sites with multiple bases the quality is `mean(bases)`;
  * The length of the quality string is equal to the fragment string, therefore if the fragment string contains missing alleles (`-`), the quality string will also have `-` on same positions. The dashes on other positions in quality string denote base quality 12 (`ord('-') - 33`).
5. *Mapping quality* 


The fragments are separated by new lines.  
**Note**: A fragment with only one allele may appear in a fragment matrix only if it has a mate with at least one *different* allele in the same genomic region  
Some examples are shown in `fragment_examples.txt`.
#### Genotypes:
(Output only if `--genotypes` flag is set)  
Start with genomic region name:  
```>scaffold2314|size104419```  
Following with genotypes information from vcf file for every polymorphic sites:  
```1    0/0/1/1/1/1:96:38:1272:57:1673```
1. *Polymorphic site index* (starting with 1 in every new genomic region)
2. *GT column* from vcf file


The genotypes are separated by new lines.
