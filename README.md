# DNA read simulation and bridging

## Configurations
you can change settings in **config.json** only the config file in the specific directory would determin the settings for the main python file 
### Definitions
* Template size : the size of the fuly constructed read
* SNP Count : the amount of SNPs in the read that are generated
* Average Depth : the amount of reads generated
* Read Size : size of pair end read
* Insert Size : read size + insert hole + read size = Insert Size
* Start , End : the amount to attempt from start to end

## Single simulation
**Basefile.py** it simulates a template, type a, type b not based on any actual information and attemps to blend then reconstruct itself

## Barchart Comparison
**Barchart.py** like single simulation but attempts to vary the SNP count and makes a chart for determining the amount of groups the output get generated

## Resplice FastA
**FastASplicer.py** uses *FastAFile.fasta* as base for information instead of randomly generating fake info

## FastQ Simulator
**FastQMaker.py** uses *FastABase.fasta* as base and creates mock fastq file for continuation in BWA for SAM files to continue on next header, outputs *output.fastq*

## Sam Bridger
**SamSplicer.py** uses *SamFile.sam* as information of the reads and the position on reference file *ReferenceFile.fasta* to bridge the reads, this is the latest version and should be used for all real world information
