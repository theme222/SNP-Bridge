# DNA read simulation and bridging

## Configurations
you can change settings in **config.json** only the config file in the specific directory would determine the settings for the main python file 
### Definitions
* Template size : the size of the fuly constructed read
* SNP Count : the amount of SNPs in the read that are generated
* Average Depth : the amount of reads generated
* Read Size : size of pair end read
* Insert Size : read size + insert hole + read size = Insert Size
* Start , End : the amount to attempt from start to end

## Directories in this project

### Single simulation
**Basefile.py** it simulates a template, type a, type b not based on any actual information and attempts to blend then reconstruct itself this is the first prototype of this project and can only simulate SNPs and doesn't handle any type of insertion deletion or broken information

### Barchart Comparison
**Barchart.py** like single simulation but attempts to vary the SNP count and makes a chart for determining the amount of groups the that the output generates

### Resplice FastA
**FastASplicer.py** uses *FastAFile.fasta* as base for information instead of randomly generating fake info and proceeds to use a similar technique in *basefile.py*

### FastQ Simulator
**FastQMaker.py** uses *FastABase.fasta* as base and creates mock fastq file for continuation in BWA for SAM files to continue on outputs *output.fastq* 

### Sam Bridger
**SamSplicer.py** uses *SamFile.sam* as information of the reads and the position on reference file *ReferenceFile.fasta* to bridge the reads, this is the latest version and should be used for all real world information

### PhasePipeline 
**Start.py** using files set in *config.json* to change filetypes through bwa and outputs the final result based on the script in *SamSplicer.py* use this if you are looking for direct from the sequencer information read *README.md* in the directory for more info 

### PhasePipelineMRef 
MRef stands for Multi Reference which basically means it is able to use multiple references which are stored in the References subdirecrtory

## How the code works
I will be describing the inner workings of the Phase Pipeline as it is the most relevant one out of all the projects listed here, and it is a combination of almost all the other directories
1. After receiving a fastQ file or generating it in FastQ Simulator and getting relevant reference info the code will run the info through bwa 
   * the script will run shell commands through `os.system`
   * it will index the reference with the command `bwa index` 
   * the script will run `bwa mem` to align the fastq file to the reference
   * BWA will produce a *bwa-output.sam* file which will get deleted in the future
2. The script will use `pysam` to read the newly created sam file
3. Using the reference file it will detect if the sam read has a "CIGAR string" and will clip away all information that bwa deems to be hard clips or soft clips as they are not usually reliable info (example: 5S95M, 70M30H)
4. Since the CIGAR string mentioned earlier can't specify SNPs as it will identify them as a match (example: 100M) we will use `Bio.Align.PairwiseAligner` to do the aligning for us and then gather the positions of the snp relative to the reference and relative to the read which will then be put into a class object called `SNP`
5. The reads collected from the sam are then translated into a class object called `DNA`
6. We will need to delete some SNPs from the list generated because of these problems
   1. When generating the list of SNPs it will go through every single read and check the differences as mentioned in step 4 but this will give duplicate SNPs that the script needs to remove
   2. It currently contains positions of SNPs of the reads that are both different when compared to the reference but when compared to each other it has the exact same base (example: ref "ATTACG" read 1: "ATGACG" read 2: "ATGACG") in position 3 the program will need to remove the SNP in that area
   3. Due to either an error with BWA or an error with `Bio.Align.PairwiseAligner` or even faulty information the script will look through all possible alleles in that position to place into the `SNP` class as `SNP.values` and collect the two most common types and in the future the code will discard any read that contains other alleles as useless information
7. It then will start "Bridging"
   1. First, it will find all possible reads that overlap a SNP
   2. It will then set the first read it detects and adds a "pivot" to be compared to afterwords 
   3. Using said pivot it will separate the reads into two groups
   4. It will then grab the two groups and put it into a bigger group which is determined by if the bigger group contains a read from the first group
   5. It will then repeat steps 1-4 until it has completed all the SNP locations 
   6. It will then clean up the information to only include unique reads using `drill_down()`
8. Finally, it displays the information using `matplotlib.pyplot` and writes to a FastA file as output
9. It will also delete any extra files generated by BWA which are just the reference index files (.amb, .ann, .bwt, .sa, .pac) 
