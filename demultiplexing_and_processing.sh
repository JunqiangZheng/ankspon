# JB's version of JM's script:
# Update Jan 30, 2013
# Ion Torrent now outputs a bam file. Convert bam to sam and use new perl script from GG to get tabs
#	file (see changes below)
#---------------------------------------------------------------------------------------


########################################################################################
#Set the absolute directories that you will be running this script from
ROOT="/Users/jbisanz/Documents/Research/as_plaques"
BIN="$ROOT/bin"
DATA_DIR="$ROOT/data"
OUT_DIR="$ROOT/workflow_out"

########################################################################################
#Set the absolute directories that you will be running this script from
BAM_FILE="$DATA_DIR/as_plaques_2PG-98.bam"
SAM_FILE="$DATA_DIR/as_plaques_2PG-98.sam" #this is what the converted file will be named
KEY_FILE="$DATA_DIR/as_key.txt"

########################################################################################
#Set paramaters for analysis
MINIMUM_READS="500" #number of reads needed for a barcode (keeps barcode sequencing errors out)
CLUSTER="0.97" #minimum identity to cluster for OTUS, 97% is standard
REM_CUTOFF="1" #the percentage at which OTUs are clustered into the remainder for tables



########################################################################################

echo "Executing Ion Torrent V6 Pipeline, be sure to read header of file and set parameters and file locations"
echo ""

#Make the directories if they don't already exist, data_dir should already exist with your bam and key file
#mkdir $DATA_DIR
mkdir $OUT_DIR

#Convert to SAM if it doesn't already exist.
echo "Converting BAM File TO SAM File"
if [ -e $SAM_FILE ]; then
	echo "-BAM ALREADY CONVERTED TO SAM"
else
	#convert BAM to SAM
	$BIN/samtools view $BAM_FILE > $SAM_FILE
fi
echo "...done"

#extract the reads
#Will also do a size-selection: keeping reads that are >= 70 to 90bp between the L and R primers
#ARGV[0] is the SAM file
#ARGV[1] is the X value: barcodes with over X reads are kept (default: 500)
#output to tabbed_reads.txt
echo "Extraction Reads from SAM file with a minimum number of reads per barcode of $MINIMUM_READS"
if [ -e $OUT_DIR/tabbed_reads.txt ]; then
	echo "-Reads already extracted"
else
	$BIN/extract_torrent_sam_to_tabs.pl $SAM_FILE $MINIMUM_READS
	mv $DATA_DIR/tabbed_reads.txt $OUT_DIR/
	mv $DATA_DIR/barcode_counts.txt $OUT_DIR/
fi
echo "...done"

echo "Changing barcodes to sample IDs"
#change the barcodes to sample ids
$BIN/rekey_merged.pl $KEY_FILE $OUT_DIR/tabbed_reads.txt > $OUT_DIR/keyed_tabbed_reads.txt
echo "...done"

#generates a file of reads that are grouped by identity (ISU) as single line format saved into data/groups.txt
#also generates an index file that correlates the ISU sequences back to each sample, not interpretable 
echo "Grouping into ISUs"
$BIN/group.pl $OUT_DIR/keyed_tabbed_reads.txt
echo "...done"

#mv data/reads_in_groups.txt workflow_out/
#mv data/groups.txt workflow_out/

echo "Making FASTA File for ISUs (groups_all.fa)"
#convert to fasta format
awk '{print$1 "\n"  $2}' $OUT_DIR/groups.txt > $OUT_DIR/groups_all.fa
echo "...done"

echo "Clustering OTUs at $CLUSTER %"
$BIN/uclust3.0.617_i86darwin32 --input $OUT_DIR/groups_all.fa --uc $OUT_DIR/all_all.uc --id $CLUSTER --usersort 2> $OUT_DIR/uc.out
echo "...done"

echo "Generating read table with ISU and OTUs"
#re-generate the data file with ISU and OTU groupings added to the read identifier
$BIN/map_otu_isu_read.pl $OUT_DIR/all_all.uc $OUT_DIR/reads_in_groups.txt $OUT_DIR/keyed_tabbed_reads.txt > $OUT_DIR/mapped_otu_isu_reads.txt
echo "...done"

echo "Making data tables and creating fastas for ISU and OTUs, minimum cut off for OTU table is $REM_CUTOFF %"
#make the table of values for each sample with a percent cutoff, use 0.5 for small numbers of samples
#use 1 to 2 for large numbers
#expects a directory called analysis

mkdir analysis

$BIN/get_tag_pair_counts.pl $OUT_DIR/mapped_otu_isu_reads.txt $REM_CUTOFF > $OUT_DIR/filter.out

#get your OTUs and their counts as a fasta file
$BIN/get_seed_otus.pl $OUT_DIR/all_all.uc $OUT_DIR/groups_all.fa analysis/OTU_tag_mapped.txt > $OUT_DIR/all_seed_OTUs.fa

#must rotate the table 90 degrees
# This puts OTU as column labels and sample IDs as row headers
cd analysis
#OTU table
R CMD BATCH $BIN/OTU_to_QIIME.R OTU_to_QIIME.out
#mv OTU_to_QIIME.out analysis/OTU_to_QIIME.out
#output is to td_OTU_tag_mapped.txt
cd ../

echo "Complete!"
echo "Run seqmatch from RDP and run taxonomy_qiime.sh to get phylogenies and files for QIIME input"
echo "Use the following settings: #Strain=BOTH Source=ISOLATES size=BOTH Quality=GOOD Taxonomy=Nomenclatural KNN Matches=20"
echo "Put seqmatch_download.txt in a folder called taxonomy"

exit
#-----------
#	RDP taxonomy
#-------------
TAX_DIR="taxonomy"

#parse the RDP output 
$BIN/parse_RDP.pl $TAX_DIR/seqmatch_download.txt > $TAX_DIR/parsed_RDP.txt

$BIN/RDP_lineage.pl $TAX_DIR/parsed_RDP.txt analysis/td_OTU_tag_mapped.txt > analysis/td_OTU_tag_mapped_RDPlineage.txt

cd analysis/
mkdir qiime

macqiime
cd qiime

#convert to new .biom format
convert_biom.py -i ../td_OTU_tag_mapped_RDPlineage.txt -o td_OTU_RDPlineage_fixed.biom --biom_table_type="otu table" --biom_type=dense --process_obs_metadata taxonomy

