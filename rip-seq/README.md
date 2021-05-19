### Pipeline
Preprocessing, alignment, and quality control of RIP-seq experiments were performed via `nf-core/rna-seq` pipeline using the following commands:
```bash
cd pipeline
export NXF_VER=20.11.0-edge
curl -s https://get.nextflow.io | bash
sudo ./nextflow -bg run nf-core/rnaseq -resume -profile docker -revision 3.0 --input=design.csv --genome=GRCm38
```
To start a pipeline from scratch, one needs to download raw fq.gz reads from the SRA, place them in the `pipeline/fq.gz` folder, and modify the design document accordingly.

We also additionally filtered the resulting BAM files with the following command:
```bash
cd pipeline/results/star_salmon
mkdir final-bams
for file in *.bam
do
    samtools view -@ $(nproc) -F 2820 -b $file > final-bams/$file
done
```
### Analysis
All software dependencies are described in a Dockerfile, which can be built like this:
```bash
cd analysis
sudo docker build -t ripseq-analysis:latest docker/
# run the container
sudo docker run --rm -it -v $(pwd):/project --name ripseq-analysis-container ripseq-analysis:latest
```
Then run any **Python** script in the `scripts` folder. Additional notes for some steps are listed below.
#### DE genes
We used DESeq2 to determine genes enriched in Z22 over IgG. To run the DE analysis, first move Salmon per-gene quantification results `pipeline/results/star_salmon/SAMPLE_R1/quant.genes.sf` to the `analysis/resources/Salmon/SAMPLE.quant.genes.sf`. Then amend the `Salmon/samples.tsv` file if you used a different design for the pipeline and run the Python script: `scripts/run-DESeq2.py`.
#### Alu editing index (AEI)
To compute the AEI index we used the following commands:
```bash
# Go to directory containing final-bams folder with filtered BAM files
cd pipeline/results/star_salmon/final-bams
# Reheader each BAM file (1->chr1)
for file in *.bam
do
    newfile=${file/.bam/_chr.bam}
    samtools view -H $file | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - $file > $newfile
    rm -rf $file
    mv $newfile $file
done
cd ...
# build docker image
sudo docker build -t rnaedits https://raw.githubusercontent.com/a2iEditing/RNAEditingIndexer/master/Dockerfile
# run the container
sudo docker run -v $(pwd):/data --rm -it ranedits bash
# launch the pipeline
RNAEditingIndex -d /data/final-bams -l /data/logs -o /data/output -os /data/editing-summary -f _chr.bam --genome mm10 --ts $(nproc)
```
The output file located in the `summary` folder was hand-formatted and then attached as a supplement table to the paper.
#### Interferome database and IFN-I stimulated genes (ISG)
To parse the interferome database, we divided all mm10 genes (GENCODE vM25 annotation) into 6 batches and manually submitted each of them to the interferome web interface. The server output for each package is located in the `results/Interferome` folder.

Run the Python script `scripts/infer-ISG-interferome.py` to merge them together and create a summary table with ISG showing the median median fold change > 1.5.

Additionaly, run `scripts/infer-ISG-MEF-NGS.py` to determine ISG in MEFs based on the NGS data. The script will: (1) download gene-counts for untreated and IFNb stimulated genes from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128110); (2) run DESeq2 to detect ISG; (3) save results.
#### Mapping fragments to exons/introns/intergenic regions
To calculate a summary table showing the percentage of fragments falling into exons / introns / intergenic regions, it is necessary to sort and place the BAM files in the `resources / BAM` folder:
```bash
cd pipeline/results/star_salmon/final-bams
for file in *.bam
do
    newfile=../../../../analysis/resources/BAM/${file/_R1.markdup.sorted.bam/.namesorted.bam}
    samtools sort -@ $(nproc) -n -o $newfile $file
done
```
Then one can run the `map-to-exons-introns-intergenic.py` script.
#### Coverage tracks
Simple coverage tracks for Eif2ak2/Dddx58/Ifih1 can be built via `coverage-tracks-Z22-genes.py` script. The script depends on depth normalized coverage tracks, which can be built using deeptools:
```bash
# create blacklist with rRNA/tRNA genes and repeats
python3 scripts/make-coverage-blacklist.py
# ..... -> reheader bam files as shown in the AEI section
cd pipeline/results/star_salmon/final-bams
blacklist=../../../../analysis/resources/coverage-blacklist.bed.gz
for file in *.bam
do
    bw=${file/.bam/.bw}
    bamCoverage -bl $blacklist --normalizeUsing CPM -p $(nproc) -b $file -o $bw
done
mv *.bw ../../../../analysis/resources/signal/
```
