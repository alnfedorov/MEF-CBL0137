In total, 5 ChIP-seq experiments were analyzed with the following antibodies: **Z22**(curaxin 14h), **IgG**(curaxin 14h/0h), and **FLAG**(curaxin 14h/0h). Each experiment had two technical replicates.

Here, curaxin 0h refers to untreated cells, **Z22** is an antibody for Z-DNA, **IgG** is a nonspecific antibody (control), and **FLAG** is an antibody to FLAG-
tagged ZBP1 protein. 

### Pipeline
Initial processing was performed via the Snakemake pipeline present in the `pipeline` folder:
```bash
cd pipeline
sudo docker build -t chipseq:latest docker/
sudo docker run --rm -it -v $(pwd):/project --name chipseq-container chipseq:latest
snakemake --cores $(nproc) --configfile config/config.yaml -d . --snakefile workflow/Snakefile  all
```
Before running the pipeline, make sure to place raw fq.gz files in `pipeline/results/pe-{SAMPLE NAME}/fq.gz` folders. Refer to the `pipeline/config/config.yaml` for a list of samples used in the research (fq.gz files can be downloaded from the SRA).

Once the pipeline finishes, called peaks, fold enrichment tracks and basic QC metrics will be available in the `pipeline/results` folder.

### Artifacts
#### Resources
The `artifacts` folder contains scripts to reproduce chip-seq relevant plots and supplements. These scripts depends on the pipeline results, which can be linked like this:
```bash
# genome
sudo ln -f $(pwd)/pipeline/resources/mm10/genome.fa $(pwd)/artifacts/resources/genome.fa
# fold enrichment
cp -R pipeline/results/signal/ artifacts/resources/signal
# peaks
cp -R pipeline/results/peaks/ artifacts/resources/peaks
```
There are additional external data (RepeatMasker annotation, ENCODE blacklist, L1Base) stored directly in the repository. To avoid running the entire pipeline, one can use signal and peak files from the GEO submission.

Some artifacts require a compiled zhunt program, which can be built with GCC: `gcc artifacts/resources/zhunt2.c -lm -o artifacts/resources/zhunt2`.
#### Reproduction of results 
The first step to reproduce artifacts is to restore the environment described in the docker file:
```bash
cd artifacts
sudo docker build -t chipseq-artifacts:latest docker/
sudo docker run --rm -it -v $(pwd):/project --name chipseq-artifacts-container chipseq-artifacts:latest
```
Then run any **Python** script in the `scripts` folder, for example: `python3 scripts/mm10-coverage-by-repeats`. 

Note that figures in the publication are hand-edited versions of graphs created with these scripts. 
