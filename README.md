# RNA-seq pipeline

__Author__ : Danny Bergeron

__Email__ :  _<danny.bergeron@usherbrooke.ca>_
## Software to install
Conda (Miniconda3) needs to be installed (https://docs.conda.io/en/latest/miniconda.html)

For Linux users :
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Answer `yes` to `Do you wish the installer to initialize Miniconda3?`


To create the Snakemake environment used to launch Snakemake, run the following. The `conda create` command can appear to be stuck on `Solving environment`. While we are actually arguably [never going to solve the environment](https://www.ipcc.ch/sr15/chapter/spm/), the command is probably not stuck. Just be patient.

```bash
exec bash
conda config --set auto_activate_base False
conda create --name smake -c bioconda -c conda-forge snakemake=5.4.5
```

Before running Snakemake, you have to initialize the environment
```bash
conda activate smake
```


If working on a cluster, either go for a local installation, or check if it is not already installed on your system.


## Run
To run the workflow locally simply run the following command in the Snakemake conda environment, where `$CORES` is the number of available cores.
```bash
snakemake --use-conda --cores=$CORES
```

To run on a Slurm cluster, one can use the following command to output all tasks at once.
```bash
snakemake -j 999 --use-conda --immediate-submit --notemp --cluster-config cluster.json --cluster 'python3 slurmSubmit.py {dependencies}'
```

If the cluster nodes do not have internet access, one can run the tasks requiring the internet locally with :
```bash
snakemake all_downloads --use-conda --cores=$CORES
```

To look at your entire workflow in svg, use the following command :
```bash
snakemake --rulegraph | dot -Tsvg | display
``` 
