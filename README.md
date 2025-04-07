<a name="top"></a>
# Ghostbuster
1. [Introduction](#intro)
2. [Installing Dependencies](#dependencies)
3. [Running Ghost Buster](#running)
4. [Arguments](#arguments)
5. [Examples](#examples)

## <ins>**Introduction**</ins> <a name="intro"></a>

Four-taxon tests (e.g. Patterson's D-statistic) are widely used to identify cases of introgression between non-sister lineages in the analysis. However, there is growing evidence that introgression from an unsampled/extinct "ghost lineage" can obscure interpretation of four-taxon tests. Ghostbuster is a computational tool that applies a novel divergence-based test to distinguish between ghost lineage introgression and ingroup introgression. Ghostbuster is intended for cases in which a prior four-taxon test indicated significant introgression.

A manuscript describing the method is currently under review. Check back here for citation information.

## <ins>**Installing Dependencies**</ins> <a name="dependencies"></a>

We recommmend using conda to create a virtual environment with the needed programs installed. We have created an environment with all the needed despendencies, which we exported in the file, ghostbuster_evn.yml. You can create this environment on your system using the following command:

```bash
conda env create -f ghostbuster_env.yml
```

## <ins>**Running Ghostbuster**</ins> <a name="running"></a>

## <ins>**Arguments**</ins> <a name="arguments"></a>
These arguments can also be referenced in your command line using the -h flag.
```bash
./Ghost_Buster.py -h
```

| Argument             | Description                                                                                                                    | Required |
|----------------------|--------------------------------------------------------------------------------------------------------------------------------|----------|
| `-h`, `--help`        | Show this help message and exit                                                                                                | No       |
| `-out`, `--outgroup`  | Unique string found in outgroup seqID                                                                                          | Yes       |
| `-P1`, `--pop1`       | Unique string found in pop1 seqID                                                                                               | Yes       |
| `-P2`, `--pop2`       | Unique string found in pop2 seqID                                                                                               | Yes       |
| `-P3`, `--pop3`       | Unique string found in pop3 seqID                                                                                               | Yes       |
| `-i`, `--indir`       | Full path to input directory (should end in `/`)                                                                               | Yes       |
| `-j`, `--job`         | Name for the job to create file names                                                                                          | Yes       |
| `-t`, `--threads`     | Number of threads to use for tree inference                                                                                     | Yes       |
| `-s`, `--skip`        | Add this flag to skip the gene tree inference step. Requires an `OUT_*` directory with a `node_depths.csv` file already present | No       |

## <ins>**Examples**</ins> <a name="examples"></a>

Here is an example command used to run ghostbuster:

```bash
python Ghost_buster.py -i Examp_files/ -j test -t 2 -P1 Bs_ -P2 Cr_ -P3 At_ -out Es_
```

This run will create a folder named "OUT_test/" and write results to the folder.


[Back to Top](#top)





