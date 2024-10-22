<a name="top"></a>
# Ghost Buster
1. [Introduction](#intro)
2. [Installing Dependencies](#dependencies)
3. [Running Ghost Buster](#running)
4. [Arguments](#arguments)
5. [Examples](#examples)

## <ins>**Introduction**</ins> <a name="intro"></a>
Ghost Buster is a python tool for testing whether introgression events have the divergence signatures most consistent with 'true'/'ingroup' introgression versus 'ghost' introgression. Note that Ghost Buster is not meant to identify whether introgression has occured in the first place (see existing tests such as the D-statistic and much more). Instead, Ghost Buster should only be applied to gain further resolution/insight once an introgression event have been detected. 

## <ins>**Installing Dependencies**</ins> <a name="dependencies"></a>

We use conda to create a virtual environment for installing all needed dependencies. We provide a yml file, which users can use to setup an identical environment. To do so, run:
```bash
conda env create -f ghost_int_env.yml
```


## <ins>**Running Ghost Buster**</ins> <a name="running"></a>

## <ins>**Arguments**</ins> <a name="arguments"></a>
These arguments can also be referenced in your command line using the -h flag.
```bash
python Ghost_Buster.py -h
```

UNDER CONSTRUCTION
| Short flag | Long flag         | Description | Required? | Default value |
|------------|-------------------|-------------|-----------|---------------|
| -i         | --input           | Path to fasta file(s) to clean. You can provide the path to a single fasta file or a directory containing multiple fasta files. Squeakuences will not search subdirectories. This can be the full path or relative to the squeakuences.py file location. | Yes | NA |
| -o         | --output          | Path to output folder where files generated by Squeakuences will be written. This can be the full path or relative to the squeakuences.py file location. If this directory path does not exist at runtime, Squeakuences will create it for you. | Yes | NA |
| -l         | --log          | When activated, Squeakuences will generate a log file with processing info from each fasta file cleaned. | No | NA |
| -f         | --addFileName          | When activated, Squeakuences will add the file name to the beginning of all sequences cleaned. | No | NA |

## <ins>**Examples**</ins> <a name="examples"></a>


[Back to Top](#top)





