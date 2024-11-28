# Cross-Beta predictor

Welcome to the Cross-Beta RF predictor, an Extratrees-based predictor trained to identify naturally occurring cross-beta-forming amyloids. While this algorithm has been trained against no-forming disordered protein regions, it also has shown very good results in making the difference between regular proteins and amyloids.

This is a version usable with command lines and uses multiprocessing. Below is a description of all the arguments you can use to run it as well as the different conditions to make it work.

Most used packages are already given when installing Python, however, some of them might need to be updated or installed.

An online version is also provided at https://bioinfo.crbm.cnrs.fr/index.php?route=tools&tool=35.


## Table of content

1. [Libraries and dependancies](#Libraries-and-dependencies)
2. [Files](#Files)
    1. [`CB_RF_pred.py`, the main file](#CB_RF_pred.py)
    2. [`utils` folder](#utils)
3. [Arguments and how to use the Cross-Beta RF pred with command lines](#Arguments)
    1. [Options](#options)
    2. [Usage](#usage)
4. [Example of utilization](#Example-of-utilisation)
    1. [Check the installation of Cross-Beta RF predictor](#check-install)
    2. [Run the program with a single protein sequence](#prot-sequence)
    3. [Run the program with a CSV file](#csv)
    4. [Run the program with a fasta file (work also for multifasta)](#fasta)
5. [Contact or cite us](#Contact-us)

  
<a  name="Libraries-and-dependencies"></a>
## 1. Libraries and dependancies

Cross-Beta predictor has been tested for `python 3.8` and `python 3.12`

*  `pickle (v 0.0.12)`  *(installed by default since python 3.7)*
*  `numpy (v 1.24.2)`
*  `pandas (v 2.0.0)`
*  `matplotlib (v 3.7.1)`
*  `scipy (1.8.0)`
*  `scikit-learn (1.3.1)` *Using a different version may invalidate your results, as the model used was built in version 1.3.1.*

This work also uses a protein disorder predictor [IUPred3](https://iupred3.elte.hu/) which you can download by requesting their [website](https://iupred3.elte.hu/download_new). 

Note that IUPred3 is under academic license.

Once you have downloaded IUPRed3, you can put it in the `utils` folder as:

```
- Cross-Beta_predictor
|
| - utils
|  |
|  | - iupred3 <------
|  |  |
|  |  | ...
|  |
|  | check_install.py
|  | fold_pred.py
|  | ...
|
| ...
```

  
<a  name="Files"></a>
## 2. Files

<a  name="CB_RF_pred.py"></a>
### i. `CB_RF_pred.py`, the main file

`CB_RF_pred.py`: The main Python script. Will parse the given arguments by the command line, verify them and use them to run the main function: `Cross_Beta_RF_pred()`. In this function, you will find the extraction of the given sequences and label them to run the prediction function: `run_prediction_fasta_longseq()`. As some used algorithms only work with the 20 essential amino acids code, a regex filter `r"^[ARNDCQEGHILKMFPSTWYV]*$"`is applied on all the given sequences. **All the sequences that don't match this expression won't be predicted**


<a  name="utils"></a>
### ii. `utils` folder

1.  `utils.py`: the file where you will find the sequence extraction function as well as the function for drawing graphs and saving data in .csv format. You will also find a function to detect if the script must run on a folder or on a file and a function to check if a given argument has a None value.

2.  `fold_pred.py`: the file handling all the disorder predictions by using IUpred3.

3.  `get_features.py`: the file where all function generating and extracting features from amino acid sequences are.

4.  `progress_bar.py`: the script where functions are about generating a progress bar and computing the estimated remaining time.

<a  name="Arguments"></a>
## 3. Arguments and how to use the Cross-Beta RF pred with command lines

<a  name="options"></a>
### i. Options

`-h`, `--help`: show this help message and exit

`-ci`, `--check_install`: Check if all folders and files are correctly named and in the right place. Check if all dependencies are importable

`-i INPUT`, `--input INPUT`: Input file path, can be file or folder or an amino acid sequence

`-it INPUT_TYPE`, `--input_type INPUT_TYPE`: Define the input type (sequence, CSV or fasta) (default: 'fasta').

`-nc NAME_COL`, `--name_col NAME_COL`: The name of the column containing the sequence ids. **Only needed if the input type is CSV.**

`-sc SEQ_COL`, `--seq_col SEQ_COL`: The name of the column containing the sequences. **Only needed if the input type is CSV.**

`-t THRESHOLD`, `--threshold THRESHOLD`: Classification threshold. Must be between 0 and 1 (Default: 0.54)

`-g DRAW_GRAPH`, `--draw_graph DRAW_GRAPH`: Active or not the creation of graph for each result in the given output path

`-o OUTPUT`, `--output OUTPUT`: Output file name, the output results will always be saved in the 'results/' folder (default: prediction_result.csv)

 
<a  name="usage"></a>
### ii. Usage

3 types of inputs are handled: path to a folder, path to a file, or an amino acid sequence. To make the program work for all those cases, you must specify the type of the given input: `sequence`, `CSV`, or `fasta`. From this, if the input is a folder, the program will run on all the `CSV` or `fasta` present in it.

If the type is `CSV`, you must enter the 2 additional arguments `-nc` (the name of the column where to find a sequence ID or equivalent) and `-sc` (the name of the column where to find the sequences to run the prediction on).

Then you can choose a threshold value going from 0 to 1 for the prediction. All results strictly superior to the threshold will be considered amyloidogenic. The default and recommended threshold is set to 0.54.

Two outputs are generated, the CSV result file name is determined by the `-o` argument (default set to prediction_result.csv). It contains the ID of the sequence, the length, the global prediction score, the amyloid regions (of minimum 15 amino acids with a score > threshold) and the detailed score for all the amino acids of the sequence. The other output optionally generated is the graphic view of the amyloidogenicity score for all the amino acids of the sequence. The name of the graphic is based on the id of the corresponding sequence. This option can be activated or not by entering `True` or `False` to the `-g` argument (Default set to `False`).

 
<a  name="Example-of-utilisation"></a>
## 4. Example of utilisation
<a  name="check-install"></a>
### i. Check the installation of Cross-Beta RF predictor

```bash
python3  CB_RF_pred.py  \
-ci
```

<a  name="prot-sequence"></a>
### ii. Run the program with a single protein sequence

```bash
python3  CB_RF_pred.py  \
-i  MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYEPRGRGGGRGGRGGMGGSDRGGFNKFGGPRDQGSRHDSEQDNSDNNTIFVQGLGENVTIESVADYFKQIGIIKTNKKTGQPMINLYTDRETGKLKGEATVSFDDPPSAKAAIDWFDGKEFSGNPIKVSFATRRADFNRGGGNGRGGRGRGGPMGRGGYGGGGSGGGGRGGFPSGGGGGGGQQRAGDWKCPNPTCENMNFSWRNECNQCKAPKPDGPGGGPGGSHMGGNYGDDRRGGRGGYDRGGYRGRGGDRGGFRGGRGGGDRGGFGPGKMDSRGEHRQDRRERPY  \
-it  sequence  \
-g
```

The given input is an amino acid sequence and a result graphic will be drawn. All results will be saved in the `/results/` directory in the `Cross-Beta_RF_pred` folder.

  
<a  name="CSV"></a>
### iii. Run the program with a CSV file

```bash
python3  CB_RF_pred.py  \
-i  example/Test_csv.csv  \
-it  'csv'  \
-nc  'Label'  \
-sc  'Sequence'
```

Run the prediction on the CSV available in `/example/Test_csv.csv` without drawing the resulting graph. For a CSV, label and sequence column names must be specified using `-nc` and `-sc`.

  
<a  name="fasta"></a>
### iv. Run the program with a fasta file (work also for multifasta)

```bash
python3  CB_RF_pred.py  \
-i  example/Test_fasta.fasta  \
-it  'fasta'  \
-g  \
-o  fasta_result.csv
```
Run the prediction on the fasta file available in `/example/Test_fasta.fasta`, create graphic results and store the CSV result file as 'fasta_result.csv'

<a  name="Contact-us"></a>
## 5. Contact or cite us
For any issue with the program, you can address a mail with `[Cross-Beta RF]` in subject to `valentin.gonay@crbm.cnrs.fr`

**Citing Cross-Beta predictor:**

Valentin Gonay, Michael P. Dunne, Javier Caceres-Delpiano, & Andrey V. Kajava. (2024). Developing machine-learning-based amyloid predictors with Cross-Beta DB. bioRxiv, 2024.02.12.579644. https://doi.org/10.1101/2024.02.12.579644

**Citing IUPred3:**

Gábor Erdős, Mátyás Pajkos, Zsuzsanna Dosztányi IUPred3 - improved prediction of protein disorder with a focus on specific user applications Nucleic Acids Research 2021, Submitted Bálint Mészáros, Gábor Erdős, Zsuzsanna Dosztányi (2018) IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding. Nucleic Acids Research;46(W1):W329-W337.

### Author
Valentin Gonay
[GitHub](https://github.com/Valentin-Gonay)

