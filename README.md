# AutophiGen
## A bioinformatics tool to automate phylogenetic tree analysis.
#### Video Demo:  <a href = "https://youtu.be/Podl1cscvpE"> Wath full tutorial here</a>
#### Description:
AutophiGen is a general purpose phylogenetic analysis program which can automate and generate phylogenetic trees for referencing purposes. It automates:
* Performing a BLAST search using Biopython module
* doing multiple sequence alignment using MUSCLE
* devise a simple phylogenetic tree using EMBL-EBI APIs
* gives output of phylogenetic tree in GUI using ete3
### DISCLAIMER:
This program is initially made as a final project for Harvard's CS50P coursework. It fulfills basic need but is in no way viable enough to be used for publication and academic writing (as of now!)
## Usage:
AutophiGen is a python program and is open to run using python interpreter. The usage format for AutophiGen is: 
```
    autophigen.py -i [input file] -e [email address]
```
You can also use ```autophigen.py -h``` to look for arguement options.
> [!NOTE]
> Input sequences should be in a fasta file.

> [!NOTE]
> AutophiGen can take multiple sequences as input in a single file as it have parsing functions to handle and process multiple sequences.

> [!IMPORTANT]
> Check requirements.txt to install and update modules and their dependencies to the required version to run the analysis smoothly.

> [!IMPORTANT]
> If facing error [SSL: CETIFICATE_VERIFY_FAILED]. It is error caused by Python 3.x for MacOS as it does not rely on MacOS SSL certificates. Kindly do the following to reolve this issue:
> 1. Navigate to Application/Python 3.x in Finder. You will see a command file shipped with python named "Install certificates.command".
> 2. Run the command file. It will simply install the certificates using pip.
## Known Errors and shortcomings:
The errors are mostly resolvable and are present due to my lack of experience as a beginner. Below mentioned are things which I would like to improve in this program:
1. Running BLAST can be a bit slow due to usage of Biopython's qblast module to access BLAST APIs.
2. Program may give some errors which I have not recieved during my testing runs or need some in-depth technical knowledge to induce and test.
> [!NOTE]
> If facing any errors in the program, kindly inform me of the same by mail on adnanraza3435@gmail.com.
## Thank You for trying and using this program!


    
