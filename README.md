EMILIEN README




# GenoRobotics Full Pipeline

## Table of Contents

- [Description](#description)


## Description

This repository contains the code for the full pipeline of the GenoRobotics project. 
The pipeline is composed of 3 main steps:
1. Sequence Preprocessing
2. Consensus Sequence Generation
3. Sequence Identification

### Sequence Preprocessing

The sequence preprocessing step is composed of 3 sub-steps:
 - quality filtering
 - primer trimming
 - length filtering

### Consensus Sequence Generation

The consensus sequence generation step highly depends on the consensus method selected.
Currently, the following consensus methods are available:
 - 80/20 consensus, which generates a consensus by first aligning the 20% longest sequences, creating a first consensus sequence, then aligning the remaining 80% sequences to the first consensus sequence, creating a second, usable consensus sequence.
 - streaming consensus, which generates a consensus by first aligning the sequences as they come in, creating a first consensus sequence, then aligning the remaining sequences to the first consensus sequence, creating a second consensus sequence.

### Sequence Identification

The sequence identification step is composed of 2 sub-steps:
 - BLASTn search : the consensus sequence is searched against the NCBI nucleotide database using BLASTn
 - BLASTn filtering : the BLASTn results are filtered to keep only the best hit for each species

## Installation and Usage

### Requirements

The following requirements are needed to run the pipeline:
 - Conda 4.8.3
 - BLAST 2.X.0
 - WSL 2 (Windows only)


### Notes for Windows Users

It is highly recommanded to install BLASTn on WSL and not on Windows, and to run the pipeline on WSL.
The pipeline was developed and tested on a Windows 11 machine using WSL 2 (Windows Subsystem for Linux).
You can find more information about WSL 2 here: https://docs.microsoft.com/en-us/windows/wsl/about
If you insist on installing BLASTn on Windows, you can then run the pipeline on Windows Powershell, and by changing the parameter 'windows' to True when calling the .py scripts, or in the .ipynb notebooks. This will still launch every command in a WSL terminal, so you will still need to have WSL installed, but the BLASTN commands will be launched in a Windows terminal.

### How to install BLASTn on Unix

1. **Download BLAST:**
   - Visit the NCBI BLAST FTP site: [ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
   - Choose the appropriate BLAST+ package for your system (look for a file ending in `-x64-linux.tar.gz` for Linux or `-x64-macosx.tar.gz` for macOS).

2. **Download via Terminal:**
   - Open your terminal.
   - Use `wget` or `curl` to download the BLAST+ package. Replace the URL with the one you found for the correct version.
     ```bash
     wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz
     ```
     or
     ```bash
     curl -O ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz
     ```

3. **Extract the Package:**
   - Once downloaded, extract the tarball using:
     ```bash
     tar -zxvf ncbi-blast-2.12.0+-x64-linux.tar.gz
     ```
   - This will create a directory with the BLAST+ executables.

4. **Add BLAST to Your PATH:**
   - To run BLASTn from any location, you need to add it to your PATH. Replace `/path/to/ncbi-blast-2.12.0+/bin` with the actual path to the `bin` directory inside the extracted folder.
     ```bash
     export PATH=$PATH:/path/to/ncbi-blast-2.12.0+/bin
     ```
   - You also need to add the db folder to your BLASTDB environment variable. Replace `/path/to/ncbi-blast-2.12.0+/db` with the actual path to the `db` directory inside the extracted folder.
     ```bash
     export BLASTDB=$BLASTDB:/path/to/ncbi-blast-2.12.0+/db
     ```
   - Consider adding these lines to your `.bashrc` or `.bash_profile` (or the equivalent for your shell) to make this change permanent.

5. **Test the Installation:**
   - To ensure BLASTn is installed correctly, run:
     ```bash
     blastn -version
     ```
   - This should return the installed version of BLASTn.

6. **Update your System (Optional):**
   - Ensure your system's package list is updated and all dependencies are satisfied:
     ```bash
     sudo apt-get update      # For Debian/Ubuntu
     sudo yum update          # For CentOS/RedHat
     ```

### How to install BLASTn on Windows

Here are the following steps to follow to install BLASTn on Windows:
 - download the BLAST+ executables from the NCBI website: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
 - unzip the downloaded file
 - add the path to the bin folder of the unzipped folder to the PATH environment variable
 - add the path to the db folder of the unzipped folder to a new BLASTDB environment variable
 - open a new terminal and type `blastn` to check if the installation was successful

#### Installation Steps:

1. **Download and install BLAST:**
   - Visit the NCBI BLAST download page: [BLAST+ executables](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
   - Look for the Windows version of the BLAST+ executables (e.g., `ncbi-blast-2.12.0+-win64.exe`). Click the link to download the installer.
   - Once downloaded, run the installer. This will typically be a `.exe` file.
   - Follow the installation wizard. It will guide you through the setup process, including where to install BLAST.

2. **Set Up Environment Variables (if necessary):**
   - The installer may give you an option to add BLAST to your system's PATH. If it does, select this option.
   - If not, you can add it manually:
     - Right-click on 'This PC' or 'Computer' on the desktop or in File Explorer.
     - Click 'Properties'.
     - Click 'Advanced system settings'.
     - In the System Properties window, click the 'Environment Variables' button.
     - In the Environment Variables window, under 'System variables', find the PATH variable and select it. Click 'Edit'.
     - In the Edit Environment Variable window, click 'New' and add the path to the BLAST executable (the `bin` directory in your BLAST installation).
     - In the Environment Variables window, click 'New' to add a new variable. Name it `BLASTDB` and set the value to the path to the BLAST database (the `db` directory in your BLAST installation).
     - Click 'OK' to close the dialogs.

3. **Test the Installation:**
   - Open Command Prompt and type:
     ```cmd
     blastn -version
     ```
   - This should return the installed version of BLASTn, indicating that the installation was successful.

### How to Install and Run the Pipeline

The Anaconda environment used to run the pipeline can be found in the file genorobotics_pipeline.yml.

Here are the following steps to follow to run it:
 - open the current folder in your terminal
 - type `conda env create -f genorobotics_pipeline.yml` in the terminal to create the environment genorobotics_pipeline (Anaconda must be installed)
 - type `conda activate genorobotics_pipeline` in the terminal to activate the environment
 - type `cd src` in the terminal
 - type `python3 pipeline.py <name_of_the_input_fastq_file> <windows (True/False, Optional)> <db (Optional)>` in the terminal to run the pipeline for a single input file. The input file must be located in the `src/assets/input` folder. The windows parameter is optional and is set to False by default. The db parameter is optional, if not specified, the pipeline will use all four databases (ITS, matK, rbcL, trnL) for the BLASTn search. If specified, the pipeline will only use the specified database for the BLASTn search. The db parameter must be one of the following: ITS, matK, rbcL, trnL.
 - type `python3 pipeline_expedition.py <name_of_the_expedition> <windows> <name_of_the_consensus_method>` in the terminal to run the pipeline for an entire expedition. The input files must be located in the `src/assets/input/name_of_the_expedition` folder, as outputted by the ONT sequencer (the pipeline will look for the `fastq_pass` subfolder in this folder). The windows parameter is optional and is set to False by default. It is to be set to True if BLASTn is installed on Windows, as mentionned in the Notes for Windows Users section. The name_of_the_consensus_method parameter is optional and is set to 80_20 by default. The name_of_the_consensus_method parameter must be one of the following: 80_20, streaming.
 - to exit the environment, type `conda deactivate` in the terminal

Alternatively, you can run the pipeline using the Jupyter Notebooks located in the `src` folder. The notebooks are named `detailed-pipeline.ipynb` and `detailed-pipeline_expedition.ipynb` and contain the same code as the pipeline.py and pipeline_expedition.py scripts, respectively. The notebooks can be run using Jupyter Lab or Jupyter Notebook.

## Output

To not overload the terminal/notebooks, everything is logged in specific logfiles. 

For the pipeline.py script, two logfiles are generated, one for the consensus and one for the identification. The logfiles are located in the `src/assets/output/consensus/<name_of_the_input_fastq_file>` and `src/assets/output/identification/<name_of_the_input_fastq_file>` folders and are named `<name_of_the_input_fastq_file>_consensus.log` and `<name_of_the_input_fastq_file>_identification.log`, respectively.

For the pipeline_expedition.py script, a single logfile is generated for each input file. The logfile is located in the `src/assets/log/name_of_the_expedition/` folder and is named `<name_of_the_expedition>.log`.

In `src/assets/output/consensus/<name_of_the_input_fastq_file>`, you can find the consensus sequence in the `<name_of_the_input_fastq_file>_final_consensus.fasta` file.

In `src/assets/output/identification/<name_of_the_input_fastq_file>`, you can find the BLASTn results in the `<name_of_the_db>.xml` file and the filtered BLASTn results in the `<name_of_the_input_fastq_file>_identification_results.csv` file. 

## Notes

Currently, a small dataset called `rbcL_Qiagen_tomato_5000.fastq` (containing 5000 sequences, ~5 MB) can be found in the `src/assets/input` folder. This is a sampled version from the entire dataset `rbcL_Qiagen_tomato.fastq` (containing ~200k sequences, ~204 MB). 

This dataset is used in the test loop of the pipeline.py script in the run example and can be used to check if the environment is correctly installed and if the pipeline is working as intended. 

You can test the single input file pipeline by running the following command in the terminal: `python3 pipeline.py rbcL_Qiagen_tomato_5000.fastq` (or `python3 pipeline.py rbcL_Qiagen_tomato_5000.fastq True` if BLASTn is installed on Windows).

### Additional Notes

- For MacOS, in your .bash_profile file, make sure to put the new lines about BLAST at the beginning of the file
- the run_command and run_bash_command functions are supposed to run well but they are subject to potential bug, if you encounter a bug when running minimap2, racon or BLASTn command, try to debug those first









AWEN README










# GenoRobotics Full Pipeline

## Table of Contents

- [GenoRobotics Full Pipeline](#genorobotics-full-pipeline)
  - [Table of Contents](#table-of-contents)
  - [Description](#description)
    - [Sequence Preprocessing](#sequence-preprocessing)
    - [Consensus Sequence Generation](#consensus-sequence-generation)
    - [Sequence Identification](#sequence-identification)
  - [Installation and Usage](#installation-and-usage)
    - [Requirements](#requirements)
    - [Notes for Windows Users](#notes-for-windows-users)
    - [How to install BLASTn on Unix](#how-to-install-blastn-on-unix)
    - [How to install BLASTn on Windows](#how-to-install-blastn-on-windows)
      - [Installation Steps:](#installation-steps)
    - [How to Install and Run the Pipeline](#how-to-install-and-run-the-pipeline)
  - [Output](#output)
  - [Notes](#notes)
    - [Additional Notes](#additional-notes)
- [GenoRobotics Full Pipeline](#genorobotics-full-pipeline-1)
  - [Table of Contents](#table-of-contents-1)
  - [Description](#description-1)
    - [Sequence Preprocessing](#sequence-preprocessing-1)
    - [Consensus Sequence Generation](#consensus-sequence-generation-1)
    - [Sequence Identification](#sequence-identification-1)
  - [Installation and Usage](#installation-and-usage-1)
    - [Requirements](#requirements-1)
    - [How to Install and Run the Pipeline](#how-to-install-and-run-the-pipeline-1)
  - [Notes](#notes-1)

## Description

This repository contains the code for the full pipeline of the GenoRobotics project. 
The pipeline is composed of 3 main steps:
1. Sequence Preprocessing
2. Consensus Sequence Generation
3. Sequence Identification

### Sequence Preprocessing

The sequence preprocessing step is composed of 3 sub-steps:
 - quality filtering
 - primer trimming
 - length filtering

### Consensus Sequence Generation

The consensus sequence generation step highly depends on the consensus method selected.

### Sequence Identification

The sequence identification step is composed of 2 sub-steps:
 - BLASTn search

## Installation and Usage

### Requirements

The following requirements are needed to run the pipeline:
 - Conda 4.8.3

### How to Install and Run the Pipeline

The Anaconda environment used to run the pipeline can be found in the file genorobotics_pipeline.yml.

Here are the following steps to follow to run it:
 - open the current folder in your terminal
 - type `conda env create -f genorobotics_pipeline.yml` in the terminal to create the environment genorobotics_pipeline (Anaconda must be installed)
 - type `conda activate genorobotics_pipeline` in the terminal to activate the environment
 - type `cd src` in the terminal
 - type `python pipeline.py` in the terminal to run full pipeline with a test input dataset
 - to exit the environment, type `conda deactivate` in the terminal

## Notes

Currently, a small dataset called `rbcL_Qiagen_tomato_5000.fastq` (containing 5000 sequences, ~5 MB) can be found in the `src/assets/input` folder. This is a sampled version from the entire dataset `rbcL_Qiagen_tomato.fastq` (containing ~200k sequences, ~204 MB). 

This dataset is used in the test loop of the pipeline.py script in the run example and can be used to check if the environment is correctly installed.