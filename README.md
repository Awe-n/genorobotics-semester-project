EMILIEN README



# GenoRobotics Project - Full Pipeline

## Table of Contents

[[toc]]

## Description

This repository hosts the code for the complete pipeline of the GenoRobotics project, offering a suite of bioinformatics pipelines for DNA sequence analysis. The pipelines are designed to cater to both single file analysis and large-scale expeditions, with options for standard and streaming processing.

### Available Scripts

- **Single File Standard Pipeline:**
  - `standard-pipeline.py`: For processing a single input file using the standard bioinformatics pipeline.
  - `standard-detailed-pipeline.ipynb`: An interactive and detailed version for single file processing.
- **Expedition Standard Pipeline:**
  - `expedition-pipeline.py`: For processing multiple files from an expedition using the standard pipeline.
  - `expedition-detailed-pipeline.ipynb`: An interactive and detailed version for expedition processing.
- **Single File Streaming Pipeline:**
  - `streaming-pipeline.py`: For processing a single input file using the streaming bioinformatics pipeline (optimized for efficiency).

To run the full pipeline, navigate to the `src` folder and choose the appropriate script based on your analysis needs.

### Pipeline Overview

#### Standard Bioinformatics Pipeline

The standard pipeline comprises three primary steps:

1. **Sequence Preprocessing:** *(Currently not utilized, but code is maintained for future implementation.)*
2. **Consensus Sequence Generation:** Relies on the selected consensus method, with options like straightforward consensus, 80/20 best sequence consensus, and 80/20 longest sequence consensus.
3. **Sequence Identification:** Utilizes BLASTn for identification, which involves a search against the NCBI nucleotide database and subsequent filtering to retain the best hit per species.

#### Streaming Bioinformatics Pipeline

The streaming pipeline parallels the standard pipeline but operates on a stream of sequences, allowing for processing during the sequencing process. It employs random sampling and block-wise iteration for enhanced speed and reduced memory usage. While not yet implementable during live sequencing, it offers significant performance improvements.

### Detailed Pipeline Description

For an in-depth understanding of each step in the pipelines, please refer to our accompanying paper.

## Installation and Usage

### Requirements

1. **Conda 4.8.3**
2. **BLAST 2.X.0**
3. **WSL 2 (Windows only)**

### Notes for Windows Users

It is highly recommended to install BLASTn on WSL rather than directly on Windows, and to run the pipeline using WSL. The pipeline was developed and tested on a Windows 11 machine with WSL 2 (Windows Subsystem for Linux). For more information about WSL 2, visit [Microsoft's official documentation](https://docs.microsoft.com/en-us/windows/wsl/about).

If you choose to install BLASTn on Windows, you can still run the pipeline in Windows Powershell. Make sure to set the 'windows' parameter to True when calling the `.py` scripts or in the `.ipynb` notebooks. Note that commands will still launch in a WSL terminal, so WSL installation remains a requirement.

### Notes for macOS Users

- A very useful tutorial to install BLASTn on macOS can be found [here](https://www.youtube.com/watch?v=sMImc_iwX4w).
- In your `.bash_profile`, place the new lines regarding BLAST at the beginning of the file.
- The `run_command` and `run_bash_command` functions should function correctly, but they might encounter issues. If you face any bugs when running minimap2, racon, or BLASTn commands, address those first.

### Installing BLASTn on Unix

1. **Download BLAST:**
   Visit the [NCBI BLAST FTP site](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and choose the appropriate BLAST+ package for your system (look for a file ending in `-x64-linux.tar.gz` for Linux or `-x64-macosx.tar.gz` for macOS).

2. **Download via Terminal:**
   Open your terminal and use `wget` or `curl` to download the BLAST+ package. Replace the URL with the correct version:
   ```bash
   wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.X.0+-x64-linux.tar.gz
   # or
   curl -O ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.X.0+-x64-linux.tar.gz
  ```

3. **Extract the Package:**
   - Use the following command to extract the tarball:
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

### How to Install BLASTn on Windows

To install BLASTn on Windows, follow these steps:

1. **Download BLAST+ Executables:**
   - Visit the [NCBI BLAST download page](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and select the Windows version of the BLAST+ executables (e.g., `ncbi-blast-2.12.0+-win64.exe`). 

2. **Install BLAST:**
   - Run the downloaded `.exe` installer and follow the setup wizard to install BLAST.

3. **Set Up Environment Variables:**
   - If the installer doesn't automatically add BLAST to your system's PATH, do it manually:
     - Right-click 'This PC' or 'Computer' on the desktop or in File Explorer, then click 'Properties'.
     - Select 'Advanced system settings', then click the 'Environment Variables' button.
     - Under 'System variables', find and edit the PATH variable to include the path to the BLAST `bin` directory.
     - Add a new system variable named `BLASTDB` and set its value to the path of the BLAST `db` directory.

4. **Verify Installation:**
   - To confirm the installation, open Command Prompt and type `blastn -version`. This should display the installed version of BLASTn, confirming a successful installation.

Note: If you encounter any issues, ensure the paths to the `bin` and `db` directories in your environment variables are correct. Updating the environment variables is crucial for BLASTn to function properly in Windows.


### How to install BLASTn on Windows

Here are the following steps to follow to install BLASTn on Windows:
1. download the BLAST+ executables from the [NCBI website](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
2. unzip the downloaded file
3. add the path to the bin folder of the unzipped folder to the PATH environment variable
4. add the path to the db folder of the unzipped folder to a new BLASTDB environment variable
5. open a new terminal and type `blastn` to check if the installation was successful

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

To not overload the terminal/notebooks, everything is logged in specific log files. 

For the pipeline.py script, two logfiles are generated, one for the consensus and one for the identification. The logfiles are located in the `src/assets/output/consensus/<name_of_the_input_fastq_file>` and `src/assets/output/identification/<name_of_the_input_fastq_file>` folders and are named `<name_of_the_input_fastq_file>_consensus.log` and `<name_of_the_input_fastq_file>_identification.log`, respectively.

For the pipeline_expedition.py script, a single logfile is generated for each input file. The logfile is located in the `src/assets/log/name_of_the_expedition/` folder and is named `<name_of_the_expedition>.log`.

In `src/assets/output/consensus/<name_of_the_input_fastq_file>`, you can find the consensus sequence in the `<name_of_the_input_fastq_file>_final_consensus.fasta` file.

In `src/assets/output/identification/<name_of_the_input_fastq_file>`, you can find the BLASTn results in the `<name_of_the_db>.xml` file and the filtered BLASTn results in the `<name_of_the_input_fastq_file>_identification_results.csv` file. 

## Notes

Currently, a small dataset called `rbcL_Qiagen_tomato_5000.fastq` (containing 5000 sequences, ~5 MB) can be found in the `src/assets/input` folder. This is a sampled version from the entire dataset `rbcL_Qiagen_tomato.fastq` (containing ~200k sequences, ~204 MB). 

This dataset is used in the test loop of the pipeline.py script in the run example and can be used to check if the environment is correctly installed and if the pipeline is working as intended. 

You can test the single input file pipeline by running the following command in the terminal: `python3 pipeline.py rbcL_Qiagen_tomato_5000.fastq` (or `python3 pipeline.py rbcL_Qiagen_tomato_5000.fastq True` if BLASTn is installed on Windows).











AWEN README










# GenoRobotics Full Pipeline

## Table of Contents

- [GenoRobotics Project - Full Pipeline](#genorobotics-project---full-pipeline)
  - [Table of Contents](#table-of-contents)
  - [Description](#description)
    - [Available Scripts](#available-scripts)
    - [Pipeline Overview](#pipeline-overview)
      - [Standard Bioinformatics Pipeline](#standard-bioinformatics-pipeline)
      - [Streaming Bioinformatics Pipeline](#streaming-bioinformatics-pipeline)
    - [Detailed Pipeline Description](#detailed-pipeline-description)
  - [Installation and Usage](#installation-and-usage)
    - [Requirements](#requirements)
    - [Notes for Windows Users](#notes-for-windows-users)
    - [Notes for macOS Users](#notes-for-macos-users)
    - [Installing BLASTn on Unix](#installing-blastn-on-unix)
- [or](#or)
    - [How to install BLASTn on Windows](#how-to-install-blastn-on-windows)
      - [Installation Steps:](#installation-steps)
    - [How to Install and Run the Pipeline](#how-to-install-and-run-the-pipeline)
  - [Output](#output)
  - [Notes](#notes)
- [GenoRobotics Full Pipeline](#genorobotics-full-pipeline)
  - [Table of Contents](#table-of-contents-1)
  - [Description](#description-1)
    - [Sequence Preprocessing](#sequence-preprocessing)
    - [Consensus Sequence Generation](#consensus-sequence-generation)
    - [Sequence Identification](#sequence-identification)
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