# GenoRobotics Project - Full Pipeline

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
    - [How to Install BLASTn on Windows](#how-to-install-blastn-on-windows)
    - [How to Install and Run the Pipeline](#how-to-install-and-run-the-pipeline)
      - [Setting Up the Environment](#setting-up-the-environment)
      - [Running the Pipeline Scripts Using Jupyter Notebooks](#running-the-pipeline-scripts-using-jupyter-notebooks)
      - [Running the Pipeline Scripts Using the Terminal](#running-the-pipeline-scripts-using-the-terminal)
  - [Output Files Locations](#output-files-locations)
    - [Log Files Locations](#log-files-locations)
    - [Consensus, Identification, and Streaming Pipeline Iteration Results](#consensus-identification-and-streaming-pipeline-iteration-results)
  - [Extra Notes](#extra-notes)
    - [Test Dataset](#test-dataset)
    - [Testing the Pipeline](#testing-the-pipeline)
  - [Authors](#authors)

## Description

This repository hosts the code for the complete pipeline of the GenoRobotics project, offering a suite of bioinformatics pipelines for DNA sequence analysis. The pipelines are designed to cater to both single file analysis and large-scale expeditions, with options for standard and streaming processing.

Our paper describing the pipelines and their performance is available [here](https://www.overleaf.com/read/mbqjxqjxqjxq).

### Available Scripts

- **Single File Standard Pipeline:**
  - `standard-pipeline.py`: For processing a single input file using the standard bioinformatics pipeline.
  - `standard-detailed-pipeline.ipynb`: An interactive and detailed version for single file processing.
- **Expedition Standard Pipeline:**
  - `expedition-pipeline.py`: For processing multiple files from an expedition using the standard pipeline.
  - `expedition-detailed-pipeline.ipynb`: An interactive and detailed version for expedition processing.
- **Single File Streaming Pipeline:**
  - `streaming-pipeline.py`: For processing a single input file using the streaming bioinformatics pipeline (optimized for efficiency).
  - `streaming-detailed-pipeline.ipynb`: An interactive and detailed version for single file streaming.

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

To install BLASTn on Windows, follow these steps carefully:

1. **Download BLAST+ Executables:**
   - Visit the [NCBI BLAST download page](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
   - Look for the Windows version of the BLAST+ executables (e.g., `ncbi-blast-2.12.0+-win64.exe`). Click to download the installer.

2. **Install BLAST:**
   - After downloading, run the `.exe` installer. 
   - Follow the installation wizard, which will guide you through the setup process, including the installation location for BLAST.

3. **Unzip the Downloaded File:**
   - If the downloaded BLAST+ package is a zip file, unzip it to your desired location.

4. **Set Up Environment Variables:**
   - Add the path to the `bin` folder of the BLAST installation to the PATH environment variable.
   - Create a new environment variable named `BLASTDB` and set its value to the path of the `db` folder in your BLAST installation.
   - To set these environment variables:
     - Right-click on 'This PC' or 'Computer' on the desktop or in File Explorer, then select 'Properties'.
     - Click 'Advanced system settings' and then the 'Environment Variables' button.
     - Under 'System variables', find and select the PATH variable, then click 'Edit' to add the BLAST `bin` directory.
     - Click 'New' to create the `BLASTDB` variable and set its value to the BLAST `db` directory.
     - Click 'OK' to save changes and close all dialogs.

5. **Test the Installation:**
   - Open Command Prompt and type `blastn -version`.
   - This should display the installed version of BLASTn, confirming the successful installation.

Note: The installation process may vary slightly based on the version of the BLAST+ executables. Always follow the instructions provided with the downloaded package. If you encounter any issues, ensure that the paths in your environment variables are correct.

### How to Install and Run the Pipeline

The pipeline is designed to run within an Anaconda environment. The required environment configuration is specified in the file `genorobotics_pipeline.yml`.

#### Setting Up the Environment

1. **Create the Anaconda Environment:**
   - Open your terminal and navigate to the project's folder.
   - Run the command `conda env create -f genorobotics_pipeline.yml`. This will create an environment named `genorobotics_pipeline`. (Ensure Anaconda is installed beforehand.)

2. **Activate the Environment:**
   - In your terminal, type `conda activate genorobotics_pipeline` to activate the newly created environment.

#### Running the Pipeline Scripts Using Jupyter Notebooks

Each pipeline script is paired with a corresponding Jupyter Notebook for interactive use, located in the `src` folder. These notebooks provide an engaging way to execute and visualize each step of the pipeline.

Available Notebooks:
- `standard-detailed-pipeline.ipynb`
- `expedition-detailed-pipeline.ipynb`
- `streaming-detailed-pipeline.ipynb`

To use these notebooks:
- Open them in Jupyter Lab or Jupyter Notebook.
- Execute the cells in order to run the pipeline interactively.
- Make sure to use the correct Anaconda environment when running the notebooks.

#### Running the Pipeline Scripts Using the Terminal

For those preferring a command-line approach, the pipeline scripts can also be executed directly from the terminal.

- **For a Single Input File:**
  1. Navigate to the `src` directory: `cd src`.
  2. Execute the pipeline with the command:
     ```bash
     python3 pipeline.py <name_of_the_input_fastq_file> <windows (True/False, Optional)> <db (Optional)>
     ```
     - The input FASTQ file should be located in `src/assets/input`.
     - `windows`: Optional. Set to `True` for Windows with BLASTn installed (default: `False`).
     - `db`: Optional. Specify the database for BLASTn search. Defaults to using all four databases (ITS, matK, rbcL, trnL) if unspecified.

- **For an Entire Expedition:**
  1. Ensure input files are in `src/assets/input/<name_of_the_expedition>`, following the ONT sequencer's output structure.
  2. Run the expedition pipeline:
     ```bash
     python3 pipeline_expedition.py <name_of_the_expedition> <windows> <name_of_the_consensus_method>
     ```
     - `windows`: Optional (default: `False`). Set to `True` if BLASTn is installed on Windows.
     - `name_of_the_consensus_method`: Optional (default: `80_20`). Choices are `80_20` or `streaming`.

- **For the Streaming Pipeline:**
  1. Place the input FASTQ file in `src/assets/input`.
  2. Run the streaming pipeline:
     ```bash
     python3 streaming_pipeline.py <input_fastq_filename> <windows (True/False)> <db> <streaming_method> <consensus_method> <identification_method>
     ```
     - `<input_fastq_filename>`: Name of the input FASTQ file.
     - `<windows>`: Required. Set to `True` for Windows with BLASTn.
     - `<db>`: Required. Specify the BLASTn search database.
     - `<streaming_method>`: Required. Choose the streaming method.
     - `<consensus_method>`: Required. Choose the consensus method.
     - `<identification_method>`: Required. Choose the identification method.
     - Example: `python3 streaming_pipeline.py rbcL_Qiagen_tomato_5000.fastq False matK basic_streaming 80_20_best_sequence blastn`.

To exit the Anaconda environment after running the pipeline, enter `conda deactivate` in the terminal.

## Output Files Locations

The pipeline is designed to log all operations into specific files, ensuring a clean interface in the terminal and notebooks.

### Log Files Locations

1. Log Files for `pipeline.py`
- **Consensus Logs:**
  - Location: `src/assets/output/post/<name_of_the_input_fastq_file>/consensus`.
  - Filename: `<name_of_the_input_fastq_file>_consensus.log`.
- **Identification Logs:**
  - Location: `src/assets/output/post/<name_of_the_input_fastq_file>/identification`.
  - Filename: `<name_of_the_input_fastq_file>_identification.log`.

2. Log Files for `pipeline_expedition.py`
- Individual log files for each input file.
- Location: `src/assets/log/name_of_the_expedition/`.
- Filename: `<name_of_the_expedition>.log`.

3. Log Files for `streaming_pipeline.py`
  - **Consensus Logs:**
    - Location: `src/assets/output/streaming/<name_of_the_input_fastq_file>`.
    - Filename: `<name_of_the_input_fastq_file>_consensus.log`.

### Consensus, Identification, and Streaming Pipeline Iteration Results

1. **Consensus Sequence:**
  - Location: `src/assets/output/post/<name_of_the_input_fastq_file>/consensus`.
  - Filename: `<name_of_the_input_fastq_file>_final_consensus.fasta`.
2. **BLASTn Results:**
  - Location: `src/assets/output/post/<name_of_the_input_fastq_file>/identification`.
  - Files:
    - Raw BLASTn results: `<name_of_the_db>.xml`.
    - Filtered BLASTn results: `<name_of_the_input_fastq_file>_identification_results.csv`.
3. **Streaming Pipeline Iteration Results:**
  - Location for each block: `src/assets/output/streaming/<name_of_the_input_fastq_file>/block_X`.
  - Contents:
    - BLASTn results in `block_X/identification`.
    - Consensus sequences in `block_X/consensus`.

## Extra Notes

### Test Dataset

Currently, a small dataset called `rbcL_Qiagen_tomato_5000.fastq` (containing 5000 sequences, ~5 MB) can be found in the `src/assets/input` folder. This is a sampled version from the entire dataset `rbcL_Qiagen_tomato.fastq` (containing ~200k sequences, ~204 MB). 

### Testing the Pipeline
- To verify the installation and functionality of the pipeline, use the `rbcL_Qiagen_tomato_5000.fastq` dataset.
- Run the following command in the terminal for a single input file test:
  ```bash
  python3 pipeline.py rbcL_Qiagen_tomato_5000.fastq
  ```
  For Windows with BLASTn installed, use:
  ```bash
  python3 pipeline.py rbcL_Qiagen_tomato_5000.fastq True
  ```

## Authors

- **Awen Kidel Peña--Albert**
- **Emilien Ordonneau**