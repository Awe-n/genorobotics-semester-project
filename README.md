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
    - [How to Install and Run the Pipeline](#how-to-install-and-run-the-pipeline)
  - [Notes](#notes)

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