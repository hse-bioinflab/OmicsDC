# Introduction

This repository serves as a comprehensive framework for conducting deep learning (DL) studies at the HSE bioinformatics
lab. Our research often involves correlating biological experiment results with other data sources, such as DNA
sequences or transcription factor binding sites, to gain insights into biological processes. To achieve this, we use
overlapping annotations or train DL models to identify learned dependencies. This repository primarily focuses on the
latter approach.

Specifically, our framework provides the following functionalities:

* Downloading genome assemblies (in fasta format) and omics data (in BED format)
* Training common types of DL models based on the downloaded data to predict target experimental data 
* Making genome-wide predictions using the trained model
* Interpreting trained models (not applicable to all)

With this framework, we aim to simplify the process of DL-based research and provide a standardized approach for
conducting experiments in our lab.

# Setup

1. Establish an SSH connection with port forwarding using the following command:

```bash
# Forwarding local port 8888 to remote port 8888
ssh -L localhost:8888:localhost:8888 -p <server-port> <username>@<server-address>
```

Replace <server-address> with the address of the server you want to connect to, and <server-port> with the port number.

3. Clone the repository to the server:

```bash
git clone git@github.com:alnfedorov/AI-center.git
```

3. Build the Docker container with all the necessary dependencies:

```bash
cd AI-center/docker
docker build -t ai-center:1.0 .
cd ..
```

4. Launch the JupyterLab frontend:

```bash
docker run -it --rm -p 8888:8888 -v $(pwd):/workspace/ ai-center:1.0
```

Check your console for a URL starting with http://127.0.0.1:8888/lab?token=.
Copy the URL and open it in your browser.

Remember to keep the SSH session running while you work. You may find tmux useful.

# Workflow

The typical workflow consists of three main steps:
* Download target genomic data using data/loader.ipynb.
* Train a DL model by choosing a Jupyter notebook from the `models` zoo.
* Perform post-training analysis, including interpretation and whole-genome prediction. See details in each notebook.
