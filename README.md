# Omics Data Loader

The Omics Data Loader is a Python library for processing omics data. It provides a convenient interface for loading and manipulating data from various omics experiments.

## Features

- Processing omics data
- File handling and management
- Filtering and matching experiments
- Compressing data into archives

## Installation

You can clone the Omics Data Loader repository using Git:

```bash
git clone https://github.com/hse-bioinflab/data-loader
```

## Usage

Import the `omics` function from the `loader` module to start processing omics data:

```python
from loader import omics

result = omics(expid=None, assembly=['hg38'], assembly_threshold='05', antigen_class=None,
               antigen=None, cell_type=None, cell=None, output_path=Path("./storage/"))
```	 
			   
## Documentation

TBA

## Contributing

Contributions to the Omics Data Loader are welcome! If you find any issues or have suggestions for improvements, please submit a pull request or open an issue in the GitHub repository.


