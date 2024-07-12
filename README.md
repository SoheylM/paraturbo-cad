# ParaturboCAD Repository

## Overview
This repository contains the ParaturboCAD library, a comprehensive Python-based tool for the automated visualization and design of gas-bearing supported turbocompressors. Developed as part of ongoing research, this tool leverages CadQuery 2, a scripting CAD platform, to streamline the design process traditionally involving multiple software tools and extensive manual intervention.

## Citation
If you use this software in your research, please cite it as follows:

Massoudi, S., Bejjani, J., Horvath, T., Üstün, D., Schiffmann, J. (2024). PARATURBOCAD: AUTOMATED PARAMETRIC GEOMETRY CONSTRUCTION FOR GAS-BEARING SUPPORTED TURBOCOMPRESSOR DESIGN WITH PYTHON. Proceedings of the ASME 2024 International Design Engineering Technical Conferences and Computers and Information in Engineering Conference, DETC2024-143696, August 25-28, Washington, District of Columbia. [Software]. Available at https://github.com/SoheylM/ParaturboCAD.

## Installation

### Prerequisites
- Python 3.9
- Anaconda (Recommended for managing Python packages and environments)

### Steps
1. **Clone the Repository**
```bash
git clone https://github.com/YourGitHub/ParaturboCAD.git
cd paraturbo-cad
```


2. **Set up a Python Environment**
```bash
conda create --name paraturbocad python=3.9
conda activate paraturbocad
pip install -r requirements.txt
```

## Usage

The ParaturboCAD library is organized into four main classes, each responsible for different aspects of the turbocompressor design:
Classes and Their Functions

    IMPELLER(): Manages the centrifugal impeller design, from parameter extraction to 3D modeling.
    SGTB(): Tailored for designing and modeling Spiral Groove Thrust Bearings.
    ROTOR(): Handles the rotor design process, incorporating Herringbone Grooved Journal Bearings (HGJB).
    HELPER(): Provides utility tools to assist in the assembly of the turbocompressor.

## License

This project is licensed under the Apache License 2.0.

## Acknowledgments

    JustinSDK for cqMore
    CadQuery team and contributors for the CadQuery 2 library




