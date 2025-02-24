# PDB Cleaner

## Overview
The **PDB Cleaner** script is designed to clean and process multiple **PDB (Protein Data Bank)** files efficiently. It allows users to perform various operations such as removing water molecules, handling alternate locations, and reporting sequence gaps. The cleaned PDB structures are then saved in a specified output directory.

This script can be used to clean many structures at a time. In my use case, I had a folder containing multiple PDB structures, and I processed them to get cleaned output files into a different folder.

## Features & Actions Performed
This script performs the following actions:
- **Removes water molecules** (`-rmw` flag)
- **Keeps or removes hydrogen atoms** (`-kh` flag)
- **Handles alternate locations** by keeping the highest occupancy (`-alt` flag)
- **Removes insertion codes** (`-ins` flag)
- **Reports sequence gaps** in the structure (`-gaps` flag)
- **Ensures negative occupancy values are set to 0.00**
- **Extracts cleaned atoms and saves them in a new PDB file**

## Requirements
The script requires the following dependencies:
- Python 3.x
- `numpy`
- `pandas`
- `biopython`

To install the dependencies, run:
```bash
pip install numpy pandas biopython
```

## Usage
Run the script with default parameters:
```bash
python clean.py input_pdbs output_pdbs
```

### Optional Arguments
| Flag | Description |
|------|-------------|
| `-rmw, --remove-waters` | Removes water molecules from PDB files |
| `-kh, --keep-hydrogens` | Keeps hydrogen atoms in the structure |
| `-alt, --handle-altloc` | Retains only the highest occupancy alternate locations |
| `-ins, --remove-insertions` | Removes insertion codes from residues |
| `-gaps, --report-gaps` | Identifies and reports sequence gaps |

### Example Usage
To clean PDB files by removing water molecules and reporting sequence gaps:
```bash
python clean.py input_pdbs output_pdbs -rmw -gaps
```

## Output
- The cleaned PDB files will be saved in the specified `output_pdbs` directory.
- If the `-gaps` flag is used, sequence gap details will be printed in the console.

## License
This project is open-source and available under the MIT License.

## Contributions
Contributions are welcome! Feel free to submit issues or pull requests to improve the script.

