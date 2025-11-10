<table>
  <tr>
    <td width="160" align="center" valign="top">
      <img src="https://github.com/user-attachments/assets/7e77bd5a-42d6-45db-b8e6-2c82cac81b9d" width="140" style="border-radius: 50%;"/>
    </td>
    <td valign="top">
      <b>For any inquiries or assistance, feel free to contact Mr. CAO Bin at:</b><br>
      📧 Email: <a href="mailto:bcao686@connect.hkust-gz.edu.cn">bcao686@connect.hkust-gz.edu.cn</a><br><br>
      Cao Bin is a PhD candidate at the <b>Hong Kong University of Science and Technology (Guangzhou)</b>, 
      under the supervision of Professor <a href="https://gbaaa.org.hk/en-us/article/67">Zhang Tong-Yi</a>. His research focuses on 
      <b>AI for science</b>, especially intelligent crystal-structure analysis and discovery. 
      Learn more about his work on his 
      <a href="https://www.caobin.asia/">homepage</a>.
    </td>
  </tr>
</table>


## Overview

XRetriever is a comprehensive tool for matching experimental X-ray diffraction (XRD) patterns to a database of crystal structures. It uses state-of-the-art algorithms to identify crystal phases from experimental XRD data.



## Project Structure

```
XRetriever/
├── process/                # Database processing and XRD calculation
│   ├── db_processor.py     # Crystal structure database reader
│   ├── xrd_calculator.py   # XRD spectrum calculator
│   └── build_database.py   # Database builder
├── retrieve/               # XRD matching and retrieval
│   ├── xrd_reader.py      # Experimental XRD file reader
│   ├── peak_detector.py   # Peak detection algorithms
│   ├── matcher.py         # Pattern matching algorithms
│   └── retriever.py       # Main Fun
├── API.ipynb              # Jupyter API
├── build_database_parallel.py # Database builder parallelly
├── requirements.txt       # Python dependencies

```

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Install Dependencies

```bash
pip install -r requirements.txt
```


### Required Packages

- numpy, scipy, pandas: Numerical computing
- ase: Atomic Simulation Environment for crystal structures
- pymatgen: Materials analysis, especially XRD calculations
- tqdm: Progress bars

## Quick Start

### 1. Build the XRD Database

see [tutorial](build_database.ipynb)

This processes the crystal structure database (db) and creates `xrd_database.pkl`. You need a input db file contains crystals (contact me for crystal data sharing if you needed)

I have generated a `xrd_database.pkl` file, which can be accessed at [huggingface](https://huggingface.co/datasets/caobin/MP_Xretriever).

### 2. Sample XRD Data 



Folder `exp_data` contains three experimental xrd. if you need more, please visit our website [Xqueryer](xqueryer.caobin.asia) or contact me.

### 3. Match XRD Pattern

see [tutorial](match_xrd_pattern.ipynb)

This matches the experimental pattern against the database and returns the top matches.

or use our [API](API.ipynb)

## Algorithm Details

### Peak Detection

- Baseline removal using morphological filtering
- Savitzky-Golay smoothing for noise reduction
- scipy.signal.find_peaks for robust peak detection
- Configurable height, prominence, and distance thresholds

### Pattern Matching

1. **Element Filtering**: Pre-filter database entries by chemical composition
2. **Peak Assignment**: Hungarian algorithm for optimal peak-to-peak matching
3. **Scoring**: Weighted combination of:
   - Position similarity (default: 70%)
   - Intensity similarity (default: 30%)
4. **Ranking**: Sort by combined score

### Tolerance Handling

- Position tolerance: ±0.1-0.3° in 2θ (configurable)
- Accounts for instrument calibration differences
- Handles peak shifts from sample displacement

## CSV File Format

Your experimental XRD file should be formatted as:

```csv
two_theta,intensity
10.00,5.2
10.02,5.5
10.04,5.1
...
```

The reader supports various delimiters (`,`, `\t`, ` `, `;`) and auto-detects the format.


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use XRetriever in your research, please cite:

```
@software{xretriever2025,
  title={XRetriever: XRD Phase Search Tool Based on Retrieval Strategy},
  author={Bin Cao},
  year={2025},
  url={https://github.com/Bin-Cao/XRetriever}
}
```
