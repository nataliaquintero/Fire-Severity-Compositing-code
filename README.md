# Fire-Severity-Compositing-code

This repository contains Google Earth Engine (GEE) code for generating Landsat-based fire severity composites. The methods and data processing techniques are detailed in the paper titled *"Optimising Regional Fire Severity Mapping using Pixel-Based Image Compositing."* It includes scripts for various compositing methods and exports related GeoJSON and TIFF data.

## Overview

The repository includes the necessary scripts to generate fire severity maps using various compositing techniques:

- **Initial_Max-Min Composite**
- **Initial_Mean-Min Composite**
- **Extended_Medoid Composite**
- **Spring_Mean Composite**
- **Summer_Max-Min Composite**

These scripts allow the export of pixel centroids and TIFF imagery containing the following information:

- **Pre- and post-fire Normalized Burn Ratio (NBR)**
- **Near-Infrared (NIR) Band**
- **Short-Wave Infrared 2 (SWIR2) Band**
- **Pre- and post-fire Sensing Times**
- **Relativized Burned Ratio (RBR)**
- **Fire Year**

## How to Use

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/nataliaquintero/Fire-Severity-Compositing-code.git
   cd Fire-Severity-Compositing-code

2. #### Run the Scripts in Google Earth Engine:
Open the .js files in the Google Earth Engine (GEE) code editor.
Modify the study area, fire perimeters, and input parameters as needed.
Run the scripts to generate the fire severity composites and export the results to your Google Drive folder.

## Additional Information
For a detailed explanation of the methodology, data collection process, and the significance of each attribute, please refer to the associated research paper:

Paper Title: "Optimising Regional Fire Severity Mapping using Pixel-Based Image Compositing."
Authors: Natalia Quintero, Olga Viedma, Sander Veraverbeke, José Manuel Moreno
Affiliation: Universidad de Castilla La Mancha, Vrije Universiteit Amsterdam, University of East Anglia
The paper provides comprehensive insights into the study's objectives, data processing techniques, and the application of these codes.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Contact
For any questions or further information, please contact:

Natalia Quintero (natalia.quintero@uclm.es)
