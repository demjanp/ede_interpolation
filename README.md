# EDE Interpolation
A QGIS 3 plugin for spatio-temporal interpolation of archaeological settlement evidence.

Created on 21.4.2019

Evidence Density Estimation (EDE) interpolation of archaeological settlement data.

Produces spatio-temporal distribution maps or summed distributions representing intensity of settlement activities within the examined area in different time periods.

Uses radiocarbon-dated as well as typo-chronologically dated evidence as input.

Download: [https://plugins.qgis.org/plugins/ede_interpolation/](https://plugins.qgis.org/plugins/ede_interpolation/)

Wiki (also contains sample data): [https://osf.io/v7ahe/wiki/home/](https://osf.io/v7ahe/wiki/home/)

#### Input:
* Shapefile with point data representing evidence of settlement activities (e.g. archaeologically dated components of excavation sites) in a projected coordinate system.

#### Required fields for every feature in the source shapefile:
* `Spatial Accuracy (m)` - radius around registered point, where the actual location of the evidence is expected.
* `Dating Mean (years BP)` - mean value of dating of the archaeological component, representing either a Uniform Probability Distribution (UPD) of calendar years, or a Normal Probability Distribution (NPD) of radiocarbon years.
* `Dating Uncertainty (years)` - half length of the UPD interval or 1 standard deviation of the radiocarbon age in case of an NPD interval.
* `Dating Type ('UPD' or 'NPD')` - UPD is a range of calendar years BP (Before Present), assigned to an archaeological period (e.g. a culture), represented here as a mean value and half length of the interval. NPD is a radiocarbon age, represented here as a mean and standard deviation in radiocarbon years BP.

#### Model parameters:
* `Expected Settlement Duration (years)` - standard length of time that a settlement in the observed time and space is expected to exist, before it moves or in case of typological dating its cultural expression changes.
* `Expected Settlement Diameter (m)` - standard size of a settlement core. A diameter of 200 m means a settlement of cca 1 ha area.

#### Output parameters:
* `Time Step (years)` - distribution maps representing probability of presence of evidence of settlement activities at different spatial coordinates will be created at regular intervals in time, specified by this value. Evidence is summed along the time axis for each time step.
* `Time From (years BP)` - optional parameter specifying the begin of the observed time period.
* `Time To (years BP)` - optional parameter specifying the end of the observed time period.
* `Raster Cell Size (m)` - specifies length of the interpolation step along the spatial axes.
* `Approximate Spatial Probability (yes/no)` - set to use a faster, but slightly less precise calculation of the spatial component of the EDE function.
* `Save Output Layers to Directory` - specify directory where the spatial distribution maps in GeoTIFF format will be saved.
* `Save Summed Probability to File` - specify location of a CSV file, where values of probability summed along the spatial axis for every time step from the Intcal13 calibration curve (the highest possible temporal resolution) will be saved.

### Author:
Peter Demján [peter.demjan@gmail.com](peter.demjan@gmail.com)

Demján, P., & Dreslerová, D. (2016). Modelling distribution of archaeological settlement evidence based on heterogeneous spatial and temporal data. Journal of Archaeological Science, 69, 100–109. [DOI:10.1016/j.jas.2016.04.003](https://doi.org/10.1016/j.jas.2016.04.003)

## License: <a name="license"></a>

This code is licensed under the [GNU GENERAL PUBLIC LICENSE](https://www.gnu.org/licenses/gpl-3.0.en.html) - see the [LICENSE](LICENSE) file for details
