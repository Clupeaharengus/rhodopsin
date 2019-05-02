# rhodopsin
analysis of rhodopsin locus in Atlantic herring

Files for reproduction of results in ecological adaptation of Atlantic Herring rhodopsin

delta_allele_frequency (directory) - Allele frequency differences between Atlantic and Baltic populations
dAF_RH1.R   -   R script for measuring delta allele frquency and chi-square significance
rh1pools_BGI_s14.AD -   Allele depth per population table output from GATK for scaffold containing rhodopsin
rh1_samplelist.txt  -   Translation table between pool names and standard population designations
rh1_poolAF_secchi_sal_abs_20190204.csv  -   Table associating ecological variable with populations

pylogeny_habitat (directory)    -   Association between habitat and phylogeny
rabodend.R  -   R script to produce phylogeny associated with habitat with allelic identiy at 261
actinopt_12k_treePL.tre -   Time callibrated phylogeny from Rabosky 2018
Enbody_fish_species_habitat_combined_formatted.txt  -   Habitat data for species included in analysis
final_alignment.selection.pro.213_261.tsv   -   Protein alignment used in 261 analysis
spp_to_drop_badAln.txt  -   Species dropped due to lack of habitat data or bad alignment
timetree_trimmed.tre    -   time calibrated phylogeny with species defined above dropped
allele_freq.R   -   R script to test for significant associations between 261 and habitat
rabo_allele_hab.tsv -   Table of 261 allele, species, and habitat associations
newr213.R   -   R script to test for significant associations between 213 and habitat
final_alignment.translated.fullrhodopsin.fasta  -   Protein alignment of fish rhodopsin

ecological_variables_pops (directory)   -   Data and analysis used for associating populations with salinity, etc.
├── data
│ ├── 10m-admin-0-countries
│ │ ├── 10m_admin_0_countries.dbf
│ │ ├── 10m_admin_0_countries.prj
│ │ ├── 10m_admin_0_countries.shp
│ │ ├── 10m_admin_0_countries.shx
│ │ └── Icon\r
#General country outlines used in all map scripts
│ ├── environmental_variables
│ │ ├── Icon\r
│ │ ├── NASA_color
│ │ │ ├── 4Feb19_full_412_from_MODIS.tif
│ │ │ └── Icon\r
#Absorbance values from NASA MODIS sattelite data used in
Figure1_map_absorbance.R
│ │
└── salinity_ESA
│ │ ├── global-analysis-forecast-phy-001-024_1548930859285.nc
│ │ ├── global-analysis-forecast-phy-001-024_1548931193550.nc
│ │ └── Icon\r
#Salinity values from Europeon science commission (marine copernicus) database
for 2017 (averaged across all months)
│ ├── Icon\r
│ ├── pop_order_df.txt
#order of populations according to frequency of r261 in population. Used in all
scripts
│
└── rh1_poolAF_secchi_sal.csv
#Input file that includes coordinates that is used in all files for map generate
and extracting environmental variables.
├── Figure1_extract_salinity.R
#Extracts salinity from coordinates defined in data/rh1_poolAF_secchi_sal.csv and
data/environmental_variables/salinity_ESA
├── Figure1_map_absorbance.R
#Produces a map using absorbance values from data in
data/environmental_variables/NASA_color
├── Figure1_map_pi.R
#Using r261 frequences from data/rh1_poolAF_secchi_sal.csv produces pie plots to
be overlayed on top of the figure produced from Figure1_map_absorbance.R (done
manually in Adobe Illustrator)
├── Icon\r
└── output
├── Icon\r
├── piechart_plotting
│ ├── grob_pie_plots.pdf
│ ├── Icon\r
│ └── point.ref.pdf
├── rh1_poolAF_secchi_sal_20190327.csv
└── rh1_poolAF_secchi_sal_abs_20190327_100km.csv
#folder containing output files
