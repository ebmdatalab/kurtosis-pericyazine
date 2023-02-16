#!/bin/bash

### Modified figures generated for publication
### Journal: JMIR
### Includes:
### - resizing
### - removing transparency
### - removing white space
### - adding panel labels
### - combining panels vertically into one plot

# Figure 1

magick ratio_kurtosis_plot_with_histograms.png -resize 55% ratio_kurtosis_plot_with_histograms_resized.png

convert -flatten ratio_kurtosis_plot_with_histograms_resized.png ratio_kurtosis_plot_with_histograms_resized+noalpha.png

# PERICYAZINE - A

magick pericyazine_map.png -resize 55% pericyazine_map_resized.png

convert pericyazine_map_resized.png -gravity North -chop 0x100 pericyazine_map_resized2.png

convert -flatten pericyazine_map_resized2.png pericyazine_map_resized2+noalpha.png

convert pericyazine_map_resized2+noalpha.png -gravity NorthWest -pointsize 45 -annotate +0+0 '(A)' pericyazine_map_resized2+noalpha+labelled.png 


# PROMAZINE - B

magick promazine_map.png -resize 55% promazine_map_resized.png

convert promazine_map_resized.png -gravity North -chop 0x100 promazine_map_resized2.png

convert -flatten promazine_map_resized2.png promazine_map_resized2+noalpha.png

convert promazine_map_resized2+noalpha.png -gravity NorthWest -pointsize 45 -annotate +0+0 '(B)' promazine_map_resized2+noalpha+labelled.png 

# COMBINE
convert pericyazine_map_resized2+noalpha+labelled.png promazine_map_resized2+noalpha+labelled.png -append pericyazine+promazine_choropleth-maps.png

magick pericyazine+promazine_choropleth-maps.png -resize 80% pericyazine+promazine_choropleth-maps_TOC.png

# REMOVE intermediate plots
rm *_resized*.png

