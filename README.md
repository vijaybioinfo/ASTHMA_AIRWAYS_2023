# ASTHMA AIRWAYS 2021

This repository contains the scripts used to analyzed our asthma airways bulk and single-cell data.

Developer: Ciro Ramírez-Suástegui (ciro@lji.org)

Vijayanand Lab
Division of Vaccine Discovery
La Jolla Institute for Immunology
La Jolla, CA 92037, USA

## Sex stratified associations of airway cytotoxic CD4+ T cells with tissue-residency and pro-inflammatory features in severe asthma

Description
---

*Demultiplexing libraries*: Cell Ranger was used to demultiplex the 10x libraries.
*Demultiplexing donors*: Demuxlet was used to derive the cells' donor source.
*Quality control*: [An in-house script](https://github.com/vijaybioinfo/quality_control/blob/main/single_cell.R) was used to explore the quality of the data and select the thresholds.
Clustering: Seurat was used to cluster the data.

For more specific information about the data generation and processing, please check citation.

Citation
---
Please cite the following manuscript if you are using this repository:

*Under review*

Contact
---
Please email Ciro Ramírez-Suástegui (ciro@lji.org) and Vijayanand Pandurangan (vijay@lji.org).
