####################################################
##
## Project: amd_aging
##
## D.Dhirachaikulpanich@liverpool.ac.uk
## cyril.lagger@liverpool.ac.uk
##
####################################################
##

## Libraries ####

library(Seurat)
library(scDiffCom)
library(CellChat)
library(data.table)
library(ggplot2)
library(future)
library(ComplexUpset)
library(escape)
library(ggVennDiagram)
library(ggExtra)
library(limma)
library(ggpubr)

## Options ####

plan(multicore, workers = 20)

## Directory paths ####

path_data <- "data/"
path_results <- "results/"
