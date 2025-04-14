#!/bin/bash
# ******************************************* #
# ------------------------------------------- #
# This is a script to run in command line GTS #
# ------------------------------------------- #
# ******************************************* #

__author__      = "Pierre Lesturgie"
__date__        = "2025-04-14"


## 1 - compute statistics for each tree contained in toy/toy.trees
gts stats -t toy/toy -Nanc 1000 -r 10 -R 2.5e-6 -D toy/scenario.txt -G toy/genome.txt

## 2 - summarize statistics
gts summary --tag toy/toy -r 10 -G toy/genome.txt 