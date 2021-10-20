#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
    stop("Need path to ipynb script.", call.=FALSE)
}

library(nbconvertR)

nbconvert(args[1], fmt="script")
