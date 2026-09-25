#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { NORMALISE_CNV_QUALITY_SCORES; NORMALISE } from '../normalise/normalise-quality.nf'
