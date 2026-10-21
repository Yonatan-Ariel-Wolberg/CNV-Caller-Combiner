#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/**
 * MODULE: common-merger.nf
 *
 * Unified abstraction layer for multi-caller SV consensus merging,
 * re-exporting both SURVIVOR and Truvari consensus routines.
 */

include { SURVIVOR; SURVIVOR_UNION } from './survivor.nf'
include { TRUVARI; TRUVARI_COLLAPSE } from './truvari.nf'
