#!/usr/bin/env bash


# $1 input fasta file
# $2 created db file
# $3 tmp folder for the search database indexing.

mmseqs createdb $1 $2

mmseqs createindex $2 $3