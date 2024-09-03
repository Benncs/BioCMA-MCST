#!/bin/bash 

type=debugoptmized
name="$1"
executable=./builddir/$type/apps/cli/biocma_mcst_cli_app

cli_args=$(python3 ./tools/cli_formater.py $name)

echo ./builddir/$type/apps/cli/biocma_mcst_cli_app $cli_args  

