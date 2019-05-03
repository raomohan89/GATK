#!/bin/bash

if [ $# -ne 1 ]
then
    echo "Usage: ./evoquer.sh interval"
    exit 1
fi

./gatk Evoquer --project-id broad-dsp-spec-ops -L "$1" --sequence-dictionary ~/src/references/Homo_sapiens_assembly38_v0/Homo_sapiens_assembly38.dict --dataset-map evoquer_dataset_map -O evoquer.vcf
exit 0
