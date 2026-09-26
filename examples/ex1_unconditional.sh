#!/usr/bin/env bash

pixi run -w PDchain rfd_chain 140-160 \
    --idealize --relax --ca_stdev 1 \
    --model_type protein_mpnn \
    --design_cycles 3 \
    --numdes 3 \
    --outprefix outputs/ex1_rand
