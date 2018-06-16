#!/bin/sh
module load python/2.7.8
module load abaqus/6.12

# If you only want to clean up all directories
# python paradis_abaq.py tests/frank_read_src_abaqus clean

 python paradis_abaq.py tests/frank_read_src_abaqus
