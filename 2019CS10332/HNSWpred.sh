#!/bin/bash
mpirun --bind-to none ./hnsw $1 $2 $3 $4
