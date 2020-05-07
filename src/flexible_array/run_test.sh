#!/bin/bash

export OMP_NUM_THREADS=1

valgrind \
	--leak-check=full --show-leak-kinds=all \
	./test
