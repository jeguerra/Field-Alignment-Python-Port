#!/bin/bash
	rm -rf foo_*
	rm -rf out*
	rm -rf memFA20mem_*
	rm -rf runqsub_*
	rm -rf memFA0.e*
	rm -rf memFA0.o*
	mv copy_foo.txt foo.txt
	mv copy_memFA20mem.m memFA20mem.m
	mv copy_runqsub.sh runqsub.sh
	rm -rf core.3*
