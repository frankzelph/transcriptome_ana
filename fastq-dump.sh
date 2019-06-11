#!/bin/bash

IFS=$'\n'
set -f
for ac in $(cat < "$1"); do
	echo $ac
	fastq-dump $ac
done

