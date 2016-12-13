#!/bin/bash
# checking for missing programs.

if ! [ -x "$(command -v wget)" ]; then
  echo 'wget is not installed.' >&2
fi

if ! [ -x "$(command -v STAR)" ]; then
  echo 'STAR is not installed.' >&2
fi

if ! [ -x "$(command -v samtools)" ]; then
  echo 'samtools is not installed.' >&2
else
	if (( $(echo "`samtools --version-only | cut -d '+' -f1` < 1.3" | bc -l) )); then
		echo "samtools 1.3 or higher is needed."
	fi
fi

if ! [ -x "$(command -v java)" ]; then
  echo 'java is not installed.' >&2
fi

if ! [ -x "$(command -v bamtools)" ]; then
  echo 'bamtools is not installed.' >&2
fi

if ! [ -x "$(command -v bedtools)" ]; then
  echo 'bedtools is not installed.' >&2
fi

if ! [ -x "$(command -v awk)" ]; then
  echo 'awk is not installed.' >&2
fi

if ! [ -x "$(command -v perl)" ]; then
  echo 'perl is not installed.' >&2
fi

if ! [ -x "$(command -v gsnap)" ]; then
  echo 'gsnap is not installed.' >&2
fi

if ! [ -x "$(command -v liftOver)" ]; then
  echo 'liftOver is not installed.' >&2
fi

if ! [ -x "$(command -v parallel)" ]; then
  echo 'parallel is not installed.' >&2
fi

# check for bioperl
perl -e "{eval \"require Bio::Seq;\"; if (\$@) { print \"BioPerl is not installed.\n\"; }}"