#!/usr/bin/bash

hmmfile="$1"
hmmname="$2"

hmmfetch $1 $2 >"$2".hmm
hmmemit "$2".hmm >"$2".faa
cat "$2".faa
