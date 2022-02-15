#!/bin/bash
#hostname

fa2gbk() {
	cat $1 | seqkit replace -s -p '(\w{10})' -r '$1 ' -w 66 | perl -ne 'if (/^>/) {print; $n=1}  else {s/ \r?\n$/\n/; printf "%9d %s", $n, $_; $n+=60;}'
}
