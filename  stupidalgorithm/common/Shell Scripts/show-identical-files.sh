#!/bin/bash
IFS='
 
	'
PATH=/usr/local/bin:/usr/bin:/bin
export PATH

md5sum "$@" /dev/null 2> /dev/null |
awk '{
	count[$1]++
	if (count[$1] == 1) first[$1] = $0
	if (count[$1] == 2) print first[$1]
	if (count[$1] > 1) print $0
}' |
sort |
	awk '{
	if (last != $1) print ""
	last = $1
	print
}'
