#!/bin/bash - 

set -o nounset                              # Treat unset variables as an error

awk '
{
 if ( $0 ~/^#.*/ || $0 ~/^\s*$/ )
   next;
 print $1,$2;
}' $1 >obj_coor

