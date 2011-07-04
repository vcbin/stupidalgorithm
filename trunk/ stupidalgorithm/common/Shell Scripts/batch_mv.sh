#!/bin/bash - 
#===============================================================================
#
#          FILE:  batch_mv.sh
# 
#         USAGE:  ./batch_mv.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR: YOUR NAME (), 
#       COMPANY: 
#       CREATED: 2011年01月21日 23时04分05秒 CST
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

for file in `ls para*`
do
	new_name=`echo $file |sed -n 's/parameters_\(.*\.txt$\)/\1/gp'`
	mv -f $file $new_name
done
