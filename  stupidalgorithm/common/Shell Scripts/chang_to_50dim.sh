#!/bin/bash - 
#===============================================================================
#
#          FILE:  chang_to_50dim.sh
# 
#         USAGE:  ./chang_to_50dim.sh 
# 
#   DESCRIPTION:  change AP parameters file to 50 dimension
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR: YOUR NAME (), 
#       COMPANY: 
#       CREATED: 2011年01月19日 23时28分01秒 CST
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

filelist=high_dim_AP_para_file
exclude_list="\(F14\|F15\|CamelBack\)"
sudo find ./ -regex ".*AP\.txt$" -print |grep -vG $exclude_list>$filelist
for file in `cat $filelist`
do
	sed -i '{s#\(run\s=\)\s60#\1 30#g;
	s#\(dimension\s=\)\s30#\1 50#g;
		s#\(stop_threshold\s=\)\s50000#\1 100000#g;
		s#\(learn_p\s=\)\s100#\1 200#g}' $file
done
