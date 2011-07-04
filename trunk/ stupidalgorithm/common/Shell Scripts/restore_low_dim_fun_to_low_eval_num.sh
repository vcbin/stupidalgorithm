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

target_list="\(F14\|F15\|CamelBack\)"
filelist=low_dim_AP_para_file
sudo find ./ -regex ".*AP\.txt$" -print |grep -G  $target_list>$filelist
for file in `cat $filelist`
do
	sed -i "{s#\(stop_threshold\s=\)\s$1#\1 $2#g;
		}" $file
done

