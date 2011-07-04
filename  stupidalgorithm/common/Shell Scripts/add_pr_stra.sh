#!/bin/bash - 
#===============================================================================
#
#          FILE:  add_pr_stra.sh
# 
#         USAGE:  ./add_pr_stra.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR: YOUR NAME (), 
#       COMPANY: 
#       CREATED: 02/28/2011 10:24:19 AM CST
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

filelist=my_sde_para_file
sudo find ./ -regex ".*my_sde.*AP\.txt$" -print >$filelist
for file in `cat $filelist`
do
	sed -i "s/\(pr_stra =\).*/\1 2/g" $file
done

