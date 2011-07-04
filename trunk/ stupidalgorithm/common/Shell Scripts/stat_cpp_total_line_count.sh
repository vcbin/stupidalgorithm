#@/bin/bash
IFS='
 
	'
PATH=/usr/local/bin:/usr/bin:/bin
export PATH

i=0;
i=expr $i +  $(xargs find $1 -regex '.*\.\(cpp\|h\)$' |wc -l|cut -f1);
echo i;
