awk '{ if (NF > max)
	{
	max = NF;
	max_line_num= NR;} 
	}
      END { 
      		NR==max_line_num;
		print;
		print "\n" max_line_num ":" max "\n";
       	}' $1
