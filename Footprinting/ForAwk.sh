file=merp.txt; echo -e "merp\nprem\nmerp" > $file

awk -F"\n" '{ for (i=1; i<=3; i++) print $i }' $file