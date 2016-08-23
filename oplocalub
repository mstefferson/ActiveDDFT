# Bash executable
#!/bin/bash

# exe: bash oplocalub.sh #_files

# Number of files is an input, if none given do 1
if [ $# -eq 1 ];then
  nfiles=$1
else
nfiles=1
fi

echo "Starting run"
echo "Making OP objects for $nfilesi files"
echo "In dir `pwd` "
matlab -nodesktop -nosplash \
  -r  "try, OPHardRod( $nfiles ) , catch, exit(1), end, exit(0);" \
  2>&1 | tee opHR.out

echo "Finished. Matlab exit code: $?" 
exit
