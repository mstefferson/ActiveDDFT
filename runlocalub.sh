# Bash executable
#!/bin/bash

echo "Starting run"
echo "In dir `pwd` "
matlab -nodesktop -nosplash \
  -r  "try, RunHardRod, catch, exit(1), end, exit(0);" \
  2>&1 | tee runHR.out

echo "Finished. Matlab exit code: $?" 
exit
