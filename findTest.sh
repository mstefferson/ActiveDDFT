#!/bin/bash/


a=`ls | grep 'poop' | grep '.mat'`
echo $a

if [ -z "$a" ];
then
  echo "no files"
  echo "files = $a"
else
  echo "found files"
  echo "files = $a"
fi

