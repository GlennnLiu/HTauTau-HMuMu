#!/bin/sh

for name in $1
do
mv ../backup/$name.h ../interface/
mv ../backup/$name.cc ../src/
done
