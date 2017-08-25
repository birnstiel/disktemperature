#!/bin/bash
n=`grep -n -m1 'END$' opacity.f | cut -f1 -d:`
sed -e "1,${n}d" < opacity.f > opacity_routine.f