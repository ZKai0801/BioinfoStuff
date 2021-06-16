#! /usr/bin/bash

# alias rm function with this script

TRANSH_DIR="$HOME/.trash"

if [[ ! -d $TRANSH_DIR  ]];
then
    mkdir $TRANSH_DIR
fi

for i in $*;
do
    stamp=`date +%y%m%d.%H%M`
    filename=`basename $i`

    mv $i $TRANSH_DIR/$stamp.$filename
done
