#! /usr/bin/bash

# alias rm function with this script

TRANSH_DIR="$HOME/.trash"

if [[ ! -d $TRANSH_DIR  ]];
then
    mkdir $TRANSH_DIR
    touch $TRANSH_DIR/filepaths
fi

for i in $*;
do
    stamp=`date +%y%m%d.%H%M`
    filename=`basename $i`
    
    echo "$stamp.$filename $(realpath $i)" >> $TRANSH_DIR/filepaths
    mv $i $TRANSH_DIR/$stamp.$filename
done
