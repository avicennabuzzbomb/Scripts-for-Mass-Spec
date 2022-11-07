#!/bin/bash
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# testing shell command line args
if [[ $# -gt 0 ]]; then 
    sleep 2
    echo "..."
    echo "BEEP BOOP input arguments is/are: "$1", is it supposed to mean something?"
else
    sleep 2
    echo "..."
    echo "Try inputing an argument like 'Hello World'"
fi