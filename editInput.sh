#!/bin/bash

# Check if a number was given as argument
if [ "$#" -ne 1 ]; then
        echo "Usage: ./edit.sh n"
            exit 1
            fi

# Check if the given argument is an integer
if ! [[ "$1" =~ ^[0-9]+$ ]]; then
        echo "Error: Please provide an integer value."
            exit 1
            fi

# Save the argument in a variable
n=$1

# Start writing to the file
echo "$((n))" > input.txt

# Loop to write the file names
for (( i=1; i<$n+1; i++ ))
        do
            #echo "h$i.mat" >> input.txt
            echo "NatOrbs-wholedensitymatrix-$i-1.mat" >> input.txt
        done
echo "File 'input.txt' has been updated."
