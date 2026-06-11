#!/bin/bash

DATE=$(date +%Y%m%d)
OUTPUT_FOLDER="${DATE}_phamerate"

RUN_MACHINE="local" # local or remote
OUTPUT_DIR=$1
CONFIG=$2
PEM_FILE=$3
DATABASE=$4
NAME=$5

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

echo "database: $DATABASE"
echo "name: $NAME"

echo "==========================================================="
echo "Updating database"
echo "==========================================================="

python3 $SCRIPT_DIR/update.py -d $DATABASE -o $OUTPUT_DIR -c $CONFIG -m $OUTPUT_FOLDER \
    --cpus 8 -e

echo "==========================================================="
echo "Phamerating database"
echo "==========================================================="

python3 $SCRIPT_DIR/phamerate.py -d $DATABASE -n $NAME -c $CONFIG \
      -o "${OUTPUT_DIR}/${OUTPUT_FOLDER}" \
      --cpus 8 \
      --run_machine "${RUN_MACHINE}"

echo "==========================================================="
echo "Uploading databases"
echo "==========================================================="

python3 $SCRIPT_DIR/upload.py -i "${OUTPUT_FOLDER}/${DATABASE}" -n $NAME \
     -c $CONFIG -k $PEM_FILE

python3 $SCRIPT_DIR/upload.py -i "${OUTPUT_FOLDER}/${DATABASE}_v6" -n $DATABASE \
     -c $CONFIG -k $PEM_FILE
