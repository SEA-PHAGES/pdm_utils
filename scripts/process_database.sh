#!/bin/bash

DATE=$(date +%Y%m%d)
OUTPUT_FOLDER="${DATE}_phamerate"

DATABASE="Actino_Draft"
RUN_MACHINE="local"
OUTPUT_DIR=$1
CONFIG=$2
PEM_FILE=$3

echo "==========================================================="
echo "Updating database"
echo "==========================================================="

python3 update.py -d $DATABASE -o $OUTPUT_DIR -c $CONFIG -m $OUTPUT_FOLDER \
      --cpus 8 -e

echo "==========================================================="
echo "Phamerating database"
echo "==========================================================="

python3 phamerate.py -d $DATABASE -c $CONFIG \
      -o "${OUTPUT_DIR}/${OUTPUT_FOLDER}" \
      --cpus 8 \
      --run_machine "${RUN_MACHINE}"

echo "==========================================================="
echo "Uploading databases"
echo "==========================================================="


python3 upload.py -i "${OUTPUT_FOLDER}/${DATABASE}" -n $DATABASE \
      -c $CONFIG -k $PEM_FILE

# python3 upload.py -i "${OUTPUT_FOLDER}/${DATABASE}_v8" -n $DATABASE \
#      -c $CONFIG -k $PEM_FILE

python3 upload.py -i "${OUTPUT_FOLDER}/${DATABASE}_v6" -n $DATABASE \
      -c $CONFIG -k $PEM_FILE
