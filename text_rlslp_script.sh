#!/bin/bash

LOG_FILE="script.log"
> $LOG_FILE
exec > >(tee -a $LOG_FILE) 2>&1

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <input_slp> <lz77_path> [<rlslp_path>]"
    exit 1
fi

INPUT_FILE=$1
LZ77_FILE="${INPUT_FILE}.lz77"
SLG_FILE="${INPUT_FILE}.lz77.slg"
SLP_FILE="${INPUT_FILE}.lz77.slp"
PRUNE_SLP_FILE="${INPUT_FILE}.lz77.prune_slp"
RSLP_FILE="${INPUT_FILE}.rlslp"


if [ ! -f "$INPUT_FILE" ]; then
    echo "Input file not found: $INPUT_FILE"
    exit 1
fi

TEXT_TO_LZ77_DIR="slp-queries/text-to-lz77"
LZ77_TO_SLP_DIR="slp-queries/lz77-to-slp"
SLG_TO_SLP_DIR="slp-queries/slg-to-slp"
PRUNE_SLP_DIR="slp-queries/prune-slp"


echo "Building in $TEXT_TO_LZ77_DIR"
make -C "$TEXT_TO_LZ77_DIR" nuclear 
make -C "$TEXT_TO_LZ77_DIR" clean
make -C "$TEXT_TO_LZ77_DIR" 

echo "Running text_to_lz with input file: $INPUT_FILE"
yes | "$TEXT_TO_LZ77_DIR"/text_to_lz "$INPUT_FILE"

echo "Completed text_to_lz."

echo "Building in $LZ77_TO_SLP_DIR"
make -C "$LZ77_TO_SLP_DIR" nuclear
make -C "$LZ77_TO_SLP_DIR" clean
make -C "$LZ77_TO_SLP_DIR" 2>/dev/null

echo "Running convert with output file: $LZ77_FILE"
"$LZ77_TO_SLP_DIR"/lz_to_grammar "$LZ77_FILE"
echo "Completed convert."

echo "Building in $SLG_TO_SLP_DIR"

make -C "$SLG_TO_SLP_DIR" nuclear
make -C "$SLG_TO_SLP_DIR" clean
make -C "$SLG_TO_SLP_DIR" 2>/dev/null

echo "Running convert with SLG file: $SLG_FILE"
"$SLG_TO_SLP_DIR"/convert "$SLG_FILE" "$SLP_FILE" 

echo "Building in $PRUNE_SLP_DIR"

make -C "$PRUNE_SLP_DIR" nuclear
make -C "$PRUNE_SLP_DIR" clean
make -C "$PRUNE_SLP_DIR" 2>/dev/null


"$PRUNE_SLP_DIR"/prune "$SLP_FILE" "$PRUNE_SLP_FILE"


echo "Building in Recompression"

make nuclear
make clean
make 2>/dev/null

./recomp "$PRUNE_SLP_FILE" "-o" "$RSLP_FILE"

printf "\n"
echo "All Done!"
printf "\n"


awk '/peak =/ {
    gsub(/.*peak = /, "", $0);

    value = substr($0, 1, length($0) - 3); 
    unit = substr($0, length($0) - 2);

    if(value > max) {
        max = value;
    }
} 
/Peak RAM usage for Construction:/ {
    gsub(/.*Peak RAM usage for Construction: /, "", $0);
    
    value = substr($0, 1, length($0) - 3);
    unit = substr($0, length($0) - 2);

    if(value > max) {
        max = value;
    }
} 
END { 
    print "Peak RAM usage:", max, "MiB"; 
}' script.log

