#!/bin/bash

LOG_FILE="text_rlslp_script.log"
> $LOG_FILE
exec > >(tee -a $LOG_FILE) 2>&1

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <input_text>"
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
BM_COMPRESSION_DIR="slp-queries/bm-text-to-lz77"
LZ77_TO_SLP_DIR="slp-queries/lz77-to-slp"
SLG_TO_SLP_DIR="slp-queries/slg-to-slp"
PRUNE_SLP_DIR="slp-queries/prune-slp"

# echo ""
# echo "Building in $TEXT_TO_LZ77_DIR"
# make -C "$TEXT_TO_LZ77_DIR" nuclear 
# make -C "$TEXT_TO_LZ77_DIR" clean
# make -C "$TEXT_TO_LZ77_DIR" 

echo ""
echo "Building in $BM_COMPRESSION_DIR"
make -C "$BM_COMPRESSION_DIR" nuclear 
make -C "$BM_COMPRESSION_DIR" clean
make -C "$BM_COMPRESSION_DIR" 

echo "Running bm_text_to_lz with input file: $INPUT_FILE"
# yes | "$TEXT_TO_LZ77_DIR"/text_to_lz "$INPUT_FILE"
"$BM_COMPRESSION_DIR"/bm-compression "$INPUT_FILE" -o "$LZ77_FILE"

echo "Completed bm_text_to_lz."

echo ""
echo "Building in $LZ77_TO_SLP_DIR"
make -C "$LZ77_TO_SLP_DIR" nuclear
make -C "$LZ77_TO_SLP_DIR" clean
make -C "$LZ77_TO_SLP_DIR" 2>/dev/null

echo "Running convert with output file: $LZ77_FILE"
"$LZ77_TO_SLP_DIR"/lz_to_grammar "$LZ77_FILE"
echo "Completed convert."

echo ""
echo "Building in $SLG_TO_SLP_DIR"

make -C "$SLG_TO_SLP_DIR" nuclear
make -C "$SLG_TO_SLP_DIR" clean
make -C "$SLG_TO_SLP_DIR" 2>/dev/null

echo "Running convert with SLG file: $SLG_FILE"
"$SLG_TO_SLP_DIR"/convert "$SLG_FILE" "$SLP_FILE" 

echo ""
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


awk '/peak =|Peak RAM usage for Construction/ {
    # print $0
    for(i = 1; i <= NF; i++) {
        if($i ~ /[0-9.]+MiB/) {
            # printf $i
            # printf "\n"

            value = substr($i, 1, length($i)-3) + 0 # convert to numeric
            if(res < value) {
                res = value
            }
        }
        else if(i < NF && $i ~ /^[0-9.]+$/ && $(i+1) == "MB") {
            # printf $i
            # printf "\n"
            
            value = $i + 0 # convert to numeric
            if(res < value) {
                res = value
            }
        }
    }
} 
END { 
    print "*** Overall Peak RAM usage:", res, "MiB ***"; 
}' $LOG_FILE

awk '
# By default awk splits the line by spaces. NF = number of fields/tokens/words
/Compute SA|Compute LZ77|Conversion time|Read SLG from file and convert to SLP|Read SLP from file|time =|Time taken for Construction/ {
    # print $0
    for(i = 1; i <= NF; i++) {
        if($i ~ /[0-9.]+s/) {
            # printf $i
            # printf "\n"

            time = substr($i, 1, length($i)-1)
            totalTime += time
            # printf "%.3f\n", time
        }
        else if(i < NF && $i ~ /^[0-9.]+$/ && $(i+1) == "seconds") {
            # printf $i
            # printf "\n"

            time = substr($i, 1, length($i)-1)
            totalTime += time
            # printf "%.3f\n", time
        }
    }
}
END {
    printf "*** Total time to convert Text to RLSLP: %.3fs ***\n", totalTime
}' $LOG_FILE
