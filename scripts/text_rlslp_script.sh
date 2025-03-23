#!/bin/bash

print_usage() {
    echo "Usage: ./text_rlslp_script.sh <input_text> [block_size]"
    echo
    echo "<input_text>  : Path to the input text file (required)"
    echo "[block_size]  : Integer block size -- (text to lz77 -> BM Compression) (optional, default: 50)"
}

LOG_FILE="text_rlslp_script.log"
> $LOG_FILE
exec > >(tee -a $LOG_FILE) 2>&1

echo "Current timestamp: $(date)"

if [ "$#" -lt 1 ]; then
    print_usage
    exit 1
fi

INPUT_FILE=$1
BLOCK_SIZE=${2:-50}
LZ77_FILE="${INPUT_FILE}.lz77"
SLG_FILE="${INPUT_FILE}.lz77.slg"
SLP_FILE="${INPUT_FILE}.lz77.slp"
PRUNE_SLP_FILE="${INPUT_FILE}.lz77.prune_slp"
RSLP_FILE="${INPUT_FILE}.rlslp"


if [ ! -f "$INPUT_FILE" ]; then
    echo "Input file not found: $INPUT_FILE"
    exit 1
fi

TEXT_TO_LZ77_DIR="../code/slp-queries/text-to-lz77"
BM_COMPRESSION_DIR="../code/slp-queries/bm-text-to-lz77"
LZ77_TO_SLP_DIR="../code/slp-queries/lz77-to-slp"
SLG_TO_SLP_DIR="../code/slp-queries/slg-to-slp"
PRUNE_SLP_DIR="../code/slp-queries/prune-slp"

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
"$BM_COMPRESSION_DIR"/bm-compression "$INPUT_FILE" -o "$LZ77_FILE" "$BLOCK_SIZE"

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

make -C "../code/" nuclear
make -C "../code/" clean
make -C "../code/" 2>/dev/null

../code/recomp "$PRUNE_SLP_FILE" "-o" "$RSLP_FILE" # "-t" "$INPUT_FILE"

printf "\n"
echo "All Done!"
printf "\n"

echo
echo "--------------- SUMMARY -----------------"
echo

awk '
BEGIN {
    stages[1] = "text_lz77"
    stages[2] = "lz77_slg"
    stages[3] = "slg_slp"
    stages[4] = "slp_prunedslp"
    stages[5] = "prunedslp_rlslp"

    pattern_map["Block Size:"] = "block_size"
    pattern_map["text length:"] = "text_length"
    pattern_map["Parsing Size:"] = "lz77_phrases"
    pattern_map["Number of nonterminals ="] = "slg_size"
    pattern_map["Size of original SLP ="] = "slp_size"
    pattern_map["Size of pruned SLP ="] = "prunedslp_size"
    pattern_map["Productions ="] = "rlslp_size"

    time_index = 1
    ram_index = 1
}



/peak =|Peak RAM usage for Construction/ {
    # print $0
    for(i = 1; i <= NF; i++) {
        if($i ~ /[0-9.]+MiB/) {
            # printf $i
            # printf "\n"

            value = substr($i, 1, length($i)-3) + 0 # convert to numeric
        }
        else if(i < NF && $i ~ /^[0-9.]+$/ && $(i+1) == "MB") {
            # printf $i
            # printf "\n"
            
            value = $i + 0 # convert to numeric
        }
        else {
            continue
        }

        if(ram_index <= 5) {
            stage = stages[ram_index]
            rams[stage] = value
            if(value > peakRAM) {
                peakRAM = value
            }
            ram_index++
        }
        break

    }
} 

# By default awk splits the line by spaces. NF = number of fields/tokens/words
/Compute SA|Compute LZ77|Conversion time|Read SLG from file and convert to SLP|Read SLP from file|time =|Time taken for Construction/ {
    # print $0
    for(i = 1; i <= NF; i++) {
        gsub(/[()=]/, "", $i)
        # printf $i, "\n"
        if($i ~ /[0-9.]+s/) {
            # printf $i
            # printf "\n"

            time = substr($i, 1, length($i)-1)
            # totalTime += time
            # printf "%.3f\n", time
        }
        else if(i < NF && $i ~ /^[0-9.]+$/ && $(i+1) == "seconds") {
            # printf $i
            # printf "\n"

            time = substr($i, 1, length($i)-1)
            # totalTime += time
            # printf "%.3f\n", time
        }
        else {
            continue
        }

        # printf time, "\n"
        if(time_index <= 5) {
            stage = stages[time_index]
            times[stage] = time
            # printf time_index, times[stage], "\n"
            totalTime += time
            time_index++
        }

        break
    }
}

{
    for(p in pattern_map) {
        # printf p, "\n"
        if(index($0, p) && !seen[pattern_map[p]]) {
            for(i = 1; i <= NF; i++) {
                if($i ~ /^[0-9]+$/) {
                    label = pattern_map[p]
                    if(label ~ /block_size|text_length/) {
                        printf "%s = %d", label, $i + 0
                        printf "\n"
                    }
                    else {
                        printf "%s = %d", label, $i + 0
                        printf "\n"
                    }
                    seen[label] = 1
                    break
                }
            }
        }
    }
}
END {
    printf "\n"
    for(i = 1; i <= 5; i++) {
        stage = stages[i]
        printf "%s_time = %.2fs\n", stage, times[stage]
        printf "%s_ram = %.2fMiB\n", stage, rams[stage]
    }

    printf "\n"
    printf "total_time = %.2fs\n", totalTime
    printf "overall_peak = %.2fMiB\n", peakRAM
}' "$LOG_FILE"
