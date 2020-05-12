
# Given a directory, or pattern  of items with CONSISENT date timestamp in rought-to-fine format,
# return the item with latest timestamp.
# This is for use to find latest metadata file, image file or latest export genome directory
# Rough-to-fine format is longer time unit first, smaller time unit later.
# For example: 2020-05-11_12-26 is %Y-%m-%d_%H-%M
# Example usage:
#    findTheLatest "/srv/rbd/covid19/genomes/*export"
#    findTheLatest "/srv/rbd/thecontainer/*.sif"
#    It is IMPORTANT to quote the pattern
findTheLatest(){
    pattern="$1"

    for d in $(printf "%s\n" $pattern); do
        timestamp=$(echo ${d##*/}|sed -s 's/[^0-9]*//g');
        echo "$timestamp":::"$d"
    done|awk -F':::' '{print $1"\t"$2}'|sort -k1,1rn|head -n 1| awk -F'\t' '{print $2}'
}

