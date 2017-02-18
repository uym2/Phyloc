# /bin/bash

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -i|--input)
    INPUT="$2"
    shift # past argument
    ;;
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    ;;
    -p|--rate)
    RATE="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

./3DAln $INPUT $OUTPUT $RATE

