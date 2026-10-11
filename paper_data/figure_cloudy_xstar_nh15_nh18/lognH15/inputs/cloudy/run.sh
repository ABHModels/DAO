#!/bin/zsh
set -eu
cd -- "${0:A:h}"
test ! -e cloudy.out
: "${CLOUDY_DATA_DIR:?Set this to your installed Cloudy data directory}"
export CLOUDY_DATA_PATH=".:$CLOUDY_DATA_DIR"
"${CLOUDY_EXECUTABLE:-cloudy.exe}" < cloudy.in > cloudy.out 2>&1
