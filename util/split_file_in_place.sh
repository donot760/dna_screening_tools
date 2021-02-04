#!/bin/zsh
# (c) whitequark 2010
# (c) dertalai 2015 (minimal modifications)

set -e

if [ $# != 2 ]; then
  echo "Usage: $0  "
  echo "  This script will split file to multiple parts, starting from"
  echo "  the end, and truncating the original file in process."
  echo "  Part size is specified in megabytes (1 MB = 1048576 bytes)."
  echo "  Use at your own risk."
  exit 0
fi

filename=$1
#partsize=$2
partsizeMB=$2
partsize=$(($2 * 1048576))

size=$(stat -c '%s' "${filename}")
parts=$(($size / $partsize))

do_split() {
  _part=$1
  _size=$2

  echo "Splitting part $_part"
  echo $(($partsize * ($_part - 1)))
  dd if="${filename}" of="${filename}.$(printf '%04d' $_part)" \
      count=$partsizeMB bs=1M skip=$((($_part - 1) * $partsizeMB))
  echo "Truncating source file"
  truncate "${filename}" --size="-$_size"
}

lastsize=$(($size % $partsize))
if [ $lastsize != 0 ]; then
  do_split $(($parts + 1)) $lastsize
fi

for i in $(seq $parts -1 1); do
  do_split $i $partsize
done

if [ $(stat -c '%s' "${filename}") == 0 ]; then
  rm "${filename}"
fi
