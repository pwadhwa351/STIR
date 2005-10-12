#! /bin/sh
# $Id$

# Script to normalise simset simulations
#
# With our current phg.rec files, Simset normalises the
# output always to the same mean counts, independent
# of the number of photons simulated.
# This script multiplies the results with the number of photons
# such that they correspond to simulations of an experiment with
# a given number of counts.

# This script relies on names of files fixed on other scripts, so beware.

# Authors:  Kris Thielemans


####################### Script inputs ##################################

if [ $# -eq 0 ]; then
    echo "usage: $0 directory_name ..."
    exit 1
fi

set -e

while [ $# -ne 0 ]; do
  pushd $1
  num_photons=`grep num_to_simulate phg.rec|awk '{print $4}'`
  for f in doubles11  doubles20  multiples  singles  trues; do
    stir_math -s --including-first --times-scalar $num_photons --divide-scalar 100000000 \
      ${f}_norm.hs ${f}.hs 2> /dev/null
  done
  popd
  shift
done

