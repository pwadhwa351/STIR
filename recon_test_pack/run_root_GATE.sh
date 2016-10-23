#! /bin/sh
# Shell script for automatic running of the tests
# see README.txt
#
#  Copyright (C) 2016, UCL
#  This file is part of STIR.
#
#  This file is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 2.1 of the License, or
#  (at your option) any later version.

#  This file is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  See STIR/LICENSE.txt for details
#      
# Author Nikos Efthimiou

# Scripts should exit with error code when a test fails:
if [ -n "$TRAVIS" ]; then
    # The code runs inside Travis
    set -e
fi

echo This script should work with STIR version ">"3.0. If you have
echo a later version, you might have to update your test pack.
echo Please check the web site.
echo

#
# Parse option arguments (--)
# Note that the -- is required to suppress interpretation of $1 as options 
# to expr
#
while test `expr -- "$1" : "--.*"` -gt 0
do

if test "$1" = "--help"
  then
    echo "Usage: run_root_GATE.sh [install_dir]"
    echo "See README.txt for more info."
    exit 1
  else
    echo Warning: Unknown option "$1"
    echo rerun with --help for more info.
    exit 1
  fi

  shift 1

done 

${INSTALL_DIR}lm_to_projdata --input-formats > my_lm_supported_inputs.log 2>&1

if ! grep "ROOT" my_lm_supported_inputs.log > /dev/null
then
echo GATE support has not been installed in this system. Aborting.
exit 1;
else
echo "GATE support detected!"
fi

echo "Executing tests on ROOT files generated by GATE simulations, with cylindrical PET scanners"
echo ""
# first delete any files remaining from a previous run
rm -f my_*v my_*s my_*S

INSTALL_DIR=$1
ThereWereErrors=0
export INPUT_ROOT_FILE=test_PET_GATE.root
export INPUT=root_header.hroot
export TEMPLATE=template_for_ROOT_scanner.hs

echo "---- Converting ROOT file to ProjData file -----"
echo ">> Running lm_to_projdata for all events"
echo ""
export OUT_PROJDATA_FILE=my_proj_from_lm_all_events
export EXCLUDE_SCATTERED=0
export EXCLUDE_RANDOM=0
${INSTALL_DIR}lm_to_projdata ./lm_to_projdata.par 2>"./my_root_log_all.log"
all_events=$(awk -F ":" '/Number of prompts/ {print $2}' my_root_log_all.log)
echo "Number of prompts stored in this time period:" ${all_events}

echo ""

echo ">> Running lm_to_projdata *ONLY* for true events"
echo ""
export OUT_PROJDATA_FILE=my_proj_from_lm_true_events
export EXCLUDE_SCATTERED=1
export EXCLUDE_RANDOM=1
${INSTALL_DIR}lm_to_projdata ./lm_to_projdata.par 2>"./my_root_log_trues.log"
true_events=$(awk -F ":" '/Number of prompts/ {print $2}' my_root_log_trues.log)
echo "Number of prompts stored in this time period:" ${true_events}

echo ""
echo "Counting all values from ROOT file ..."
all_root_num=$(root -l ${INPUT_ROOT_FILE} << EOF | grep "Number of prompts stored in this time period" | grep -o -E '[0-9]+'
Coincidences->Draw(">>eventlist","","goff");
Int_t N = eventlist->GetN();
cout<<endl<<"Number of prompts stored in this time period:"<< N<<endl;
EOF)
echo "All events in ROOT file: "${all_root_num}
if [ "$all_events" -ne "$all_root_num" ]; then
echo "The total number of events used by STIR did not match the events in the root file. Error. "
exit 1
fi

echo ""
echo "Counting true values from ROOT file ..."
true_root_num=$(root -l ${INPUT_ROOT_FILE} << EOF | grep "Number" | grep -o -E '[0-9]+'
Coincidences->Draw(">>eventlist","eventID1 == eventID2 && comptonPhantom1 == 0 && comptonPhantom2 == 0","goff");
Int_t N = eventlist->GetN();
cout<<endl<<"Number of trues stored in this time period:"<< N<<endl;
EOF)

echo "True events in ROOT file :" ${true_root_num}
if [ "$true_events" -ne "$true_root_num" ]; then
echo "The total number of true events used by STIR did not match the events in the root file. Error. "
exit 1
fi

echo
echo '--------------- End of tests -------------'
echo
echo "Everything seems to be fine !"
echo 'You could remove all generated files using "rm -f my_* *.log"'
exit 0
