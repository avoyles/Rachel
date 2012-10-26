#!/bin/bash

echo "ELAST_EXECUTABLE = `pwd`/elast_source/elast" > .rachel_setup
echo "GOSIA_EXECUTABLE = `pwd`/gosia_source/gosia_20110524.2" >> .rachel_setup
echo "BRICC_IDX_FILE = `pwd`/bricc_files/BrIccFOV22.idx" >> .rachel_setup
echo "BRICC_ICC_FILE = `pwd`/bricc_files/BrIccFOV22.icc" >> .rachel_setup
echo "RACHEL_DIRECTORY = `pwd`" >> .rachel_setup
echo "POPUPS = True" >> .rachel_setup
echo "POPUP_TIPS = True" >> .rachel_setup
cp -i .rachel_setup ~/

