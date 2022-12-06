#!/bin/bash
# File              : exeScript.sh
# Author            : Anton Riedel <anton.riedel@tum.de>
# Date              : 14.07.2021
# Last Modified Date: 25.11.2022
# Last Modified By  : Anton Riedel <anton.riedel@tum.de>

Index="0"

while read -r DataFile; do
    echo "Working on $DataFile"
    # root -l -q -b postProcessing.C+\(\"StandardCuts.json\",\"HistConfig.json\",\"$DataFile\"\,\"Output_StandardCuts_OUT-${Index}.root\"\)
    # root -l -q -b postProcessing.C+\(\"StandardCuts_NoPid.json\",\"HistConfig.json\",\"$DataFile\"\,\"Output_StandardCuts_NoPid_OUT-${Index}.root\"\)
    # root -l -q -b postProcessing.C+\(\"OpenCuts.json\",\"HistConfig.json\",\"$DataFile\"\,\"Output_OpenCuts_OUT-${Index}.root\"\)
    # root -l -q -b postProcessing.C+\(\"OpenCuts_NoPid.json\",\"HistConfig.json\",\"$DataFile\"\,\"Output_OpenCuts_NoPid_OUT-${Index}.root\"\)
    ./postProcessing.py "$DataFile" "Output_StandardCuts_OUT-${Index}.root" "HistConfig.json" "StandardCuts.json" &
    ./postProcessing.py "$DataFile" "Output_StandardCuts_NoPid_OUT-${Index}.root" "HistConfig.json" "StandardCuts_NoPid.json" &
    ./postProcessing.py "$DataFile" "Output_OpenCuts_OUT-${Index}.root" "HistConfig.json" "OpenCuts.json" &
    ./postProcessing.py "$DataFile" "Output_OpenCuts_NoPid_OUT-${Index}.root" "HistConfig.json" "OpenCuts_NoPid.json" &
    wait
    ((Index++))
done <$1

exit 0
