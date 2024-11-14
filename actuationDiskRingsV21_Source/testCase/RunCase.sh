#TO BE RUN INSIDE caseRunned FOLDER DURING RunFullCase.sh

#---GET INTO CASE RUNNED FOLDER
cd caseRunned

#---CLEAN DIRECTORY AND COPY NEW FILES
rm -rf *
cp -rf ../problemMesh/* .
rm -rf 0/*
cp -rf 0.org.$BoundaryCond/* 0/.

#---EXPORT VARIABLES FOR PYTHON TO READ
export Uinf Uref TI D H UrefDiff

#---CALCULATE VARIABLES DEPENDING ON VELOCITY WITH PYTHON SCRIPT
#cd baseProblem
python3 calculatedSetup.py
#---RUN OpenFOAM
rm -rf log.decomposePar
rm -rf log.reconstructPar
rm -rf log.reconstructParMesh
if $tupac
then
  bash ./AllrunCaseParallelTUPAC
else
  bash ./AllrunCaseParallelLOCAL
fi

#---SAMPLE
postProcess -func sample -latestTime > log.sample
postProcess -func surfaces -latestTime > log.surfaces

#---PARAVIEW SCRIPT


#---SAVE FILES FOR LATER
if [[ $UrefDiff = 'no' ]]; then
  mkdir ../files_Uinf${Uinf}_TI${TI}_cD${cD}

  # cp -f log.simpleFoam ../files_Uinf${Uinf}_TI${TI}_cD${cD}/log.simpleFoam_Uinf${Uinf}_TI${TI}_cD${cD} # no está puesto que lo copie porque para casos de mallas muy finas es un archivo muy pesado y no da info que valga la pena
  cp -f outTurbines.csv ../files_Uinf${Uinf}_TI${TI}_cD${cD}/outTurbines_Uinf${Uinf}_TI${TI}_cD${cD}.txt
  cp -f outRings.csv ../files_Uinf${Uinf}_TI${TI}_cD${cD}/outRings_Uinf${Uinf}_TI${TI}_cD${cD}.txt
  if test -f outListFvOptions.txt; then
    cp -f outListFvOptions.txt ../files_Uinf${Uinf}_TI${TI}_cD${cD}/outListFvOptions_Uinf${Uinf}_TI${TI}_cD${cD}.txt
  fi
  cp -f log.checkSnappy ../files_Uinf${Uinf}_TI${TI}_cD${cD}/log.checkSnappy_Uinf${Uinf}_TI${TI}_cD${cD}
  cp -f log.simpleTime ../files_Uinf${Uinf}_TI${TI}_cD${cD}/log.simpleTime_Uinf${Uinf}_TI${TI}_cD${cD}
  cp -f includeSettings ../files_Uinf${Uinf}_TI${TI}_cD${cD}/includeSettings_Uinf${Uinf}_TI${TI}_cD${cD}
  cp -r system ../files_Uinf${Uinf}_TI${TI}_cD${cD}/system_Uinf${Uinf}_TI${TI}_cD${cD}


  mkdir ../files_Uinf${Uinf}_TI${TI}_cD${cD}/postProcessing_Uinf${Uinf}_TI${TI}_cD${cD}
  cp -rf postProcessing/*/*/* ../files_Uinf${Uinf}_TI${TI}_cD${cD}/postProcessing_Uinf${Uinf}_TI${TI}_cD${cD}/.
  cd ../files_Uinf${Uinf}_TI${TI}_cD${cD}/postProcessing_Uinf${Uinf}_TI${TI}_cD${cD}
  #-------RENAME POSTPROCESS FILES WITH PROBLEM MAIN DATA
  for FILENAME in *; do
    mv $FILENAME Uinf${Uinf}_TI${TI}_cD${cD}_$FILENAME;
  done
else
  mkdir ../files_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}

  # cp -f log.simpleFoam ../files_Uinf${Uinf}_TI${TI}_cD${cD}/log.simpleFoam_Uinf${Uinf}_TI${TI}_cD${cD} # no está puesto que lo copie porque para casos de mallas muy finas es un archivo muy pesado y no da info que valga la pena
  cp -f outTurbines.csv ../files_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}/outTurbines_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}.txt
  cp -f outRings.csv ../files_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}/outRings_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}.txt
  if test -f outListFvOptions.txt; then
    cp -f outListFvOptions.txt ../files_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}/outListFvOptions_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}.txt
  fi
  cp -f log.checkSnappy ../files_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}/log.checkSnappy_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}
  cp -f log.simpleTime ../files_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}/log.simpleTime_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}
  cp -f includeSettings ../files_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}/includeSettings_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}
  cp -r system ../files_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}/system_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}


  mkdir ../files_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}/postProcessing_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}
  cp -rf postProcessing/*/*/* ../files_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}/postProcessing_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}/.
  cd ../files_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}/postProcessing_Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}
  #-------RENAME POSTPROCESS FILES WITH PROBLEM MAIN DATA
  for FILENAME in *; do
    mv $FILENAME Uref${Uref}_Uinf${Uinf}_TI${TI}_cD${cD}_$FILENAME;
  done
fi
