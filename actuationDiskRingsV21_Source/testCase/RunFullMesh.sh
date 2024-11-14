#----------------------------------------
# created by Juan Ignacio Teich
# Facultad de Ingeniería, UBA
# jteich@fi.uba.ar
#----------------------------------------

#---Call the OpenFOAM functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

source baseSettings
for problemNum in $(seq $problemN1 1 $finalProblem | tr , .); do
  echo problemNumber: $problemNum
  #---READ problemSettings FILE.
  source problem${problemNum}Files/problem${problemNum}Settings

  #---REMOVE OLD PROBLEM DIRECTORY AND CREATE NEW ONE.
  #---CREATE problemMesh DIRECTORY
  #---COPY problemReadme FILE
  rm -rf problem$problemNumber
  mkdir problem$problemNumber
  cp -f problem${problemNumber}Files/problem${problemNumber}Readme.md problem$problemNumber/
  cd problem$problemNumber
  mkdir problemMesh
  cd problemMesh


  #---COPY BASE FILES
  cp -r ../../baseProblem/* .
  cp -r ../../problem${problemNumber}Files/fvOptions$problemNumber system/fvOptions
  cp -r ../../problem${problemNumber}Files/topoSetDict$problemNumber system/topoSetDict
  cp -r ../../problem${problemNumber}Files/snappyHexMeshDict$problemNumber system/snappyHexMeshDict
  cp -r ../../problem${problemNumber}Files/includeSettings$problemNumber includeSettings
  cp -r ../../problem${problemNumber}Files/sample$problemNumber system/sample
  cp -r ../../problem${problemNumber}Files/surfaces$problemNumber system/surfaces
  cp -r ../../problem${problemNumber}Files/problem${problemNumber}Settings system/problemSettings


  #---APPLY BOUNDARY CONDITIONS TO 0 FOLDER
  rm -rf 0.org
  cp -r 0.org.$BoundaryCond 0.org

  if [ $problemMeshNumber == 0 ]; then
      #---RUN OPENFOAM MESHING SCRIPT IN PARALLEL
      if $tupac
        then
          bash ./AllrunMeshParallelTUPAC
        else
          bash ./AllrunMeshParallelLOCAL
      fi
      #---REMOVE PROCESSORS DIRECTORIES
      rm -rf processor*
  else
      #---COPY DATA FROM SAME MESH IN OTHER PROBLEM
      echo Copying mesh from problem$problemMeshNumber
      cp -r ../../problem$problemMeshNumber/problemMesh/constant/polyMesh constant/.
      cp -r ../../problem$problemMeshNumber/problemMesh/log* .
      cp -r ../../problem$problemMeshNumber/problemMesh/0 .
      echo Done copying mesh
  fi


  #---END --> EXIT TO MAIN DIRECTORY
  cd ../..
done
