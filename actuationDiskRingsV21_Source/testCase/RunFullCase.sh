#----------------------------------------
# created by Juan Ignacio Teich
# Facultad de Ingeniería, UBA
# jteich@fi.uba.ar
#----------------------------------------

#---Call the OpenFOAM functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

source baseSettings
for problemNum in $(seq $problemN1 1 $finalProblem | tr , .); do
  echo directorio "$PWD"
  echo problemNumber: $problemNum
  #---READ problemSettings FILE.
  source problem${problemNum}Files/problem${problemNum}Settings

  #---REMOVE OLD DIRECTORY AND CREATE caseRunned
  cp -r RunCase.sh problem$problemNumber
  cd problem$problemNumber
  rm -rf caseRunned
  mkdir caseRunned

  #---CHECK IF STEP IS USED: IF USED --> FOR, IF NOT USED --> GO ON TO NEXT POSSIBLE STEP VARIABLE
  if [ $UrefDiff == 'no' ]; then
    if [ $Ustep != 'off' ]; then
      for Uinf in $(seq $Umin $Ustep $Umax | tr , .); do
          if [ $TIstep != 'off' ]; then
            for TI in $(seq $TImin $TIstep $TImax | tr , .); do
              #-Ustep=on + ITstep=on-----------------------------
              Uref=$Uinf
              export Uinf Uref TI D H BoundaryCond cD UrefDiff tupac
              bash ./RunCase.sh
              echo Run with Uinf=$Uinf and TI=$TI finished
            done
          else
            #-Ustep=on + ITstep=off------------------------------
            Uref=$Uinf
            export Uinf Uref TI D H BoundaryCond cD UrefDiff tupac
            bash ./RunCase.sh
            echo Run with Uinf=$Uinf and TI=$TI finished
          fi
        done
    else
      if [ $TIstep != 'off' ]; then
        for TI in $(seq $TImin $TIstep $TImax | tr , .); do
          #-Ustep=off + ITstep=on--------------------------------
          Uref=$Uinf
          export Uinf Uref TI D H BoundaryCond cD UrefDiff tupac
          bash ./RunCase.sh
          echo Run with Uinf=$Uinf and TI=$TI finished
        done
      else
        #-Ustep=off + ITstep=off---------------------------------
        Uref=$Uinf
        export Uinf Uref TI D H BoundaryCond cD UrefDiff tupac
        bash ./RunCase.sh
        echo Run with Uinf=$Uinf and TI=$TI finished
      fi
    fi
  else
    if [ $UrefStep != 'off' ]; then
      for Uref in $(seq $UrefMin $UrefStep $UrefMax | tr , .); do
        if [ $Ustep != 'off' ]; then
          for Uinf in $(seq $Umin $Ustep $Umax | tr , .); do
              if [ $TIstep != 'off' ]; then
                for TI in $(seq $TImin $TIstep $TImax | tr , .); do
                  #-Ustep=on + ITstep=on-----------------------------
                  export Uinf Uref TI D H BoundaryCond cD UrefDiff
                  bash ./RunCase.sh
                  echo Run with Uinf=$Uinf and TI=$TI finished
                done
              else
                #-Ustep=on + ITstep=off------------------------------
                  export Uinf Uref TI D H BoundaryCond cD UrefDiff
                bash ./RunCase.sh
                echo Run with Uinf=$Uinf and TI=$TI finished
              fi
            done
        else
          if [ $TIstep != 'off' ]; then
            for TI in $(seq $TImin $TIstep $TImax | tr , .); do
              #-Ustep=off + ITstep=on--------------------------------
              export Uinf Uref TI D H BoundaryCond cD UrefDiff
              bash ./RunCase.sh
              echo Run with Uinf=$Uinf and TI=$TI finished
            done
          else
            #-Ustep=off + ITstep=off---------------------------------
            export Uinf Uref TI D H BoundaryCond cD UrefDiff
            bash ./RunCase.sh
            echo Run with Uinf=$Uinf and TI=$TI finished
          fi
        fi
        done
    else
      if [ $Ustep != 'off' ]; then
        for Uinf in $(seq $Umin $Ustep $Umax | tr , .); do
            if [ $TIstep != 'off' ]; then
              for TI in $(seq $TImin $TIstep $TImax | tr , .); do
                #-Ustep=on + ITstep=on-----------------------------
                export Uinf Uref TI D H BoundaryCond cD UrefDiff
                bash ./RunCase.sh
                echo Run with Uinf=$Uinf and TI=$TI finished
              done
            else
              #-Ustep=on + ITstep=off------------------------------
                export Uinf Uref TI D H BoundaryCond cD UrefDiff
              bash ./RunCase.sh
              echo Run with Uinf=$Uinf and TI=$TI finished
            fi
          done
      else
        if [ $TIstep != 'off' ]; then
          for TI in $(seq $TImin $TIstep $TImax | tr , .); do
            #-Ustep=off + ITstep=on--------------------------------
            export Uinf Uref TI D H BoundaryCond cD UrefDiff
            bash ./RunCase.sh
            echo Run with Uinf=$Uinf and TI=$TI finished
          done
        else
          #-Ustep=off + ITstep=off---------------------------------
          export Uinf Uref TI D H BoundaryCond cD UrefDiff
          bash ./RunCase.sh
          echo Run with Uinf=$Uinf and TI=$TI finished
        fi
      fi
    fi
  fi
  #---END --> EXIT TO MAIN DIRECTORY
  cd ..

done
