#----------------------------------------
# created by Juan Ignacio Teich
# Facultad de Ingeniería, UBA
# jteich@fi.uba.ar
#----------------------------------------

#---Call the OpenFOAM functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

#---RUN MESHING
rm -rf log.RunFullMesh
echo running RunFullMesh
bash RunFullMesh.sh > log.RunFullMesh

#---RUN CASE
rm -rf log.RunFullCase
echo running RunFullCase
bash RunFullCase.sh > log.RunFullCase
