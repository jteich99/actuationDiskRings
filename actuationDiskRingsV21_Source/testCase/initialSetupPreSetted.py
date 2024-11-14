#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
----------------------------------------
created by Juan Ignacio Teich
facultad de Ingeniería, UBA
jteich@fi.uba.ar
----------------------------------------
"""

"""
version 2023-10-22
    incorporates sample creation from initalSettings file
"""

def is_float(string):
    try:
        # float() is a built-in function
        float(string)
        return True
    except ValueError:
        return False

#import csv
import os
from shutil import copy
import math
from initialSettings import *
#from pathlib import Pathº


#---ASK FOR AMOUNT OF PROBLEMS TO SETUP
print(f"Cantidad de problemas a setupear = {len(problems)}")
# print("Numero inicial de problema = " + str(problemN1))

finalProblem= len(problems) + 1
baseSettings=open('baseSettings','w')
baseSettings.write(
    f'problemN1=1\n' +
    f'totalProblems={len(problems)}\n' +
    f'finalProblem={len(problems)}\n'
)
if tupac:
    baseSettings.write('tupac=true\n')
else:
    baseSettings.write('tupac=false\n')


#---ASK FOR PROCESSOR PARAMETERS. THEY REMAIN CONSTANT IN ALL PROBLEMS
# print("\nCantidad de nodos (si computadora personal-->1) = " + str(nodes))
# print("Cantidad total de cores = " + str(cores))

#---for THROUGH ALL PROBLEM NUMBERS
for problemNumber in range(len(problems)):
    problemDict = problems[problemNumber]
    problemNumber += 1
    print('\nSETUP PROBLEM ' + str(problemNumber) + '---------------------\n')
    #---OPEN FILES AND APPEND NUMBER TO FILE NAME
    if os.path.exists(f'problem{problemNumber}Files')!=True:
        os.makedirs(f'problem{problemNumber}Files')
    includeSettings=f'problem{problemNumber}Files/includeSettings{problemNumber}'
    includeSettingsFile = open(includeSettings ,"w")
    problemSettings=f'problem{problemNumber}Files/problem{problemNumber}Settings'
    problemSettingsFile = open(problemSettings, "w")
    
    includeSettingsFile.write(f'//----PROBLEM {problemNumber} CONFIG FILE---------\n\n')
    problemSettingsFile.write(f'problemNumber={problemNumber}\n')
    
    #---SETUP OF README
    src_path = 'problemReadme.md'
    destination_path = f'problem{problemNumber}Files/problem{problemNumber}Readme.md'
    if os.path.isfile(destination_path)==True:
        os.remove(destination_path)
    copy(src_path, destination_path) 
    readme=open(destination_path,'a')
    
    #---BASIC SETUP OF problemSettings
    #---SETUP OF PROCESSOR PARAMETERS
    problemSettingsFile.write("\n# NUMBER OF NODES & NUMBER OF CORES (TOTAL)\n")
    problemSettingsFile.write(f"nodes={problemDict['nodes']}\n")
    problemSettingsFile.write(f"cores={problemDict['cores']}\n")
    includeSettingsFile.write(f"coresTotal {problemDict['cores']};\n")
    
    #---PROBLEM SETUP
    #---SETUP OF DIAMETER. SAME FOR ALL ACTUATORS
    # D = globals()['D%s' % problemNumber] # --> THIS RESULTS IN D = THE VARIABLE IN initalSettings THAT IS NAMED D#problemNumber
    diameter = problemDict["diameter"]
    print(f"\nDiametro de actuadores = {diameter} mts")
    includeSettingsFile.write(
        "\n// Diameter \n" +
        f"D {diameter};\n")
    problemMeshNr = problemDict["problemMeshNumber"]
    problemSettingsFile.write(
        "\n# problem Mesh Number\n" + 
        "# if != 0 then it will copy the mesh from the problem with the number specified\n" + 
        f"problemMeshNumber={problemMeshNr}\n")
    problemSettingsFile.write(f"\n# actuator diameter\nD={diameter}\n")
    readme.write(f"\nActuator diameter - D={diameter}\n")
    
    #---DOMAIN SIZE
    domainSize = problemDict["domain"]
    print('\nTamaño del dominio.')
    Dx = domainSize[0] 
    print(f"Cantidad de diametros aguas abajo = {Dx}")
    Lx=Dx*diameter
    
    Dy = domainSize[1] 
    print(f"Cantidad de diametros de ancho de tunel = {Dy}")
    Ly=Dy*diameter
    
    Dz = domainSize[2] 
    print(f"Cantidad de diametros de alto de tunel = {Dz}")
    Lz=Dz*diameter
    
    readme.write(f"\nDomain size in diameters per axis: {Dx}x{Dy}x{Dz}\n")
    
    #---MESHING CORNERS CALCULATION
    includeSettingsFile.write("// domain corners x=x1, y=x2, z=x3----------\n")
    for i in range(0,8):
        if i==0 or i==3 or i==4 or i==7:
            x = 0
        elif i==1 or i==2 or i==5 or i==6:
            x = Lx
        if i==0 or i==1 or i==4 or i==5:
            y = 0
        elif i==2 or i==3 or i==6 or i==7:
            y = Ly
        if i==0 or i==1 or i==2 or i==3:
            z = 0
        elif i==4 or i==5 or i==6 or i==7:
            z = Lz
        includeSettingsFile.write(
            f'x_{i} {x};\n' +
            f'y_{i} {y};\n' +
            f'z_{i} {z};\n' )
    includeSettingsFile.write("\n")
    
    #---MESHING SETUP FOR BLOCKMESH
    c_D=problemDict["cellsPerDiameter"]
    print(f"Cantidad de celdas por diámetro, en zona menos densa = {c_D}")
    readme.write(f"\nCells per diameter: {c_D}\n")
    problemSettingsFile.write(f"\n# Cells per diameter:\ncD={c_D}\n")
    l_celda = diameter/c_D
    print(f"Longitud de celda: {l_celda}")
    includeSettingsFile.write("// cells per axis (c_x c_y c_z)----------\n")
    c_x = Dx*c_D
    c_y = Dy*c_D
    c_z = Dz*c_D
    includeSettingsFile.write(
        f"c_x {c_x} ;\n" +
        f"c_y {c_y} ;\n" +
        f"c_z {c_z} ;\n" +
        "\n"
    )
    
    #---BOUNDARY CONDITION
    print("Condicion de contorno a utilizar [ABL 'logaritmica' o SLIP 'uniforme']:")
    BC = problemDict["boundaryCondition"]

    if BC == 'ABL':
        height= problemDict["nacelleHeight"]
        ce_z=5
    elif BC == 'SLIP':
        height= Lz/2 
        ce_z=1

    print('Se considera actuador a H=' + str(height))
    includeSettingsFile.write(
        "// cell expansion rate in each axis (ce_x ce_y ce_z)\n" +
        "ce_x 1;\n" +
        "ce_y 1;\n" +
        f"ce_z {ce_z};\n\n" +
        f"// turbine height\n" +
        f"H {height};\n\n" )
    
    problemSettingsFile.write(
        "\n# boundary condition to apply\n" + 
        f"BoundaryCond={BC}\n" +
        "\n# turbine height \n" + 
        f"H={height}\n")

    readme.write(
        f"\nBoundary condition: {BC}\n" +
        f"Height of actuator [mts]: {height}\n")
    
    #---AMOUNT OF ACTUATORS IN PROBLEM
    amountActuators = len(problemDict["actuators"]) + 1
    readme.write(f'\nNumber of actuators in problem: {amountActuators}')
    
    #---COPY BASE fvOptions FILE AND RENAME IT ACCORDING TO PROBLEM NUMBER
    fvOptionsfilename=f'problem{problemNumber}Files/fvOptions{problemNumber}'
    copy('baseFiles/fvOptions.base', fvOptionsfilename) 
    fOptfile=open(fvOptionsfilename, 'a+')
    
    #---COPY BASE topoSetDict FILE AND RENAME IT ACCORDING TO PROBLEM NUMBER
    tSDfilename=f'problem{problemNumber}Files/topoSetDict{problemNumber}'
    copy('baseFiles/topoSetDict.base', tSDfilename) 
    tSDfile=open(tSDfilename, 'a+')
    
    #---COPY BASE snappyHexMeshDict FILE AND RENAME IT ACCORDING TO PROBLEM NUMBER
    sHMDfilename=f'problem{problemNumber}Files/snappyHexMeshDict{problemNumber}'
    copy('baseFiles/snappyHexMeshDict.base', sHMDfilename) 
    sHMDfile=open(sHMDfilename, 'a+')
    sHMDaux=""
    
    #---LOOP TO SET UP EACH ACTUATOR
    for actNum in range(1,amountActuators):
        print(f"\nSetup actuador numero {actNum}---------\n")
        actuatorInfo = problemDict["actuators"][actNum-1]
        actuatorModel = actuatorInfo["ADmodel"]
        
# A SER USADO CUANDO TENGA MAS DE UN TIPO DE ACTUADOR---------------------------
        #---ACTUATOR TYPE. SELECT FROM LIST.
        # actuatorOptions = ['actuationDiskRingsV11']
        # print("Tipo de actuador:")
        # contActOpt=0
        # for item in actuatorOptions:
        #     contActOpt+=1
        #     print(str(contActOpt) + '. ' + str(item))
        # actType=input('Eleccion: ')
        # while actType not in range(1,len(actuatorOptions)+1):
        #     actType=input('Elegir numero de la lista: ')
        # actuatorType=actuatorOptions[int(actType)-1]
#-------------------------------------------------------------------------------
        #ACTUATOR OPTIONS ARE
        #   - actuationDiskRingsV11_calib_Source
        #   - actuationDiskRingsV11_Source
        #   - actuationDiskRingsV2_Source
        #   - actuationDiskRingsV2_calib_Source
        if actuatorModel == 99:
            actuatorType = "actuationDiskRingsV2_Source"
        else:
            actuatorType = "actuationDiskRingsV21_Source"

        includeSettingsFile.write(f'//----------actuador nro {actNum} ----------\n')
        readme.write(f'- actuator n{actNum} setup:\n')
        
        #---DETERMINE POSITION OF ACTUATOR. IF NOT FIRST, POSITION WILL BE RELATIVE TO 1ST
        x_act = actuatorInfo["xPosition"]

        if actNum == 1:
            if actuatorInfo["yPosition"] == 0:
                wt0_yPosition = Ly / (2 * diameter)
            else:
                wt0_yPosition = 0
        y_act = actuatorInfo["yPosition"] + wt0_yPosition
        
        readme.write(f'   - diameters downstream: {x_act}\n')
        readme.write(f'   - diameters sideways from center: {y_act-Ly/(2*diameter)}\n')
        
        #---CALCULATE DISK CENTER
        diskPoint_x=diameter*x_act
        diskPoint_y=diameter*y_act
        diskPoint_z=height
        includeSettingsFile.write(
            f'// center point of actuator N{actNum}\n' +
            f"diskPoint_x_{actNum} {diskPoint_x};\n" +
            f"diskPoint_y_{actNum} {diskPoint_y};\n" +
            f"diskPoint_z_{actNum} {diskPoint_z};\n" )
        readme.write(f'   - center: ({diskPoint_x},{diskPoint_y},{diskPoint_z})\n')
        
        #---CALCULATE CORNERS OF WAKE AND TURBINE BOXES (WAKE-->SNAPPYHEXMESHDICT AND TURBINE-->TOPOSETDICT)
        #FOR snappyHexMeshDict
        x1b_0=(x_act-1)*diameter
        x1b_1=Lx
        x2b_0=y_act*diameter-1.5*diameter
        x2b_1=y_act*diameter+1.5*diameter
        x3b_0=height-1.5*diameter
        x3b_1=height+1.5*diameter
        #---WRITE THE INPUT DATA FOR snappyHexMeshDict
        includeSettingsFile.write(
            '\n// input para snappyHexMeshDict\n' +
            f"x1b_0_{actNum} {x1b_0};\n" +
            f"x1b_1_{actNum} {x1b_1};\n" +
            f"x2b_0_{actNum} {x2b_0};\n" +
            f"x2b_1_{actNum} {x2b_1};\n" +
            f"x3b_0_{actNum} {x3b_0};\n" +
            f"x3b_1_{actNum} {x3b_1};\n" )     
        
        #FOR topoSetDict
        c_Dref=c_D*2
        l_celda_ref=diameter/c_Dref
        x1d_0=x_act*diameter
        x1d_1=x1d_0 + 3*l_celda_ref
        x2d_0=y_act*diameter-0.5*diameter
        x2d_1=y_act*diameter+0.5*diameter
        x3d_0=height-0.5*diameter
        x3d_1=height+0.5*diameter
        #---WRITE THE INPUT DATA FOR topoSetDict
        includeSettingsFile.write(
            '\n// input para topoSetDict\n' +
            f"x1d_1_{actNum} {x1d_0};\n" +
            f"x1d_2_{actNum} {x1d_1};\n" +
            f"x2d_1_{actNum} {x2d_0};\n" +
            f"x2d_2_{actNum} {x2d_1};\n" +
            f"x3d_1_{actNum} {x3d_0};\n" +
            f"x3d_2_{actNum} {x3d_1};\n" )     
        
        #FVOPTIONS            
        #---OPEN FILE OF ACTUATOR TO INPUT TO fvOptions
        actFile1Dir=f'baseFiles/{actuatorType}.part1'
        actuatorFile1=open(actFile1Dir, 'r')
        actFile2Dir=f'baseFiles/{actuatorType}.part2'
        actuatorFile2=open(actFile2Dir, 'r')
        
        #---WRITE INTO fvOptions FILE
        fOptfile.write('\ndisk' + str(actNum) + '\n')
        fOptfile.write(actuatorFile1.read())
        fOptfile.write(
            f'	      diskPoint         ($diskPoint_x_{actNum} $diskPoint_y_{actNum} $diskPoint_z_{actNum}); //disk center point\n' +
            f'        cellSet           actuationDisk{actNum};\n'
        )
        if actuatorModel != 99: # check that it is the V21 and not the V2
            fOptfile.write(
                '\n' +
                f'        ADmodel                       {actuatorInfo["ADmodel"]};\n' +
                f'        nodesCellsRatio               {actuatorInfo["nodesCellsRatio"]};\n' +
                f'        gradInterpolation             {actuatorInfo["gradInterpolation"]};\n' +
                f'        lambda                        {actuatorInfo["lambda"]};\n' +
                f'        centerRatio                   {actuatorInfo["centerRatio"]};\n' +
                f'        rootDistance                  {actuatorInfo["rootDistance"]};\n' +
                f'        UdCellsMethod                 {actuatorInfo["UdCellsMethod"]};\n' +
                f'        UdCenterToggle                {actuatorInfo["UdCenterToggle"]};\n' +
                f'        UdCorrection                  {actuatorInfo["UdCorrection"]};\n' +
                f'        rootFactor                    {actuatorInfo["rootFactor"]};\n' + 
                f'        tipFactor                     {actuatorInfo["tipFactor"]};\n' +
                f'        forceDistributionMethod       {actuatorInfo["forceDistributionMethod"]};\n' +
                f'        averageChordLength            {actuatorInfo["averageChordLength"]};\n'
            ) 

        fOptfile.write(actuatorFile2.read())
        
        
        #---WRITE INTO topoSetDict FILE
        #tSDfile.write()
        tSDfile.write('    {' + f'\n        name    actuationDisk{actNum}CellSet;\n')
        tSDpart1=open('baseFiles/topoSet.part1','r')
        tSDfile.write(tSDpart1.read())
        tSDpart1.close()
        tSDfile.write(
            f'        box     ($x1d_1_{actNum} $x2d_1_{actNum} $x3d_1_{actNum}) ($x1d_2_{actNum} $x2d_2_{actNum} $x3d_2_{actNum});\n' +
            '   }\n' +
            '    {\n'+
            f'        name    actuationDisk{actNum};\n'
        )
        tSDpart2=open('baseFiles/topoSet.part2','r')
        tSDfile.write(tSDpart2.read())
        tSDpart2.close()
        tSDfile.write(f'        set     actuationDisk{actNum}CellSet;\n    ' + '}\n')
        
        
        #---WRITE INTO snappyHexMeshDict FILE
        sHMDfile.write(
            f'    windTurbine{actNum}\n' + 
            '    {\n' + 
            '        type box;\n' +
            f'        min  ( $x1b_0_{actNum} $x2b_0_{actNum} $x3b_0_{actNum} );\n' + 
            f'        max  ( $x1b_1_{actNum} $x2b_1_{actNum} $x3b_1_{actNum} );\n' + 
            '    }\n')
        sHMDaux=sHMDaux + f'        windTurbine{actNum} \n'+'        {\n            mode inside;\n            levels ((1 1));\n        }\n'
        
        includeSettingsFile.write('//---------------------------------\n')
    
    #CLOSE fvOptions FILE CREATION
    fOptfile.close()
    
    #CLOSE topoSetDict FILE CREATION
    tSDfile.write('\n);')
    tSDfile.close()
    
    #CLOSE snappyHexMeshDict FILE CREATION
    sHMDmiddle=open('baseFiles/snappyHexMeshDict.middle','r')
    sHMDfile.write(sHMDmiddle.read())
    sHMDmiddle.close()
    sHMDfile.write(sHMDaux)
    sHMDend=open('baseFiles/snappyHexMeshDict.end','r')
    sHMDfile.write(sHMDend.read())
    sHMDend.close()
    sHMDfile.close()
        
    ## variables para el fvOptions
    includeSettingsFile.write("\n// variables para el fvOptions\n")
    
    ## Uinf
    if problemDict["inletVelocities"][1] > problemDict["inletVelocities"][0]:
        varU = True
    else:
        varU = False
        
    #VARIAS VELOCIDADES
    if varU:
        Umin=problemDict["inletVelocities"][0]
        print(f"\nVelocidad minima = {Umin} m/s")

        Umax=problemDict["inletVelocities"][1]
        print(f"Velocidad maxima = {Umax} m/s")
        
        Ustep=problemDict["inletVelocities"][2]
        print(f"Step de velocidad = {Ustep}m/s")
        
        problemSettingsFile.write(
            "\n# min, max and step velocities\n" + 
            f"Umin={Umin}\n" +
            f"Umax={Umax}\n" + 
            f"Ustep={Ustep}\n")
        readme.write(
            "\nUstep=on\n" + 
            f"Umin={Umin}\n" +
            f"Umax={Umax}\n" + 
            f"Ustep={Ustep}\n")
    #UNICA VELOCIDAD
    else:
        Uinf=problemDict["inletVelocities"][0]
        print(f"\nVelocidad de entrada = {Uinf} m/s")
        problemSettingsFile.write(
            "\n# velocity with step off\n" +
            f"Uinf={Uinf}\n" +
            f"Ustep='off'\n")
        readme.write("\nUstep=off\nUinf={Uinf}\n")

    ## Uref
    if problemDict["referenceVelocities"] != []:
        Urefs = problemDict["referenceVelocities"]
        problemSettingsFile.write(f'UrefDiff=\'yes\'\n')
        if Urefs[0] < Urefs[1]: 
            varUref = True
        else:
            varUref = False
            
        #VARIAS VELOCIDADES
        if varUref:
            UrefMin=Urefs[0]
            print(f"\nVelocidad de referencia minima = {UrefMin} m/s")

            UrefMax=Urefs[1]
            print(f"Velocidad de referencia maxima = {UrefMax} m/s")
            
            UrefStep=Urefs[2]
            print(f"Step de velocidad de referencia = {UrefStep}m/s")
            
            problemSettingsFile.write(
                "\n# min, max and step reference velocities\n" +
                f"UrefMin={UrefMin}\n" +
                f"UrefMax={UrefMax}\n" +
                f"UrefStep={UrefStep}\n")
            readme.write(
                "\nUrefStep=on\n" + 
                f"UrefMin={UrefMin}\n" + 
                f"UrefMax={UrefMax}\n" + 
                f"UrefStep={UrefStep}\n")
        #UNICA VELOCIDAD
        else:
            Uref=Urefs[0]
            print(f"\nVelocidad de referencia = {Uref} m/s")
            problemSettingsFile.write(
                "\n# reference velocity with step off\n" +
                f"Uref={Uref}\n" +
                f"UrefStep='off'\n")
            readme.write(
                "\nUrefStep=off\n" + 
                f"Uref={Uref}\n")
    else:
        problemSettingsFile.write(f'UrefDiff=\'no\'\n')

    yaw=0
    print(f"se considera yaw={yaw}")
    includeSettingsFile.write(f"yaw {yaw};\n")

    # pitch=0
    # print(f"se considera pitch={pitch}")
    # includeSettingsFile.write(f"pitch {pitch};\n")

    # omega=1
    # print(f"se considera omega={omega}")
    # includeSettingsFile.write(f"omega {omega};\n")

    cellSize=l_celda/2 #considero que el refinado de la zona del actuador y la estela es el doble de la zona no refinada
    print(f"se considera cellSize={cellSize}")
    includeSettingsFile.write(f"cellSize {cellSize};\n")

    Area=math.pi*0.25*diameter**2
    includeSettingsFile.write(f"Area {Area};\n")
    

    ## variables para el initalConditions
    if problemDict["turbulentIntensities"][0] < problemDict["turbulentIntensities"][1]:
        varTI = True
    else:
        varTI = False
        
    #VARIAS VELOCIDADES
    if varTI:
        TImin=problemDict["turbulentIntensities"][0]
        print(f"Intensidad turbulenta minima = {TImin}")
        
        TImax=problemDict["turbulentIntensities"][1]
        print(f"Intensidad turbulenta maxima = {TImax}")
        
        TIstep=problemDict["turbulentIntensities"][2]
        print(f"Step de intensidad turbulenta = {TIstep}")
        
        problemSettingsFile.write(
            "\n# min, max and step turbulent intensities\n" +
            f"TImin={TImin}\n" +
            f"TImax={TImax}\n" +
            f"TIstep={TIstep}\n")

        readme.write(
            "\nTIstep=on\n" + 
            f"TImin={TImin}\n" + 
            f"TImax={TImax}\n" + 
            f"TIstep={TIstep}\n")
    #UNICA VELOCIDAD
    else:
        TI=problemDict["turbulentIntensities"][0]
        print(f"Intensidad turbulenta = {TI}")
        
        problemSettingsFile.write(
            "\n# TI with step off\n" +
            f"TI={TI}\n" +
            f"TIstep='off'\n")
        readme.write(
            "\nTIstep=off\n" + 
            f"TI={TI}\n")

    #---SAMPLE SETUP
    samplefilename=f'problem{problemNumber}Files/sample{problemNumber}'
    copy('baseFiles/sample.base', samplefilename) 
    samplefile = open(samplefilename, 'a+')
    x_samples = problemDict["lineSamples"][0]
    y_samples = problemDict["lineSamples"][1]
    z_samples = problemDict["lineSamples"][2]
    
    for sampleNr in range(len(x_samples["positions"])):
        xD = x_samples["positions"][sampleNr]
        samplefile.write(
            f'    xD{xD}\n' +
            '    {\n' +
            '        type            uniform;\n' +
            '        axis            x;\n' + 
            '        start           (0 #eval{' + f' $y_2 * 0.5 + $D * {xD}' + ' } $H);\n' + 
            '        end             ($x_2 #eval{' + f' $y_2 * 0.5 + $D * {xD}' + ' } $H);\n'+ 
            '        nPoints         1000;\n    }\n')

    for sampleNr in range(len(y_samples["positions"])):
        yD = y_samples["positions"][sampleNr]
        samplefile.write(
            f'    yD{yD}\n' +
            '    {\n' +
            '        type            uniform;\n' +
            '        axis            y;\n' + 
            '        start           (#eval{' + f' $D * {yD}' + ' } 0 $H);\n' + 
            '        end             (#eval{' + f' $D * {yD}' + ' } $y_2 $H);\n'+ 
            '        nPoints         1000;\n    }\n')

    for sampleNr in range(len(z_samples["positions"])):
        zD = z_samples["positions"][sampleNr]
        samplefile.write(
            f'    zD{zD}\n' +
            '    {\n' +
            '        type            uniform;\n' +
            '        axis            z;\n' + 
            '        start           (#eval{' + f' $D * {zD}' + ' }' + ' #eval{ $y_2 / 2 } 0);\n' + 
            '        end             (#eval{' + f' $D * {zD}' + ' }' + ' #eval{ $y_2 / 2 } $z_7);\n'+ 
            '        nPoints         1000;\n    }\n')
    samplefile.write(');')
    samplefile.close()    

    #---SURFACES SETUP
    surfacesfilename=f'problem{problemNumber}Files/surfaces{problemNumber}'
    copy('baseFiles/surfaces.base', surfacesfilename) 
    surfacesfile = open(surfacesfilename, 'a+')
    for surfaceSamples in problemDict["surfaceSamples"]:
        if surfaceSamples["normalAxis"] == "x":
            x_surfaceSample = surfaceSamples
            counter = 1
            for position in x_surfaceSample["positions"]:
                surfacesfile.write(
                    f'        xNormal{counter}\n' +
                    '        {\n' +
                    '            type            cuttingPlane;\n' +
                    '            planeType       pointAndNormal;\n' +
                    '            pointAndNormalDict\n' +
                    '            {\n' +
                    f'                point   ({position} 0 0);\n' +
                    '                normal  (1 0 0);\n' +
                    '            }\n' +
                    '            interpolate    true;\n' +
                    '        }\n'
                )
                counter += 1
        elif surfaceSamples["normalAxis"] == "y":
            y_surfaceSample = surfaceSamples
            counter = 1
            for position in y_surfaceSample["positions"]:
                surfacesfile.write(
                    f'        yNormal{counter}\n' +
                    '        {\n' +
                    '            type            cuttingPlane;\n' +
                    '            planeType       pointAndNormal;\n' +
                    '            pointAndNormalDict\n' +
                    '            {\n' +
                    f'                point   (0 {position} 0);\n' +
                    '                normal  (0 1 0);\n' +
                    '            }\n' +
                    '            interpolate    true;\n' +
                    '        }\n'
                )
                counter += 1
        elif surfaceSamples["normalAxis"] == "z":
            z_surfaceSample = surfaceSamples
            counter = 1
            for position in z_surfaceSample["positions"]:
                surfacesfile.write(
                    f'        zNormal{counter}\n' +
                    '        {\n' +
                    '            type            cuttingPlane;\n' +
                    '            planeType       pointAndNormal;\n' +
                    '            pointAndNormalDict\n' +
                    '            {\n' +
                    f'                point   (0 0 {position});\n' +
                    '                normal  (0 0 1);\n' +
                    '            }\n' +
                    '            interpolate    true;\n' +
                    '        }\n' 
                )
                counter += 1

    surfacesfile.write('    }\n}')
    surfacesfile.close()    



    #---FINAL CLOSE
    includeSettingsFile.close()
    problemSettingsFile.close()
    readme.close()

