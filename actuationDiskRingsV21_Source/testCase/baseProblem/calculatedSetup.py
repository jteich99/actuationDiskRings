#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
----------------------------------------
created by Juan Ignacio Teich
Facultad de Ingeniería, UBA
jteich@fi.uba.ar
----------------------------------------
"""

def is_float(string):
    try:
        # float() is a built-in function
        float(string)
        return True
    except ValueError:
        return False

def look_Ct_Cp(Cp_Ct_table_csv,Uinlet):
    #Search for Ct and Cp according to Uinlet
    inputcsv=csv.reader(open(Cp_Ct_table_csv,'r'))
    inputlist=list(inputcsv)
    upperCheck = True
    UinletLower = 0
    for row in inputlist:
        if is_float(row[0])==True:
            if float(row[0])<=Uinlet:
                UinletLower = float(row[0])
                CtLower = float(row[2])
                CpLower = float(row[3])
            elif (
                    float(row[0]) > Uinlet and
                    upperCheck
                ):
                UinletUpper = float(row[0])
                CtUpper = float(row[2])
                CpUpper = float(row[3])
                upperCheck = False 
    if UinletLower == Uinlet:
        return CtLower,CpLower
    else:
        Ct = CtLower + (UinletUpper - UinletLower) * (CtUpper - CtLower)
        Cp = CpLower + (UinletUpper - UinletLower) * (CpUpper - CpLower)
        return Ct,Cp

def look_omega_pitch(referenceTable, Uinlet):
    inputcsv=csv.reader(open(referenceTable,'r'))
    inputlist=list(inputcsv)
    upperCheck = True
    UinletLower = 0
    for row in inputlist:
        if is_float(row[0])==True:
            if float(row[0])<=Uinlet:
                UinletLower = float(row[0])
                omegaLower = float(row[1])
                pitchLower = float(row[2])
            elif (
                    float(row[0]) > Uinlet and
                    upperCheck
                ):
                UinletUpper = float(row[0])
                omegaUpper = float(row[1])
                pitchUpper = float(row[2])
                upperCheck = False 
    if UinletLower == Uinlet:
        return omegaLower,pitchLower
    else:
        omega = omegaLower + (UinletUpper - UinletLower) * (omegaUpper - omegaLower)
        pitch = pitchLower + (UinletUpper - UinletLower) * (pitchUpper - pitchLower)
        return omega,pitch
    

import csv
import os
import math

Uinf=float(os.environ["Uinf"])
Uref=float(os.environ["Uref"])
TI=float(os.environ["TI"])
D=float(os.environ["D"])
H=float(os.environ["H"])

includeSettingsFile=open("includeSettings","a")

includeSettingsFile.write("\n// variables propias de este problema\n")
includeSettingsFile.write("Uinf "+ str(Uinf) + ";\n")

TI *= 0.8

includeSettingsFile.write("TI "+ str(TI) + ";\n")

# TURBULENT KINETIC ENERGY
TKE=1.5*(Uinf*TI)**2
includeSettingsFile.write("TKE " + str(TKE) + ";\n")

C_mu=0.09
l=0.07*D
Kappa=0.408
TEpsilon=(C_mu*TKE**(1.5))/l
includeSettingsFile.write("TEpsilon " + str(TEpsilon) + ";\n")

TOmega=math.sqrt(TKE)/l
includeSettingsFile.write("TOmega " + str(TOmega) + ";\n")

Z0_TI = H/(math.exp((Kappa*C_mu**(-1/4)*math.sqrt(2/3))/(0.8*TI)))
includeSettingsFile.write("Z0_TI " + str(Z0_TI) + ";\n")

Ct,Cp = look_Ct_Cp('Tabla_UdAvg_Cp_Ct.csv',Uref)
includeSettingsFile.write("Ct " + str(Ct) + ";\n")
includeSettingsFile.write("Cp " + str(Cp) + ";\n")
omega,pitch = look_omega_pitch('NRELtable.csv',Uref)
includeSettingsFile.write("omega " + str(omega) + ";\n")
includeSettingsFile.write("pitch " + str(pitch) + ";\n")
