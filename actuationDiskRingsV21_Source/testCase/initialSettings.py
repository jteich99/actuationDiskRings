"""
----------------------------------------
created by Juan Ignacio Teich
Facultad de Ingeniería, UBA
jteich@fi.uba.ar
----------------------------------------
"""
tupac = False
problems = [
    {
        # problem 1
        "nodes": 1,
        # "cores": 64,
        "cores": 6,
        "diameter": 126,
        "nacelleHeight": 90, # if ABL BC, then real height for nacelle = 90m
        "problemMeshNumber": 0, # if =0 then it creates the mesh with the data, if = another problemNumber it copies the mesh from than problem
        "domain": [25, 9, 9],
        "cellsPerDiameter": 10,
        "boundaryCondition": "SLIP", # SLIP or ABL
        "actuators": [
            {
                "ADmodel": 0, #  0=uniform, 1=numeric, 2=analytic original, 3=analytic non uniform flow, 4=generalized analytic, 5=eliptic, 99=airfoil
                "xPosition": 5, # x position of actuator in [D]
                "yPosition": 0, # y position of actuator in [D] relative to 1st AD
                "nodesCellsRatio": 2,
                "rootFactor": 1,
                "tipFactor": 1,
                "forceDistributionMethod": 1,
                # "rThicknesCellsizeRatio": 0.5,
                "gradInterpolation": 1,
                "rootDistance": 0.14,
                "lambda": 7,
                "centerRatio": 0.3,
                "UdCellsMethod": 2,
                "UdCenterToggle": 0,
                "UdCorrection": 0,
                "averageChordLength": 3.46
            },
            {
                "ADmodel": 4, #  0=uniform, 1=numeric, 2=analytic original, 3=analytic non uniform flow, 4=generalized analytic, 5=eliptic, 99=airfoil
                "xPosition": 10, # x position of actuator in [D]
                "yPosition": 1, # y position of actuator in [D] relative to 1st AD
                "nodesCellsRatio": 2,
                "rootFactor": 1,
                "tipFactor": 1,
                "forceDistributionMethod": 1,
                # "rThicknesCellsizeRatio": 0.5,
                "gradInterpolation": 1,
                "rootDistance": 0.14,
                "lambda": 7,
                "centerRatio": 0.3,
                "UdCellsMethod": 2,
                "UdCenterToggle": 0,
                "UdCorrection": 0,
                "averageChordLength": 3.46
            },
        ],
        "inletVelocities": [10,10,1], #minU, maxU, stepU
        "referenceVelocities": [],
        "turbulentIntensities": [0.05, 0.05, 0.01],
        "surfaceSamples": [
            {
                "normalAxis": "z",
                "positions": [9 * 126 / 2]
            },
            {
                "normalAxis": "x",
                "positions": [126 * 5, 126 * 8, 126 * 10, 126 * 13]
            },
        ],
        "lineSamples": [
            {
                "axis": "x",
                "positions": [0, 1]
            },
            {
                "axis": "y",
                "positions": [1, 3, 5, 7, 9, 11, 13, 15, 20, 24]
            },
            {
                "axis": "z",
                "positions": []
            }
        ]
    }
]
