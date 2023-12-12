Acá define `U_dPointCells` como un vector que **no se que es**, pero aparte, en el `if`, se define el vector `dU` como un producto interno entre 2 vectores (`dx` definido antes, y `gradU` de la celda).  ([link de user guide que dice que & es producto interno](https://www.openfoam.com/documentation/guides/v2112/doc/openfoam-guide-expression-syntax.html)).

**Esto no daría un escalar como resultado en lugar de un vector?**
```cpp
// line 457 of actuationDiskRingsV11_calib_SourceTemplates.C
            vector U_dPointCells = vector(1000,1000,1000);
            if (nodeCellID_[total_nodes_counter] != -1) //if the closer cell is in this procesor
            {
	            U_dPointCells =  U[nodeCellID_[total_nodes_counter]];
        		if (gradInterpolation_ == 1) // ?? condition
        		{     		
            		vector dx = Bi - mesh().cellCentres()[nodeCellID_[total_nodes_counter]];
		            vector dU = dx & gradU[nodeCellID_[total_nodes_counter]]; // & = internal product of 2 vectors in OF --> scalar...
		            U_dPointCells += dU;
                }
            }
            reduce(U_dPointCells, minOp<vector>()); //take only normal values of U

            if (mag(U_dPointCells) > 1000) // We add a flag in case it does not find a cell near
            {
	            U_dPointCells =  vector(10,0,0);
                Info<<"OpenFOAM cell Not found"<<endl;
                Info<<"ring: "<<ring<<endl;
                Info<<"node: "<<total_nodes_counter<<endl;
                Info<<"radius: "<< radius <<endl;
            
            }
```
