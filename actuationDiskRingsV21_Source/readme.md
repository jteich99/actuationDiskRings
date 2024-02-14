# README to actuationDiskRingsV21_Source
> The main idea behind this AD is to compile the features to be compared in Juan Ignacio Teich's thesis into 1 actuator and to have toggles to decide what option to use.

## Header file

## Source file

## SourceTemplate file
Order of code execution:
- [Save cells in disc](#save-cells-in-disc)
- [Disc direction from yaw angle](#disc-direction-form-yaw-angle)
- [Ud calculation](#ud-calculation)
- [Forces calculation](#forces-calculation)
    - [Numeric AD](#numeric-ad-method)
    - [Force distribution](#force-distribution)
    - [Tip and root factor](#tip-and-root-factor)
        - [Tip factor](#tip-factor)
        - [Root factor](#root-factor)

---

### Save cells in disc
>Pending
>- [ ] move to begining of file as first thing to do
>- [ ] move to Source file to avoid calculatin in each time step
>   - [ ] leave commented since in AL it wouldn't be possible to store it in the first time step.

Loop through all cells and if within sphere of 1.15 R from disc center, and within $3\cdot E$ from disc plane, append cell to **cellsDisc** list.

---

### Disc direction from yaw angle
>Pending
>- [ ] ahora esta implementado antes de **Save cells in disc**, pero deberia implementarse despues ya sobre la lista de las celdas que pertencen al disco en lugar de a todas las celdas del dominio.

If *yaw angle* is 360, then:
- Calculate average velocity in center of disc using a *centerRatio*, **U_dCenterCells_orient**.
- Calculate self orientation depending on velocity on the disc by normalizing average velocity vector, **uniDiskDir**.
- Calculate angle in between self orientation direction (*uniDiskDir*) and disc direction (*diskDir_*), **Alpha**.

Else:
- Calculate direction of yawed disc as **diskYawed**, and make **uniDiskDir** equal to it.

---

### Ud calculation
>Pending
>- [ ] verify that the weight in sphere is what is expected. Plot and verify.

Depending on **UdCellsMethod** chosen in *fvOptions* the average velocity in the disc is calculated. The options change the weighting method used (see below).

Then, according to selection of **UdCenterToggle**, the **U_dCells** is calculated considering weight in full disc or in disc center. Additionaly **U_dCenterCells** is always calculated with the center weight.

#### Method 0 - no weighting
All weights are 1, then all cells are weighted equal, and the direct average is obtained.

#### Method 1 - weight according to distance to AD plane
Weight to AD plane is computed as:
$$
\text{weightADplane} = \frac{1}{E \cdot \sqrt{\pi}} \cdot e ^ {- \left( \frac{d}{E} \right)^ 2}
$$
Being $E$ the smearing factor and $d$ the distance from de cell center to the AD plane.

#### Method 2 - weight according to distance to AD plane and distance to AD center
Weight to AD plane is computed as in [Method 1](#method-1-weight-according-to-distance-to-ad-plane). Then the weight for the distance to de AD center is computed in 2 ways, one for the full disc and another for the center of the disc. Afterwards, when calculating the average velocity, through the *fvOptions* toggle **UdCenterToggle** it is chosen if the center of the disc or the full disc is used for the calculation.

The weights are computed as:
$$
\text{weightSphereCenter} = \frac{1}{(\text{centerRatio} \cdot R) \cdot \sqrt{\pi}} \cdot e ^ {- \left( \frac{dSphere}{(\text{centerRatio} \cdot R)} \right)^ 2}
$$
$$
\text{weightSphereAD} = \frac{1}{R \cdot \sqrt{\pi}} \cdot e ^ {- \left( \frac{dSphere}{R} \right)^ 2}
$$

---

### Forces calculation
>Pending
>- [ ] Add a variable from *fvOptions* that allows to chose what basic AD method will be used. For now only allow for 1: Numeric AD.

#### Numeric AD method
>Method from Navarro-Diaz 2019.

>Pending
>- [ ] test that loops going over central node works as expected 
>- [ ] test if it can be generalized the operations done for each position of the table into 1 block of code that is in a loop and runs for position 1 and 2.

1. The positions for interpolation in table 1 (*UdAvgList*) are calculated using **U_dCells**.
2. From **U_dCells** and the positions in table 1, interpolate:
    - Uref
    - omega
    - pitch
    - Ct
    - Cp
3. Calculate thrust, power and torque.
4. Initialization of multiple parameteres used afterwards.
5. Loop over rings, and then in each ring over nodes.
    1. Calculate node position.
    2. Calculate vector from disc center to node and tangential vector.
    3. Perform coordinate system change from xyz (cartesian) system to ntr (cylindric) system. 
    4. Calculate velocity in node. If **gradInterpolation** is on, then use *grad(U)* in cell and distance to cell center for better result.
    5. Loop over cells in disc for weighting according to distance to node.
        1. Calculate distance from cell to node.
        2. According to **forceDistributionMethod** selected in *fvOptions*, weigh according to distance.
    6. Based on positions in table1 calculated before, get from inteprolation of table 2 (*UdiList*) the values of **Uinf**, **fn** and **ft** for each position.
    7. Inteprolate the values of **Uinf**, **fn** and **ft** from the values of position 1 and 2.
    8. From **fn** and **ft** (forces over area) calculate the total forces corresponding to the node **F_n_Bi** and **F_t_Bi**.
    9. Loop over cells in disc for force application 
6. Write output files *outRings* and *outTurbines*.

#### Force distribution
>Pending
>- [ ] In the way that the velocity is calculated for the node using the velocity and gradient in the cell it belongs, a similar approach can be used when applying the forces directly on the cell the node belongs to. The total_nodes_counter is used for this operation.
##### Method 0 - no distribution
##### Method 1 - force distribution originally used in Navarro Diaz code
##### Method 2 - force distribution proposed in Mikkelsen 2003


<!-- #### Uinf estimation -->

---

#### Tip and root factor
>Pending
>- [ ] Apply tip and root factor coded in V11_calib_Source. (when using the Numeric AD they aren't directly used since they are used in the calibration phase) 
##### Tip factor
##### Root factor

