# implementación actual
```cpp
//Shen tip correction factor:
scalar fcorr=0.0;
scalar tipfactor = 1;
scalar tipfactor_f = (nblades_/2)*(maxR - radius)/(radius*sin(phi));
if (tipFactor_ == 1) //If tip factor is on
{
    scalar c1=0.125;
    scalar c2=27; //27
    scalar c3=0.1;
    scalar g =1;
    //scalar tipfactor_f = (nblades_/2)*(maxR - radius)/(radius*sin(phi)); Movemos antes del if
    //Info << "tipfactor_f " << tipfactor_f << endl;
    g=exp(-c1*(nblades_*tsr-c2))+c3;
    if (tipfactor_f > 0)
            {
            if(((exp(-g*tipfactor_f))>-1) and ((exp(-g*tipfactor_f))<1))
                    {
                    tipfactor = (2/(M_PI))*acos(exp(-g*tipfactor_f));
                    }
            }
}

//Info << "tipfactor: "<<tipfactor<<endl;
//Glauert root correction factor:
scalar rootfactor = 1;
if (rootFactor_ == 1) //If root factor is on
{
    if (radius <=root)
    {
            rootfactor = 0;
    }
    
    if ( (radius > root) and (tipfactor_f>0) and (radius< maxR/2.0)   )
    {
            scalar rootfactor_f = (nblades_/2)*(radius - 0.1*maxR)/(radius*sin(phi));
            scalar g=1;
            
            if ( ((exp(-g*rootfactor_f))>-1) and ((exp(-g*rootfactor_f))<1))
            {
                    rootfactor = (2/(M_PI))*acos(exp(-g*rootfactor_f));
            }
            
    }
}
```
$$
\begin{align*}
\text{tip factor} =
\begin{cases}
\mathcal{F}_{\text{tip}} = \cfrac{2}{\pi} \arccos{ \left( e^{-g\cdot f_{\text{tip}}} \right)} \\
f_{\text{tip}} = \cfrac{n_\text{blades}}{2} \cdot \cfrac{R - R_{i}}{R_i \cdot \sin{\phi}} \\
-1 < g \cdot f_\text{tip} < 1
\end{cases}&& &&
\text{root factor}=
\begin{cases}
\mathcal{F}_{\text{root}} = \cfrac{2}{\pi} \arccos{ \left( e^{- f_{\text{root}}} \right)} \\
f_\text{root} = \cfrac{n_\text{blades}}{2} \cdot \cfrac{R_i - 0.1 \cdot R}{R_i \cdot \sin{\phi}} \\
\begin{cases}
R_\text{root} < R_i < \cfrac{R}{2} \\
f_\text{tip} > 0
\end{cases} \\
R_i =< R_\text{root} \Longrightarrow \mathcal{F}_\text{root} = 0
\end{cases} \\
\end{align*} \\
\begin{gather*}
\phi = \arctan{\cfrac{U_n}{ R_i \cdot \omega - U_t}} \\
g = e^{ -c_1 \cdot (n_\text{blades} \cdot \text{tsr} - c_2) } + c_3 \\
c_1 ; c_2 ; c_3 = \text{constants} \\
\text{tsr} = \cfrac{R \cdot \omega}{U_\text{ref, Yaw}}
\end{gather*}
$$

# implementación de paper
$$
\begin{align*}
\mathcal{F}_\text{tip} = \frac{2}{\pi} \arccos{\Big( e^{ \left[- \cfrac{B\cdot (1-x)}{2\cdot x\cdot \sin{\phi}} \right]} \Big)}
&& &&
\begin{cases}
\mathcal{F}_\text{root} = 1 - e^{\normalsize \left[ -a \left( x / \overline{\delta} \right)^{\normalsize b}\right]} \\
a = 2.335 \quad b = 4 \\
\overline{\delta} = \cfrac{\delta}{R} = \text{vortex core} = \cfrac{R_\text{root}}{R} \\
x = R_i
\end{cases}
\end{align*}
$$