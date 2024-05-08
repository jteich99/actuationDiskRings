# priority

# backlog
## actuadores
- [ ] incorporar AD de van Der Laan 2015
- [ ] incoporar computo de fuerzas en la góndola por separado

## distribución de fuerzas
- [ ] hacer que sea una función como con los tip y root factor para tener menos clutter en el código

## velocidades
- [ ] incorporar medicion de velocidad aguas abajo y aguas arriba

## funcional
- [ ] migrar el actuador a fuera de las carpetas de instalación y a la carpeta de usuario
- [ ] implementar Ct de funcion analítica como alternativa a la power curve



# done
## 06-05
- [x] incorporar correción de velocidad promedio en disco de Li 2022
## 02-05
- [x] incorporar distribución de fuerzas de Martínez Tossas 2017 
    - [x] agregar en fvOptions que se ponga la distancia promedio de cuerda
    - [x] poner kernel anisotrópico
## 29-04
- [x] revisar AD analitico original por potencial error
    - [x] corregir que use las nuevas funciones para root y tip factor y no las viejas gFunction y FFunction
    - lo que traia el error era usar U0 del nodo en las cuentas y no la U0 global o Uref
## 04-25
- [x] corregir generación de tablas del V11_Source a partir d ela tabla de calibración 
- [x] agregar factor de corrección para tener empuje total al AD eliptico
## 04-19
- [x] incorporar AD con fuerzas elípticas de Sorensen 1992 
- [x] limpiar del V21_Source.C las partes que son del numérico poniendo un condicional
## 04-18
- [x] incorporar AD uniforme
## 04-17
- [x] incorporar AD analitico generalizado de Sorensen 2024
- [x] renombrar las opciones UrefCalculationMethod y forceCalculationMethod a ADmodel para que sea más sencilla la configuración del fvOptions
## 04-16
- [x] hacer funciones de tip y root factor para reemplazar g y F en el actuador analítico, que tomen el tip y root factor especificado en el fvOptions
## 04-13
- [x] generalizar integral de Cp para flujo no uniforme
## 04-10
- [x] incorporar función $delta$ de Li 2022
- [x] hacer que guarde las velocidades en los nodos después de la primera vez que lo corre
## 04-09
- [x] migrar a funciones getNodePosition y getNodeVelocity a todos los bloques donde se hacen esas acciones
- [x] probar poner en include el archivo de funciones y agregar al archivo files del fvOptions para la compilación 
