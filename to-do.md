# priority

# backlog
## actuadores
- [ ] incorporar AD analitico generalizado de Sorensen 2023
- [ ] incorporar AD de van Der Laan 2015
- [ ] incorporar AD con fuerzas elípticas de Sorensen 1992 

## distribución de fuerzas
- [ ] incorporar función $delta$ de Li 2022
- [ ] incorporar distribución de fuerzas de Martínez Tossas 2017 

## velocidades
- [ ] incorporar correción de velocidad promedio en disco de Li 2022


## funcional
- [ ] hacer que guarde las velocidades en los nodos después de la primera vez que lo corre
    - capaz se puede hacer que lo mida de una en cada iteración y después la cantidad de veces que lo necesite simplemente lo busca en una lista
- [ ] migrar el actuador a fuera de las carpetas de instalación y a la carpeta de usuario


# done
## 04-09
- [x] migrar a funciones getNodePosition y getNodeVelocity a todos los bloques donde se hacen esas acciones
- [x] probar poner en include el archivo de funciones y agregar al archivo files del fvOptions para la compilación 
