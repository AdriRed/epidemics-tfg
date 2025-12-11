#!/bin/bash
# Limpia la solución Fortran del proyecto epidemics-tfg

# Elimina carpeta builds si existe
if [ -d "builds" ]; then
    echo "Eliminando carpeta builds/"
    rm -rf builds
fi

# Borra todos los archivos .mod del directorio raíz
echo "Eliminando archivos .mod de la raíz"
rm -f ./*.mod

# Borra todos los archivos .mod dentro de include/
if [ -d "include" ]; then
    echo "Eliminando archivos .mod de include/"
    rm -f ./include/*.mod
fi

echo "Limpieza completada."