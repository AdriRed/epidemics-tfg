#!/bin/bash

# Directorios
INCLUDE_DIR="include"
BUILD_DIR="builds"

# Función para compilar un archivo si es más reciente que su .o
compile_file() {
    local file=$1
    local base_name=$(basename "$file" .f90)
    local base_name=$(basename "$base_name" .f)
    local obj_file="$BUILD_DIR/${base_name}.o"
    
    # Si el archivo objeto no existe, necesitamos compilar
    if [ ! -f "$obj_file" ]; then
        echo "Compilando: $file (archivo objeto no existe)"
        gfortran -c -g -Og -fimplicit-none -fcheck=all -fbacktrace "$file" -o "$obj_file"
        if [ $? -ne 0 ]; then
            echo "Error compilando $file"
            return 1
        fi
        # Copiar archivos .mod a include para el linter
        cp "$BUILD_DIR"/*.mod "$INCLUDE_DIR/" 2>/dev/null || true
        return 0
    fi
    
    # Compilar si el archivo fuente es más reciente que el .o
    if [ "$file" -nt "$obj_file" ]; then
        echo "Compilando: $file (fuente más reciente)"
        gfortran -c -g -Og -fimplicit-none -fcheck=all -fbacktrace "$file" -o "$obj_file"
        if [ $? -ne 0 ]; then
            echo "Error compilando $file"
            return 1
        fi
        # Copiar archivos .mod a include para el linter
        cp "$BUILD_DIR"/*.mod "$INCLUDE_DIR/" 2>/dev/null || true
        return 0
    else
        echo "Saltando: $file (ya compilado)"
        return 2  # Código diferente para "ya compilado"
    fi
}

# Función para compilar con resolución de dependencias
compile_with_dependencies() {
    local files=("$INCLUDE_DIR"/*.f90 "$INCLUDE_DIR"/*.f)
    local compiled_count=0
    local max_iterations=${#files[@]}
    local iteration=0
    
    # Filtrar solo archivos que existen
    local existing_files=()
    for file in "${files[@]}"; do
        if [ -f "$file" ]; then
            existing_files+=("$file")
        fi
    done
    
    files=("${existing_files[@]}")
    
    if [ ${#files[@]} -eq 0 ]; then
        echo "No hay archivos Fortran en $INCLUDE_DIR"
        return 0
    fi
    
    echo "Archivos a compilar: ${files[@]}"
    
    while [ ${#files[@]} -gt 0 ] && [ $iteration -lt $max_iterations ]; do
        local remaining_files=()
        local compiled_this_round=0
        
        echo "--- Iteración $((iteration + 1)) ---"
        
        for file in "${files[@]}"; do
            echo "Procesando: $file"
            if compile_file "$file"; then
                # Éxito en compilación (código 0)
                ((compiled_this_round++))
                ((compiled_count++))
                echo "✓ Compilado: $file"
            elif [ $? -eq 2 ]; then
                # Ya compilado (código 2) - consideramos éxito pero no cuenta para compiled_this_round
                ((compiled_count++))
                echo "✓ Ya compilado: $file"
            else
                # Error o dependencia faltante (código 1)
                remaining_files+=("$file")
                echo "→ Pendiente: $file"
            fi
        done
        
        files=("${remaining_files[@]}")
        ((iteration++))
        
        echo "Iteración $iteration: Compilados $compiled_this_round, Pendientes: ${#files[@]}"
        
        # Si no compilamos nada pero aún hay archivos pendientes, hay un problema
        if [ $compiled_this_round -eq 0 ] && [ ${#files[@]} -gt 0 ]; then
            echo "Error: No se progresó en la iteración $iteration. Archivos pendientes:"
            for file in "${files[@]}"; do
                echo "  - $file"
            done
            return 1
        fi
    done
    
    if [ ${#files[@]} -gt 0 ]; then
        echo "Error: No se pudieron compilar todos los archivos después de $iteration iteraciones:"
        for file in "${files[@]}"; do
            echo "  - $file"
        done
        return 1
    fi
    
    echo "✅ Compilación completada: $compiled_count archivos en $iteration iteraciones"
    return 0
}

# Limpiar archivos .mod viejos en builds (opcional, comenta si no quieres)
# rm -f "$BUILD_DIR"/*.mod

# Ejecutar compilación con resolución de dependencias
if compile_with_dependencies; then
    # Copiar todos los .mod a include para el linter
    cp "$BUILD_DIR"/*.mod "$INCLUDE_DIR/" 2>/dev/null || true
    echo "✅ Compilación de módulos completada exitosamente"
else
    echo "❌ Error en la compilación de módulos"
    exit 1
fi