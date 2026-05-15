# Epidemics Engine

Adrià Rojo. 2026.

---

Aquest codi està preparat per executar-se amb Visual Studio Code a Linux.

Té com a requisits tenir el `gfortran` instal·lat.

Un cop estigui tot preparat, només s'ha de compilar i executar amb els diversos script de `build.sh`, `build-main2.sh` o `build-parallel.sh`.

## Dependències

### [fhash](https://github.com/jl2922/fhash/blob/master/fhash.f90)

Implementació de diccionari (_hashmap_) a Fortran 90. És una implementació no genérica i amb clau valor `integer`.

Aquesta dependència s'utilitza pel cas on els nodes de la xarxa no tinguin nombres consecutius.

### mt19937

Generador de números aleatoris