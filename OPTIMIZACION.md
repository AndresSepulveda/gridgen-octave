# gridgen para Octave (prefijo `oct_`) — notas de la optimización

Versión optimizada de `gridgen`, probada con el caso `examples/create_glo_15m.m`
(grilla global de 15', 1441 x 721 = 1.038.961 celdas).

Los resultados son **idénticos** a los del código original, no aproximados: las
máscaras, la batimetría, las obstrucciones y los archivos de texto para
WAVEWATCH III salen bit a bit iguales. Eso se verifica con el script
`examples/oct_verify_gridgen.m` (ver más abajo).

## El prefijo `oct_`

Todas las **funciones** de este paquete llevan el prefijo `oct_`:
`oct_generate_grid`, `oct_clean_mask`, `oct_create_obstr`, y así las 28 (las 26
de `bin/`, más `oct_setup_gridgen` y `examples/oct_verify_gridgen`). El nombre
del archivo y la línea `function` van siempre en par, y todas las llamadas
internas y las de los scripts de ejemplo apuntan a los nombres nuevos.

Esto permite tener **las dos versiones en la ruta al mismo tiempo**, sin que una
tape a la otra: el original sigue respondiendo a `generate_grid` y esta versión
a `oct_generate_grid`. Por eso `oct_verify_gridgen` ya no necesita agregar y
quitar rutas para comparar: elige la versión con un prefijo y `str2func`.

Los **scripts** conservan su nombre, porque no son funciones:
`examples/create_glo_15m.m`, `examples/bench_create_glo_15m.m`,
`examples/create_grid.m`, `examples/create_grid_Reloncavi_reg.m`,
`examples/create_grid_arctic_40km.m`, `examples/maskmod.m` y `create_grid.m`.
Lo que cambió dentro de ellos son las llamadas, que ahora son a `oct_*`. Si
pones los dos directorios `examples/` en la ruta a la vez, `create_glo_15m`
queda ambiguo: conviene tener solo uno, o renombrar también los scripts.

Para instalar las rutas, la función de siempre, con su nombre nuevo:

```matlab
oct_setup_gridgen
```


## Qué se cambió

| Archivo | Cambio | Efecto |
|---|---|---|
| `bin/oct_compute_cellcorner_all.m` | **Nuevo.** Versión vectorizada de `oct_compute_cellcorner`: calcula las esquinas de todas las celdas en una llamada en vez de una por celda. | Elimina ~1 millón de llamadas y el arreglo de estructuras de 1 millón de elementos en `oct_generate_grid` y `oct_create_obstr`. |
| `bin/oct_inpolygon_fast.m` | **Nuevo.** `inpolygon` para polígonos con muchos vértices: agrupa los puntos por latitud y prueba cada grupo solo contra las aristas que cruzan esa latitud, en vez de todos los puntos contra todas las aristas. | Hasta ~900× en los polígonos costeros grandes. Es el cambio que domina `oct_clean_mask`. |
| `bin/oct_inpolygon_box.m` | **Nuevo.** `inpolygon` para cajas rectangulares alineadas a los ejes (las cajas de dominio y las sub-cajas de `oct_split_boundary`): el test se reduce a comparaciones de intervalos. | ~7× en el recorte de polígonos. |
| `bin/oct_generate_grid.m` | El promedio por celda se hace con sumas acumuladas sobre tiras de latitud de la grilla base, procesando una fila de celdas a la vez, en vez de un `inpolygon` por celda. Se omite el `unique`/`intersect` sobre la batimetría completa cuando las longitudes ya son crecientes. | 60–575× según el caso. |
| `bin/oct_clean_mask.m` | Las celdas regulares se agrupan y sus 36 puntos de muestreo se arman por columna/fila; se descartan de antemano los puntos fuera de la caja del polígono; usa `oct_inpolygon_fast`. | 8–33× |
| `bin/oct_remove_lake.m` | El relleno por inundación celda a celda se reemplaza por componentes conexas con matriz dispersa (`dmperm`), con el mismo criterio de numeración de los cuerpos de agua. | 20–300× |
| `bin/oct_create_obstr.m` | Los puntos de costa se asignan a celdas candidatas por binning (una vez), en vez de probar todos los puntos de costa contra cada celda húmeda; la estructura `cell` pasa a ser compacta (solo las celdas con costa) en vez de `Nx*Ny`. | 4–15× |
| `bin/oct_split_boundary.m` | Los trozos se acumulan en un `cell` y se concatenan una sola vez, en vez de hacer crecer un arreglo de estructuras elemento por elemento (copia O(n²)). | 11× |
| `bin/oct_compute_boundary.m` | Se seleccionan por adelantado los polígonos del nivel pedido cuyos límites cruzan el dominio, en vez de recorrer la base completa; usa `oct_inpolygon_box`. | ~1.1× (y hace mucho más rápido `oct_split_boundary`) |
| `bin/oct_write_ww3file.m`, `bin/oct_write_ww3obstr.m` | Un solo `fprintf` por arreglo en vez de uno por fila. | Menor, pero el texto generado es idéntico. |
| `examples/create_glo_15m.m` | Dos correcciones de robustez: `mkdir out_dir` creaba un directorio llamado literalmente `out_dir` (ahora crea el que indica la variable, si no existe), y `Nu` se usaba sin inicializar cuando `opt_poly = 0`. | — |

## Dos errores del código original, corregidos

### 1. `oct_write_ww3meta` abortaba en el caso global

`examples/create_glo_15m.m` llamaba a `write_ww3meta` (el original, sin
prefijo) **sin el argumento `gtype`**, que la función necesita como segundo parámetro:

```matlab
% original (todos los demás ejemplos sí pasan 'RECT')
write_ww3meta([out_dir,'/',fname],lon,lat,1/depth_scale,1/obstr_scale);
```

Así `gtype` recibía el vector `lon`, el `switch` no reconocía el tipo de
grilla, la función imprimía `!!ERROR!!: Unrecognized Grid Type`, cerraba el
archivo y retornaba. El resultado es que **`glo_15m.meta` salía solo con el
encabezado de comentarios, sin ningún dato numérico**: ni dimensiones, ni
resolución, ni origen, ni los registros de los archivos. Es el archivo que
alimenta `ww3_grid.inp`, de modo que el caso global entregaba un `.meta`
inservible, en silencio y con los otros tres archivos correctos.

Había además un segundo problema en la misma llamada: con `'RECT'` la función
toma `Nx`, `Ny`, `dx` y `dy` de las matrices de coordenadas indexándolas en 2D
(`size(lon)`, `lon(1,2)-lon(1,1)`, `lat(2,1)-lat(1,1)`), así que necesita las
matrices de `meshgrid`, no los ejes 1D. Los otros ejemplos le pasan las 2D; en
`create_glo_15m.m` `lon` y `lat` son los ejes 1D.

Corregido en `examples/create_glo_15m.m` y en
`examples/bench_create_glo_15m.m`, que ahora reconstruyen las 2D a partir de
los ejes ya recortados:

```matlab
[lon2d,lat2d] = meshgrid(lon,lat);
oct_write_ww3meta([out_dir,'/',fname],'RECT',lon2d,lat2d,1/depth_scale,...
                                1/obstr_scale);
```

El `.meta` del caso global queda con sus registros completos:

```
   'RECT'  FLAGLL CSTRNG
1440 	 721
15.00 	 15.00 	 60.00
  0.0000 	 -90.0000 	  1.00
```

Los demás ejemplos no estaban afectados y no se tocaron.

### 2. `oct_create_obstr` abortaba con dominios enteramente en tierra

Durante las pruebas apareció un error que **también está en el código
original**: cuando el dominio pedido queda completamente dentro de un polígono
de tierra, `oct_compute_boundary` guarda ese trozo (la caja del dominio) en `.x` y
`.y` como vector **fila**, mientras que en todas las demás ramas los guarda como
columna. Después, `oct_create_obstr` concatena esos campos y aborta:

```
error: vertical dimensions mismatch (1x5 vs 40x1)
error: called from
    create_obstr at line 105 column 4
```

Es reproducible con el código original sobre cualquier caja que caiga entera en
tierra. En la versión optimizada se corrigió en el origen (`oct_compute_boundary`
guarda columnas, como en el resto de las ramas) y además `oct_create_obstr` acepta
ambas orientaciones.

Consecuencia para la verificación: en los casos donde aparece ese trozo, el
campo `.x`/`.y` sale como columna en vez de fila. Los **valores** son los
mismos, y en una prueba con 4.463 trozos, 32 difieren solo en orientación y
ninguno en valores ni en metadatos. Como el original aborta en ese escenario, no
hay ningún caso que antes funcionara y ahora dé otro resultado.

## Rutas de respaldo

Las rutinas optimizadas conservan el algoritmo original como respaldo y lo usan
automáticamente cuando sus supuestos no se cumplen (grillas curvilíneas, celdas
que cruzan el meridiano 0/360, grillas de menos de 3 celdas por lado, máscaras
que no son 0/1). Los originales quedan como `bin/oct_clean_mask_orig.m`,
`bin/oct_create_obstr_orig.m` y `bin/oct_remove_lake_orig.m`; `oct_generate_grid.m` lleva la
ruta original incorporada y la aplica celda por celda a las que no entran en la
ruta rápida.

Un detalle importante en `oct_generate_grid`: cuando una arista de celda cae
numéricamente sobre un nodo de la grilla base, qué celda se queda ese nodo lo
decide el redondeo dentro de `inpolygon` (y su tolerancia, que depende de la
versión de MATLAB/Octave). Esas columnas y filas se dejan a propósito en la ruta
original, para que el resultado sea exactamente el que da `inpolygon`. Por eso el
caso global de 15' (aristas en minutos exactos) tiene 1442 celdas que siguen el
camino lento.

## Verificación

```matlab
cd examples
% oct_verify_gridgen(dir_original, dir_de_esta_version, dir_datos_referencia, caja, dx)
oct_verify_gridgen('../../gridgen-original','..','../reference_data',[-4 104 2 110],0.25)
```

Corre las mismas etapas del caso global sobre una región chica con ambas
versiones, informa el tiempo de cada etapa y verifica que todo salga idéntico
(batimetría, las tres máscaras, el mapa de cuerpos de agua, `sx`, `sy` y los
archivos de texto). La caja es `[lat_min lon_min lat_max lon_max]` con longitud
en 0..360. Conviene mantenerla chica: el código original es lento.

## Medición de tiempos

```matlab
cd examples
bench_create_glo_15m      % mismo cálculo que create_glo_15m.m, con tiempos por etapa
```

## Resultados medidos

Las mediciones se hicieron en **Octave 10.3** sobre 2 núcleos, con una
batimetría base sintética de 1' (21601 x 10801) y una base costera sintética con
la misma estructura que GSHHS (21.540 polígonos, 3,6 millones de puntos). **No
son tiempos de MATLAB con los datos reales**; sirven para comparar las dos
versiones entre sí, no para predecir tu máquina. Las razones entre versiones
deberían trasladarse; los valores absolutos, no.

Pipeline completo, global a 2.25° (13.041 celdas) — original contra optimizado,
mismos datos de entrada:

| Etapa | Original | Optimizado | Aceleración |
|---|---|---|---|
| `oct_generate_grid` | 99.1 s | 4.5 s | 22× |
| `oct_compute_boundary` | 0.24 s | 0.13 s | 1.8× |
| `oct_split_boundary` | 0.20 s | 0.14 s | 1.4× |
| `oct_clean_mask` | 7.24 s | 1.05 s | 6.9× |
| `oct_remove_lake` | 1.01 s | 0.00 s | ~300× |
| `oct_create_obstr` | 85.9 s | 0.98 s | 88× |
| escritura | 0.02 s | 0.02 s | 1× |
| **TOTAL** | **193.6 s** | **6.8 s** | **28×** |

Todas las salidas idénticas.

Caso global completo a 15' (1.038.961 celdas, 20.098 polígonos en el dominio),
versión optimizada:

| Etapa | Optimizado |
|---|---|
| Lectura de costas | 0.6 s |
| `oct_generate_grid` | 13.7 s |
| `oct_compute_boundary` | 34.4 s |
| `oct_split_boundary` | 28.6 s |
| `oct_clean_mask` | 96.7 s |
| `oct_remove_lake` | 0.3 s |
| `oct_create_obstr` | 71.2 s |
| Escritura | 1.0 s |
| **TOTAL** | **246.6 s** |

El original no se corrió completo en este caso: solo `oct_generate_grid` tomaría del
orden de 2 horas con esta grilla, y el total estaría en varias horas. Esa cifra
es una extrapolación, no una medición.

## Pruebas de equivalencia

Cada rutina se comparó con la original sobre casos construidos para ejercitar
las ramas y los casos límite. Todas dieron resultados idénticos:

- `oct_compute_cellcorner_all`: 1.050.652 celdas en 5 grillas (global 15', dx≠dy, no
  uniforme, curvilínea, con saltos de longitud) — 0 diferencias.
- `oct_generate_grid`: 7 grillas (regional, aristas sobre nodos de la base, dx≠dy, no
  uniforme, más fina que la base, global 2° con celdas que envuelven y polos,
  global 1°, mixta) — 0 diferencias.
- `oct_clean_mask`: 2 regiones (dx=dy y dx≠dy), máscaras intermedia y final.
- `oct_remove_lake`: 9 máscaras (máscara real, aleatorias de 3 densidades con y sin
  envoltura en x, con `lake_tol` positivo y negativo, manchas suaves).
- `oct_create_obstr`: 2 regiones reales y 5 pruebas de carga con costas densas
  (incluyendo `offset` 2, dx≠dy y grilla no uniforme).
- `oct_compute_boundary`: 8 cajas, incluyendo las dos mitades del globo, el globo
  completo, nivel 2 (lagos) y los polígonos de usuario.
- `oct_split_boundary`: 4.082 polígonos de entrada → 4.463 trozos.
- `oct_inpolygon_fast`: 480.000 puntos contra los 12 polígonos costeros más grandes,
  más 6 pruebas adversarias con puntos exactamente sobre aristas, sobre vértices
  y en una retícula de 0.25° — 0 diferencias en `in` y en `on`.
- `oct_inpolygon_box`: 1.174.760 pruebas de punto, incluyendo esquinas exactas,
  aristas exactas, desplazamientos de ±4 eps y polígonos que no son cajas (que
  deben delegar a la rutina original) — 0 diferencias.

## Qué no se cambió

- `oct_compute_boundary` sigue siendo cell-by-cell en su lógica de cruces de
  dominio; es la etapa que menos ganó (~1.1×) y en el caso global pesa 34 s.
  Vectorizarla implicaría reescribir el seguimiento de segmentos, con riesgo
  alto de cambiar resultados en casos límite.
- Los bucles de `oct_create_obstr` que fusionan segmentos entre celdas vecinas
  (sombras en x e y) se mantienen tal cual: operan sobre arreglos muy chicos por
  celda y su semántica depende del orden.
- La sección de visualización de `create_glo_15m.m` quedó intacta.
