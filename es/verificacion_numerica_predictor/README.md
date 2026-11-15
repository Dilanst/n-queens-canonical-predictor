# Reproducción numérica del predictor canónico

El programa calcula desde cero el máximo canónico para diez tamaños cuadrados y resuelve por separado las ecuaciones del perfil de masa uno y la discretización de Nyström. Los enteros de referencia del conteo proceden de OEIS A000170. La referencia de actividad finita se conserva en un archivo distinto de la referencia canónica para n=5.

## Ejecutar

Se requieren Node.js 20 o posterior y Python 3 con NumPy. Instale las dependencias en su entorno con `python -m pip install -r requirements.txt`. Desde esta carpeta, `node reproducir.mjs` crea un directorio nuevo con todas las salidas y los tiempos. En Unix también se puede ejecutar `sh reproducir.sh`. La variable de entorno `PYTHON` permite elegir el intérprete. Se utiliza CPU. Los resultados se calculan en doble precisión.

Para verificar las cifras que se entregan sin repetir los solucionadores, use `node comprobar.mjs`. Se comprueban individualmente los valores de las dos tablas, los ajustes, ambas paridades de la cuadratura, la precisión de los decimales del resumen y la introducción y los controles empíricos de la discusión. El programa devuelve un código de error si encuentra una discrepancia. La reproducción integral ejecuta también `barrera_4_5.py`, conserva su salida y comprueba sus dos transiciones de borde. El código encuentra por separado el valor canónico tabulado para n=5 en su referencia de masa uno.

## Tiempos medidos

Tiempos medidos en el entorno JavaScript V8 de la revisión, con CPU, redondeados a milisegundos. El tiempo de Python depende del intérprete y se escribe en el archivo de tiempos al ejecutar el reproductor integral.

| Cálculo | Casos | Segundos |
|---|---:|---:|
| Máximo canónico, pseudodeterminante y brecha | 10 | 3,251 |
| Perfil y Nyström | 50 pares de mallas | 19,541 |
| Actividad finita | 13 × 61 | 2,506 |
| Barrera alternativa con simetrías, JavaScript | 3 × 41 | 0,182 |

## Resultados y archivos

`canonico.csv` conserva los extremos, los residuos, las iteraciones y las brechas antes del redondeo. Los diez valores de log Q_pred coinciden a siete decimales con la tabla de referencia, y las brechas a cinco. El test de la brecha usa iteración inversa con residuo y Cholesky desplazado sobre el complemento del modo nulo. Las verificaciones de definitud e interioridad son numéricas.

| n | log Q_pred | min p | max p | 1/4 − max p | min h | brecha |
|---:|---:|---:|---:|---:|---:|---:|
| 20 | 24.3667849 | 0.02473820 | 0.06915949 | 0.18084051 | 0.17854349 | 0.14177800 |
| 27 | 39.9860346 | 0.01787269 | 0.05152626 | 0.19847374 | 0.17749523 | 0.14130756 |
| 40 | 73.4844735 | 0.01177350 | 0.03486465 | 0.21513535 | 0.17636131 | 0.14083134 |
| 60 | 132.9170795 | 0.00771761 | 0.02329709 | 0.22670291 | 0.17560694 | 0.14047948 |
| 80 | 199.0846876 | 0.00573905 | 0.01749172 | 0.23250828 | 0.17522429 | 0.14029676 |
| 100 | 270.2725639 | 0.00456768 | 0.01400209 | 0.23599791 | 0.17499290 | 0.14018495 |
| 130 | 384.3475269 | 0.00349688 | 0.01077685 | 0.23922315 | 0.17477810 | 0.14008030 |
| 160 | 505.3805043 | 0.00283273 | 0.00875917 | 0.24124083 | 0.17464325 | 0.14001418 |
| 200 | 675.3684475 | 0.00226030 | 0.00700937 | 0.24299063 | 0.17452600 | 0.13995644 |
| 240 | 853.3899390 | 0.00188032 | 0.00584225 | 0.24415775 | 0.17444764 | 0.13991772 |

La referencia central es α₀=1,944001. Las dos columnas extremas siguientes usan respectivamente 1,944000752 y 1,944001082. Los ajustes emplean mínimos cuadrados por transformación de Householder. Los extremos agregados impresos en el artículo están redondeados.

| Rango n | Grado en 1/n | C con α inferior | C con α₀ | C con α superior |
|---|---:|---:|---:|---:|
| 40–130 | 2 | 0.94011554 | 0.94017293 | 0.94019191 |
| 40–130 | 3 | 0.94005476 | 0.94013367 | 0.94015977 |
| 40–240 | 2 | 0.94007420 | 0.94015487 | 0.94018155 |
| 40–240 | 3 | 0.94002363 | 0.94013272 | 0.94016878 |
| 60–240 | 2 | 0.94004855 | 0.94014309 | 0.94017435 |
| 60–240 | 3 | 0.94000433 | 0.94013200 | 0.94017421 |
| 20–240 | 2 | 0.94036043 | 0.94042196 | 0.94044231 |
| 20–240 | 3 | 0.93968175 | 0.93976215 | 0.93978873 |

`ny_resultados.csv` contiene las 50 mallas. En las parejas 192/384 y 193/385 se obtiene respectivamente 0,94014395698 y 0,94014395711. Los factores nodales del perfil en m=385 están en `perfil_m385.json`. La extrapolación de I₂ y log det₂ a partir de la pareja par, con exponenciación posterior de log det₂, da I₂≈1,02638534, det₂≈0,005162694. La diferencia entre las dos paridades es una medida de estabilidad del cálculo. La estimación carece de cota de error certificada.

`referencias/anunciados.json` guarda el valor de Q(28) publicado en noviembre de 2026 (Cherepanov, Zenodo; OEIS A000170). No entra en las tablas ni en la figura: sólo sirve para comprobar la cifra informativa C'_0(28) de la Sección 9. `actividades_793.csv` documenta la malla de actividad. `referencias/barrera_archivada.txt` contiene la salida histórica del programa Python. `barrera_alternativa.json` registra los 123 puntos de una implementación independiente en JavaScript: en n=4 y 5 el máximo alcanza el umbral 0,24999 por primera vez en Λ=100 y Λ=10 000, respectivamente. En n=6, en Λ=10⁹, el máximo vale 0,2165479. El umbral pertenece a la malla ensayada. Estos resultados a actividad finita y el valor canónico estimado para n=5 tienen interpretaciones distintas.

`comprobar.mjs` compara los valores esperados con el texto LaTeX de las secciones que contienen resultados experimentales. La primera comprobación incluida con la entrega utiliza la salida archivada de Python para la barrera y el cálculo canónico nuevo para la Tabla 1. El reproductor integral sustituye el diagnóstico de actividad finita por los valores de la nueva ejecución Python.

Los códigos se distribuyen bajo la licencia MIT. La información de autoría figura en CITATION.cff.
