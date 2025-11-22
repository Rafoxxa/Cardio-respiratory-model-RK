# Modelo computacional de control cardiorrespiratorio: Simulación y ajuste

Este repositorio contiene archivos para ejecutar la simulación de un modelo de control cardiorrespiratorio y su ajuste con datos experimentales. Los archivos se dividen en cuatro categorías: simulación, ajuste, procesamiento de datos y análisis de resultados.

---

## Simulación

Archivos necesarios:

- **`ForwardModel`**: Archivo principal de simulación. Contiene:
  - **`run_ode`**: Driver del solver, que ejecuta:
    - **`model_basic`**: Contiene todas las ecuaciones del modelo utilizadas por el solver

### Configuración importante

Para ejecutar `ForwardModel`, es fundamental configurar correctamente los hiperparámetros en las funciones `set_up`. 

**Ejemplo de configuración por defecto:**
```matlab
[setup] = set_up("simulation", patient_idx, state, "mix", "dt", 0.1, ...
    'pars_from_fitting', 1, ...
    'fitting_mat_file', {'07-09-2025', '30-09-2025', '20-10-2025', '24-10-2025', '02-11-2025'}, ...
    'time_from_data', 1);
```

**Parámetros clave:**
- `'pars_from_fitting'`: Permite cargar parámetros desde archivos de ajuste
- `'fitting_mat_file'`: Fechas de los ajustes desde donde se cargarán los parámetros
- `'time_from_data'`: Permite tomar tiempos respiratorios desde los datos experimentales

---

## Análisis de sensibilidad

Utiliza la función **`parameter_analysis_fun`** con los siguientes argumentos:

- `patient_number`: char ('1', '2', '3', '4')
- `pindex`: char - índice del parámetro a perturbar
- `load_rb`: char ('1', '0') - '1' para cargar `results_base` existente, '0' para computarlo y guardarlo como nuevo

Esta función llama a:
- **`sens_functions`**: Driver principal para cálculo de sensibilidades. Ejecuta simulación base, calcula simulaciones perturbadas y sensibilidades, y construye el tensor de sensibilidad (STensor). Llama a:
  - `data_processing`
  - `run_ode`

---

## Ajuste de parámetros

**`ForwardModelFitting`**: Ejecuta el proceso de ajuste. Contiene:

- **`data_preprocessing`**: Preprocesa datos experimentales para el ajuste
- **`obj_fun`**: Función objetivo que calcula el error entre simulación y datos experimentales. Ejecuta:
  - `run_ode`

---

## Procesamiento de datos

**`data_preprocessing`**: Recibe datos originales, los lee y procesa para generar un dataset adecuado para el ajuste. También calcula estimaciones polinómicas para VO2, VCO2 y fiO2:

- **`bestPolynomialFit`**: Calcula la estimación polinómica y permite graficar datos con la función ajustada

---

## Otros archivos

### Datos

- **`estimate_newton_ohm`**: Lee archivo de parámetros y estima parámetros vasculares usando la ecuación de Newton-Ohm
- **`readBeatscopeData20`**: Archivo especializado para leer datos de resultados FINAPRES


### Configuración

- **`vectorize_dicts`**: Transforma las estructuras basadas en diccionarios (solver y modelo) a estructuras basadas en vectores. Utiliza:
  - **`transformation_imported_dict`**: Usa expresiones regulares para actualizar líneas en `run_ode` y `model_basic` a tipo vector
- **`Optimize_percentages`**: Define porcentajes de diferentes condiciones vasculares para estimar parámetros correctamente usando ecuaciones de Newton-Ohm. Ya no es necesario después de su uso inicial


