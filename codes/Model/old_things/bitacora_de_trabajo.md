Bitácora de trabajo:

# 06-03-2024

## P_sa [PENDIENTE]

Hoy día intentaré seguir corrigiendo el modelo cardiovascular. Las observaciones que encontré el día de ayer fueron:
P_sa está bajando demasiado brusco después de la sístole. Esto se debe primordialmente a que:

$dP_{sa} \sim Q_{lv} - Q_{sa}$

Donde el primero depende netamente de la apertura valvular del ventrículo izquierdo (válvula aórtica). El flujo de _lv_ solo ocurrirá si la presión dentro del ventrículo supera a la presión arterial sistémica ($P_{sa}$). Y el flujo $Q_{sa}$ depende de sí mismo y de la diferencia de presión entre la circulación periférica y la arterial.

## V_total_s_v y Pabd [CORREGIDOS]

Tuve problemas con las presiones venosas. En particular con la presión esplácnica (del sistema digestivo), y encontré que el __error__ se debía a que esta dependía de la presión abdominal (Pabd), que estaba comportándose de manera extraña, debido a que el volumen que estaba recibiendo como entrada no era el volumen tidal. Esto ya está __corregido__.

## Cambios

* Algunas cosas que cambié por el momento fueron las __inertancias__, que las estoy considerando unitarias ya que inestabilizan demasiado el sistema, lo que hace inviable correrlo hasta que este llegue a sus condiciones de equilibrio.
* También detecté que el volumen total de la arteria pulmonar partía en 0, debido a que por ser arteria no tiene un volumen no estresado (ya que todos los volumenes totales los estoy definiendo en base a los Vu).

# 07-03-2024

## Simulación respirtaria

El controlador respiratorio subió levemente la fercuencia respiratoria frente al aumento de consumo metabólico. Pero por alguna razón cuando este input entra al sistema no, el controlador óptimo no logra disminuir los valores de las variables.

## P_sa

Me di cuenta que el problema está netamente en P_sp. Por ende es ahi donde tenenmos que picar. La profe también mencionar que R_sa podría estar afectando pero dado que R_sa vale 0.06, significa que habría que variarlo casi 100 veces respecto al valor reportado en literatura como para que su efecto cambie en la edo.

# 08-03-24

Al revisar las fuentes de error me di cuenta que el problema solamente puede estar en:

* P: Parámetros
* C: Condiciones iniciales
* I: Implementación de las ecuaciones
* E: Errores estructurales de las ecuaciones (mal modelamiento o typos de los autores)

Por lo tanto tendré que plantear errores para cada una de estas fuentes de error.
P: Creo que las resistencias que son moduladas por el controlador podrían estar afectando, dado que no conozco su valor correcto. De todas formas al variarlas no logro obtener la presión arterial que deseo. Traté de variar las resistencias periféricas y no se logró un cambio significativo hasta que se llevaron todas a 1.
El otro parámetro que afecté fue la resistencia de la arteria (la aumenté) para que se generara un peso contra el mismo Q_sa. También cambié R_sa a 1.
También multipliqué por un factor el sum(Q_p) para ver si efectivamente eso generaba un cambio en P_sa. Al menos no cae a 0, pero cae como a 20 mmHg, eso está bastante extraño todavía
C: Cambié un poco las condiciones iniciales de Q_sa, pero no hubo mucha mejora
I: Creo que esto esta bien
E: No he encontrado nada por acá

Ahora la película es que la presión la puedo controlar o con los valores de P_sp y por ende de Q_p, desde las condiciones iniciales de Q_sa o desde los valores de Q_lv. Lo extraño es que a pesar incluso de que P_sp se haga pequeño aun así P_sa baja de forma muy radical. Estoy empezando a pensar que sí podría tener que ver con el flujo ventricular, tal vez este cae MUUY RÁPIDO, o tal vez los volúmenes que le llegan son muy pequeños o se vacía muy rápido. No hay chance de hacer que ese flujo sea más lento, porque hay un factor de bibliografía que hace muuuy chica la resistencia de la válvula aórtica.

Otra cosa que me di cuenta, es que si acortamos el ciclo cardíaco, la respuesta mejora significativamente. Creo que al final puede ser que todo tenga que ver más con la duración diástole-sístole y duración del ciclo cardíaco completo. Al menos de esta forma podríamos intentar dejar el sistema fisiológicamente un poco bradicárdico pero con el control natural podrían establecerse los equilibrios de flujo, presión y volumenes.

Ya, cambiando Tsys y dejando las resistencias periféricas como antes se llega a esto:

![1709913088829](image/bitacora_de_trabajo/1709913088829.png)

Lo que implica que a pesar de que no se llega a valores fisiológicos, al menos no estamos en 0.

Después hablé con Martin y me comentó que la principal forma de revisar estos problemas es chequeando los volumenes.

# 09-03-24

No tengo el problema resuelto al 100% pero al menos la presión arterial baja a 65 ml, no es bueno pero no es horrible. Al menos eso se podrá corregir una vez que el control haga su efecto en las resistencias periféricas y en los volumenes venosos no estresados y Tsys. De momento lo que sí hice fue chequear los volumenes y hallé que había una mala implementación en el cálculo de los volumenes periféricos, ya lo corregí. Además revisando un poco de fisiología

![1709983410372](image/bitacora_de_trabajo/1709983410372.png)

Y conversando con Martin, claramente los volumenes venosos TIENEN que ser mayores que los periféricos (arteriales periféricos), ya que los primeros SON por defecto los acumuladores de sangre del cuerpo. La única excepción podrían ser los volumenes venosos del músculo activo dado que ellos están en reposo. Y eso intenté arreglar un poco.

nota* Igual hay algo extraño con los volumenes venosos e_v (como que hay una tendencia a bajar, pero se hacen más chicos que los periféricos, cuando el volumen venosos e)

Para el lunes debería intentar:

- Implementar el control cardiovascular para mandárselo a la profe
- Hacer una documentación apropiada del código
- Revisar si el cardiovscular se puede pasar un ciclo respiratorio
- (para más adelante): chequear J en respiratorio para entender por qué NO baja


# 20/04/24

Últimamente he estado corriendo el modelo en steady state, intentando llegar a las siguientes gráficas

![1713650870952](image/bitacora_de_trabajo/1713650870952.png)

Lo he estado corriendo en local (sin control respiratorio) y en el servidor (con control respiratorio). Con el detalle de que lo estoy haciendo con el consumo en pasos de a 0.8 lt/min, lo que no es correcto, ya que como se ve en las gráficas debería ser desde 0.3 lt/min a 1 lt/min. En local no he obtenido algo muy representativo, pero al menos las presiones tienen algo de sentido a excepción de las diastólicas.

Dejo el resultado de la simulación en el server (con control respiratorio):

![1713651143805](image/bitacora_de_trabajo/1713651143805.png)

(hay que leer las variables en el orden de los plots: dVE, VT, BF, TI, PACO2, HR, PS, PM, PD, PAO2).

Al menos el PAO2 tanto como el PACO2 muestran una tendencia similar pero con ciertos valores un poco extremos (al comienzo). El VT está bajo, pero este depende netamente del valor a1 y a2 de Pmusc, lo que es altamente dependiente de lo que el optimizador encuentra.

Luego se corrieron las mismas simulaciones para steady state en local con las correcciones de las ecuaciones del sistema respiratorio (corrigiendo principalmente Heart Rate y el delay variable debido al flujo Qla). Se observaron cambios principalmente en la parte vascular del steady state (PS, PD y HR cercanas a los valores de interés, aunque HR estaba más cercano a 50 bmp).

Con esto las principales interrogantes que prevalecen son:

* Dónde está el problema en que el sistema no asimile la cantidad que se está consumiendo/produciendo de gas (la solución debería estar en las **ecuaciones de intercambio gasesoso**, **parámetros de control** y en el siguiente punto que es el **optimizador**).
* Por qué el optimizador no logra subir el a1 o a2 cuando TI baja. (la solución a esto debería residir en el **optimizador**, y en el espacio de búsqueda o en los **límites de las variables**, y también en **dVE**. **J** ya está corregido.  Hay que ahondar en ello)

Algo intereseante que se observó hace poco fue que la ecuación de disociación la estaba tomando desde la presión parcial arterial de O2 y CO2, cuando el paper la toma de PA. Así que revisar bien esas ecuaciones.
