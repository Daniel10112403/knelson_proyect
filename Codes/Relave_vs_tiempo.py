import pandas as pd
import matplotlib.pyplot as plt
#Función para determinar la respectiva gráfica
def graf_relvstime(Archive):
# Lectura de datos
       data = pd.read_csv(Archive)
# content = (Atomic Absorption result [mg/L] * 1g / 1000 mg * 1 L/1000 mL * 10mL) / (sample mass g * 1 ton /1e6 g )
# Extraer datos
       gold_law = [i / 3 for i in list(data.loc[:, 'Concentration'])]
       dic = {'Time (min)': data.loc[:, 'Time'], 'Concentration_value (mg/L)': data.loc[:, 'Concentration'],
              'Gold_law (g/ton)': gold_law}
# Generar nuevo set de datos
       set_datos = pd.DataFrame(dic)
       return set_datos


set_datos_BP = graf_relvstime('Concentration_Relave_vs_tiempo_rawdata_BP.csv')
set_datos_P6 = graf_relvstime('Concentration_Relave_vs_tiempo_rawdata_P6.csv')
#Graficar resultados obtenidos
figure = plt.figure(figsize=[10,10])
ax = plt.axes()
ax.plot(set_datos_BP['Time (min)'], set_datos_BP['Gold_law (g/ton)'], label= 'BP')
ax.plot(set_datos_P6['Time (min)'], set_datos_P6['Gold_law (g/ton)'], label= 'P6')
ax.set_xlabel('Tiempo (min)')
ax.set_ylabel('Contenido de Oro (g/ton)')
ax.legend(loc='best')
ax.annotate('Para BP se realizó un corte en el caudal de fluidización\n'
            'que posteriormente se reinició',
            xy=(0.1, .8), textcoords='axes fraction ')
plt.grid()
plt.show()