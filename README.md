# pyDynamics

PyDy es un software para la simulación de la dinámica no lineal de sistemas discretos. El software está escrito en Python y es utilizado por los investigadores del Grupo de Investigación en Multifísica Aplicada para la simulación de recolectores de energía. El software lee un archivo de entrada en formato ASCII que contiene los parametros del sistema dinámico. 

Hasta el momento el software permite la simulación de 6 sistemas dinámicos: vehiculo de 7 DOF, vehiculo de 5 DOF con dinamica axial, recolectores de energía magnéticos levitantes, recolectores de energía magnéticos sin levitación y sistemas genericos de 1 GDL. El software permite el modelado de sistemas con no linealidades de cualquier tipo. El software procesa la clase correspondiente al sistema a simular e integra las ecuaciones de movimiento mediante el metodo de Runge-Kutta de cuarto orden. 

Esta nueva versin de PyDy que admite sistemas dinámicos que experimenten movimientos de cuerpo rígido arbitrarios. Esto permitela utilización del software para simular recolectores de energía que operan acoplados a palas de aerogeneradores.
