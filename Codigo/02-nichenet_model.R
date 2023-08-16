
## Paquetes necesarios
library(OmnipathR)
library(nichenetr)
library(mlrMBO)
library(parallelMap)
library(igraph)

## NicheNet
# Testear el pipeline para comprobar que todo esté correcto
require(OmnipathR)
nichenet_workarounds()
nichenet_test()

# NicheNet requiere tres tipos de redes de interacción. La función nichenet_networks se encarga de la construcción de las tres
# La opción only_omnipath = FALSE hace que no coja solo los datos de Omnipath
# Se obtiene una lista con 3 dataframes que se guarda de forma automática como RDS en una carpeta llamada nichenet_results
networks <- nichenet_networks(only_omnipath = FALSE)

# Cargar información de experimentos de ligandos perturbados contenida en NicheNet y la red de interacción ligando-receptor
# Eliminar los ligandos que no se encentren en la red de interación
expression <- nichenet_expression_data()
lr_network <- networks$lr_network
expression <- nichenet_remove_orphan_ligands(
  expression = expression,
  lr_network = lr_network
  )

# Optimización del modelo
optimization_results <- nichenet_optimization(
    networks = networks,
    expression = expression,
    mlrmbo_optimization_param = list(ncores = 4)
    )

# Construcción del modelo, el resultado se guarda de forma automática como RDS
nichenet_model <- nichenet_build_model(
  optimization_results = optimization_results,
  networks = networks
  )

# Generar la matriz ligando-target
lt_matrix <- nichenet_ligand_target_matrix(
  nichenet_model$weighted_networks,
  networks$lr_network,
  nichenet_model$optimized_parameters
  )
