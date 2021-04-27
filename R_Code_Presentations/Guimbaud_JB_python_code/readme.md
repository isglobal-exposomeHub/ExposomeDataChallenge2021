
# Leveraging machine learning and explainable AI to better understand exposomic data

###  Exposome Data Challenge 2021

  
  

Show command line arguments:

```sh

python main.py [--cross] [-c] [-h]

```
With no argument given the programme will predict IQ in child with a random forest. Then it will use SHAP to explain the predictions.

Optional arguments:
| <!-- --> |<!-- --> | <!-- --> | 
|-|-|-|
|  |  *-*-cross | Compute score with cross validation and exit. | 
| -c |  *-*-correlation | Print the correlation matrix |
| -h |  *-*-help | Show this help message and exit. |

You can easily change parameters like the outcome, model, feature reduction method or input data in the main.py script.

To generate preprocessed data run:
```sh

python preprocessing.py

```

The jupyter notebook files ([hp_tuning_mlp.ipynb](hp_tuning_mlp.ipynb) and [hp_tuning_xgb.ipynb](hp_tuning_xgb.ipynb)) were used for hyperparameter tuning. Hyperparameters are written in the config folder [here](config/).

*Written by Jean-Baptiste Guimbaud, currently data engineer at Meersens.*