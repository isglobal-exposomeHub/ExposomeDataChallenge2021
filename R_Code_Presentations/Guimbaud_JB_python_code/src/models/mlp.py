# Author Jean-Baptiste Guimbaud
# Meersens

import tensorflow as tf


class MLP(tf.keras.Model):
    def __init__(self,
                 layers_size,
                 num_features,
                 dropout=0.2,
                 regularizer=None,
                 activation="elu",
                 classifier=False):
        num_layers = len(layers_size)
        assert num_layers > 1
        super(MLP, self).__init__()
        self.classifier = classifier
        self.num_layers = num_layers
        self.hidden_layers = []
        # print("num features:", num_features)

        # input layer
        self.hidden_layers.append((
            tf.keras.layers.Dense(layers_size[0],
                                  input_shape=(num_features,),
                                  activation=activation,
                                  kernel_regularizer=regularizer,
                                  name="dense_0"),
            tf.keras.layers.Dropout(dropout,
                                    name="dropout_0"))
        )

        # hidden_layers
        for i in range(1, num_layers):
            self.hidden_layers.append((
                tf.keras.layers.Dense(layers_size[i], 
                                      activation=activation,
                                      kernel_regularizer=regularizer, 
                                      name=f"dense_{i}"),
                tf.keras.layers.Dropout(dropout, name=f"dropout_{i}"))
            )

        # output layer
        if classifier:
            self.output_layer = tf.keras.layers.Dense(4, activation='softmax')
        else:
            self.output_layer = tf.keras.layers.Dense(1, activation='linear')

    def call(self, x, training=False):
        for layer in self.hidden_layers:
            x = layer[0](x)
            x = layer[1](x, training=training)

        x = self.output_layer(x)

        if self.classifier:
            return x
        return x[:, 0]