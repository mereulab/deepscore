import os
import scanpy as sc
import numpy as np
import tensorflow as tf
import pandas as pd
from sklearn.model_selection import train_test_split
from tensorflow import keras
from tensorflow.keras.constraints import max_norm
from tensorflow.keras.regularizers import l1_l2
import imageio


class PlotInference(keras.callbacks.Callback):
    
    def __init__(self, ds, adata, gif_name='mygif'):
        super(PlotInference, self).__init__()
        self.adata = adata
        self.filenames = []
        self.gif_name = gif_name
        self.ds = ds
        
    def on_epoch_end(self, epoch, logs=None):
        with tf.device('/CPU:0'):
            self.adata = self.ds.annotate(self.adata)
            filename = f'_ds_training{epoch}.png'
            sc.pl.umap(self.adata, color='Deepscore', legend_loc = 'on data',
                      save=filename, show=False)
            self.filenames.append(f'figures/umap{filename}')

    def on_train_end(self, logs=None):
        
        # build gif
        with imageio.get_writer(f'{self.gif_name}.gif', mode='I') as writer:
            for filename in self.filenames:
                print(filename)
                image = imageio.imread(filename)
                writer.append_data(image)

        # Remove files
        for filename in self.filenames:
            os.remove(filename)

            

class DeepScore():
    
    def __init__(self, hidden_nodes, n_features, n_labels, epochs=30, 
                 batch_size=32, activation="relu", dropout=True, 
                 dropout_rate=0.2, batchnorm=True, lr=0.001,
                 weight_reg=True, l1=0, l2=0):
        
        super(DeepScore, self).__init__()
        self.epochs = epochs
        self.batch_size = batch_size
        self.activation = activation
        self.dropout = dropout
        self.dropout_rate = dropout_rate
        self.lr = lr
        
        if weight_reg:
            self.w_norm = max_norm(2.)
        else:
            self.w_norm = None
    
        model = keras.Sequential(name="deepscore")
        model.add(keras.Input(shape=n_features,))
        if batchnorm:
            model.add(keras.layers.BatchNormalization())

        for i, n in enumerate(hidden_nodes):
            model.add(keras.layers.Dense(n, activation=activation, 
                    kernel_constraint=self.w_norm, name=f"dense{n}",
                    kernel_regularizer=l1_l2(l1=l1, l2=l2)))
            if dropout:
                model.add(keras.layers.Dropout(dropout_rate))
            if batchnorm:
                model.add(keras.layers.BatchNormalization())

        model.add(keras.layers.Dense(n_labels, activation='softmax', 
                                     name="output"))

        model.compile(
            optimizer=keras.optimizers.Adam(learning_rate=lr),
            loss=keras.losses.CategoricalCrossentropy(),
            metrics=[keras.metrics.CategoricalAccuracy()])

        model.summary()
        self.model = model
    

            
    def set_reference(self, adata, label_by, test_prop=0.25):
    
        # Turn each cell label into a hot one encoding
        labels = pd.get_dummies(adata.obs[label_by])
        self.label_dict = {i: k for i, k in enumerate(list(labels.columns))}

        xtrain, xtest, ytrain, ytest = train_test_split(adata.X, labels,
                                                        test_size = test_prop)

        xtrain = tf.constant(xtrain)
        xtest = tf.constant(xtest)
        ytrain = tf.constant(ytrain)
        ytest = tf.constant(ytest)

        self.splitted_data = xtrain, xtest, ytrain, ytest
    
    
    
    
    def train(self, val_split=0.1, earlystopping=True, lr_scheduler=False, 
              patience=2, create_gif=False, gif_dataset=None):
    
        xtrain, xtest, ytrain, ytest = self.splitted_data

        my_callbacks = [tf.keras.callbacks.TensorBoard(
            log_dir='./deepscore_logs')]
        
        if earlystopping:
            my_callbacks.append(
                tf.keras.callbacks.EarlyStopping(patience=patience))
        if lr_scheduler:
            my_callbacks.append(
                tf.keras.callbacks.LearningRateScheduler(lr_scheduler))
        if create_gif:
            if gif_dataset is None:
                print('Please specify a dataset in gif_dataset')
            else:
                my_callbacks.append(PlotInference(self, gif_dataset))

        self.model.fit(xtrain, ytrain, epochs=self.epochs, 
                       batch_size=self.batch_size, validation_split=val_split, 
                       verbose = True, callbacks=my_callbacks)

        print("\nEvaluating model performance on unseen data (test data):\n")
        evaluation = self.model.evaluate(xtest, ytest, batch_size=self.batch_size)
        print(f"\ntest loss: {round(evaluation[0], 5)}, test accuracy:' \
              '{round(evaluation[1], 5)}")
    
    
    
    def annotate(self, adata, pred_key='Deepscore'):

        tf.data.Dataset.from_tensor_slices(adata.X)
        test = self.model.predict(adata.X)

        indeces = np.argmax(test, axis=1)
        maxs = np.max(test, axis=1)
        predictions = []

        for m, i in zip(maxs, indeces):
            if m > 0.5:
                predictions.append(self.label_dict[i])
            else:
                predictions.append('Unclassified')

        adata.obs[pred_key] = predictions

        return(adata)        
        

    