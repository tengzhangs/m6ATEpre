# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 16:07:37 2025
@author: tengz
"""
import numpy as np
import pandas as pd
from keras.layers import Input, Dense
from keras.models import Model

def get_embedding_by_autoencoder(data_np, encoding_dim=100, epochs=100, batch_size=32, validation_split=0.3):
    """
    # 1. Input dimension
    input_dim = data_np.shape[1]
    
    # 2. Construct autoencoder
    input_layer = Input(shape=(input_dim,))
    encoded = Dense(encoding_dim, activation='relu')(input_layer)  # encoder layer(generate embedding)
    decoded = Dense(input_dim, activation='sigmoid')(encoded)      # decoder layer
    
    autoencoder = Model(input_layer, decoded)  # training model
    encoder = Model(input_layer, encoded)      # encoder（etract embedding）
    
    # 3. Traing autoencoder
    autoencoder.compile(optimizer='adam', loss='mean_squared_error')
    autoencoder.fit(
        data_np, data_np,
        epochs=epochs,
        batch_size=batch_size,
        shuffle=True,
        validation_split=validation_split
    )
    
    # 4. return embedding
    encoded_data = encoder.predict(data_np)
    return encoded_data


# 1. read  samples
pos_data = pd.read_csv(
    "F:\\m6A_translation\\dataset\\HeLa_cell\\data\\single_base_level\\new_pos_neg_samples\\cdhit_proces\\pos_motif_vec.csv",
    header=None
)
pos_datas = pos_data.to_numpy()  

neg_data = pd.read_csv(
    "F:\\m6A_translation\\dataset\\HeLa_cell\\data\\single_base_level\\new_pos_neg_samples\\cdhit_proces\\neg_motif_vec.csv",
    header=None
)
neg_datas = neg_data.to_numpy()  
# 2. obtain embedding for postive samples
pos_encoded_data = get_embedding_by_autoencoder(pos_datas)  # embedding for postive samples

neg_encoded_data = get_embedding_by_autoencoder(neg_datas)  # embedding for negative samples

# 3. save embedding for postive and negative samples
# save embedding for postive
pos_motif_encoded = pd.DataFrame(pos_encoded_data)
pos_motif_encoded.to_csv(
    'F:\\m6A_translation\\dataset\\HeLa_cell\\data\\single_base_level\\cdhit_proces\\pos_motif_encoded.csv',
    index=False,
    header=False
)

# save embedding for negative
neg_motif_encoded = pd.DataFrame(neg_encoded_data)
neg_motif_encoded.to_csv(
    'F:\\m6A_translation\\dataset\\HeLa_cell\\data\\single_base_level\\cdhit_proces\\neg_motif_encoded.csv',
    index=False,
    header=False
)
