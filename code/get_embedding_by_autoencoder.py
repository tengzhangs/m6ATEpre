# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 16:07:37 2025
@author: tengz
"""
import numpy as np
import pandas as pd
from keras.layers import Input, Dense
from keras.models import Model
import argparse  # 仅新增：用于命令行传参
import os        # 仅新增：用于创建输出目录

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

if __name__ == "__main__":
 
    parser = argparse.ArgumentParser()
    parser.add_argument('--pos_path', default='./single_base_level/cdhit_proces/pos_motif_vec.csv')
    parser.add_argument('--neg_path', default='./single_base_level/cdhit_proces/neg_motif_vec.csv')
    parser.add_argument('--output_path', default='./embedding_output/')  
    args = parser.parse_args()

    pos_data = pd.read_csv(args.pos_path, header=None)
    pos_datas = pos_data.to_numpy()  

    neg_data = pd.read_csv(args.neg_path, header=None)
    neg_datas = neg_data.to_numpy()  

    pos_encoded_data = get_embedding_by_autoencoder(pos_datas)
    neg_encoded_data = get_embedding_by_autoencoder(neg_datas)

    pos_motif_encoded = pd.DataFrame(pos_encoded_data)
    neg_motif_encoded = pd.DataFrame(neg_encoded_data)
  
    os.makedirs(args.output_path, exist_ok=True) 
    combined_emb = np.vstack([pos_encoded_data, neg_encoded_data])
 
    combined_csv_path = os.path.join(args.output_path, 'embedding.csv')
    pd.DataFrame(combined_emb).to_csv(combined_csv_path, index=False, header=False)
