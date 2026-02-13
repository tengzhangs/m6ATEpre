import torch.nn as nn
import torch
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc, precision_score, recall_score, f1_score
import torch.optim as optim
from sklearn.model_selection import KFold
import numpy as np
import pandas as pd
import random
import argparse  
import os       

class MLP(nn.Module):
    def __init__(self, input_dim, num_classes):
        super().__init__()
        self.layers = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Dropout(0.1)
        )
        self.output = nn.Linear(32, num_classes)  # 分离输出层
    
    def forward(self, x):
        features = self.layers(x)
        logits = self.output(features)
        return torch.log_softmax(logits, dim=1)

def train_mlp_kfold(
    feat_path,          # embedding features' path
    label_path,         # samples label's path
    output_path,        # output path
    seed=42,         
    n_splits=5,        
    epochs=200,        
    lr=0.001,         
    weight_decay=5e-4   
):
    try:
        # read features and label
        feat_df = pd.read_csv(feat_path)  
        label_df = pd.read_csv(label_path)
  
        candidate_gene_feat = feat_df.copy()
        candidate_gene_feat['label'] = label_df.values  
        
        shuffled_feat = candidate_gene_feat.sample(frac=1, random_state=seed)  
        # setting random seed
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
            torch.backends.cudnn.deterministic = True
        np.random.seed(seed)
        random.seed(seed)
      
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        y = torch.from_numpy(shuffled_feat["label"].to_numpy()).long().to(device)
        # extract features
        node_feature = shuffled_feat.drop(columns=['label'])
        X = torch.tensor(node_feature.values, dtype=torch.float32).to(device)
        input_dim = X.shape[1]
        num_classes = 2
        
        accuracies = []
        auc_value = []
        precision = []
        recall = []
        f1 = []
        pr_auc_value = []
        
        # 5-CV training
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
        final_model = None
        
        for fold, (train_idx, val_idx) in enumerate(kf.split(X)):
            print(f"\n===== Fold {fold+1}/{n_splits} =====")
            X_train, X_val = X[train_idx], X[val_idx]
            y_train, y_val = y[train_idx], y[val_idx]
            
            model = MLP(input_dim=input_dim, num_classes=num_classes).to(device)
            criterion = nn.NLLLoss()
            optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
            
            # training model
            for epoch in range(epochs):
                model.train()
                optimizer.zero_grad()
                outputs = model(X_train)
                loss = criterion(outputs, y_train)
                loss.backward()
                optimizer.step()
            
            # evaluate model
            model.eval()
            with torch.no_grad():
                pred_probs = torch.exp(model(X_val))[:, -1]
                pred_labels = (pred_probs >= 0.5).float()
                
                # caculate evaluted metrics
                accuracy = (pred_labels == y_val.float()).float().mean().item()
                y_true_np = y_val.cpu().numpy()
                y_prob_np = pred_probs.cpu().numpy()
                y_pred_np = pred_labels.cpu().numpy()
                
                valid_auc = roc_auc_score(y_true_np, y_prob_np)
                precs, recalls, _ = precision_recall_curve(y_true_np, y_prob_np)
                pr_auc = auc(recalls, precs)
                prec_value = precision_score(y_true_np, y_pred_np, zero_division=0)
                rec_value = recall_score(y_true_np, y_pred_np, zero_division=0)
                f1_value = f1_score(y_true_np, y_pred_np, zero_division=0)
                
                accuracies.append(accuracy)
                auc_value.append(valid_auc)
                precision.append(prec_value)
                recall.append(rec_value)
                f1.append(f1_value)
                pr_auc_value.append(pr_auc)
                
                print(f"Fold {fold+1} evaluation result：")
                print(f"Accuracy: {accuracy:.4f} | AUC: {valid_auc:.4f} | PRAUC: {pr_auc:.4f}")
                print(f"Precision: {prec_value:.4f} | Recall: {rec_value:.4f} | F1: {f1_value:.4f}")
            
            final_model = model
        
        # caculate the mean metrics
        metrics_df = pd.DataFrame({
            'AUC': [np.mean(auc_value)],
            'PRAUC': [np.mean(pr_auc_value)],
            'Accuracy': [np.mean(accuracies)],
            'Precision': [np.mean(precision)],
            'Recall': [np.mean(recall)],
            'F1 Score': [np.mean(f1)]
        })
        print(metrics_df)
        
        # predictiong all data
        final_model.eval()
        with torch.no_grad():
            full_pred_probs = torch.exp(final_model(X))[:, -1].cpu().numpy()
            full_pred_labels = (full_pred_probs >= 0.5).astype(int)
        
        # ---------------------- 3.save the result to output path ----------------------
        os.makedirs(output_path, exist_ok=True)
        # save 5-cv evaluate result
        metrics_save_path = os.path.join(output_path, "5CV_performance_metrics.csv")
        metrics_df.to_csv(metrics_save_path, index=False)
        # save predicting result
        pred_result_df = feat_df.copy()
        pred_result_df['true_label'] = label_df.values
        pred_result_df['pred_label'] = full_pred_labels
        pred_result_df['pred_prob'] = full_pred_probs
        pred_save_path = os.path.join(output_path, "final_prediction_results.csv")
        pred_result_df.to_csv(pred_save_path, index=False)
        
        print(f"\n5-CV performance evalutation is saved to：{metrics_save_path}")
        print(f"predict result is saved to：{pred_save_path}")
        
        return metrics_df, full_pred_labels, full_pred_probs

# ---------------------- 命令行参数解析 + 函数调用（核心新增） ----------------------
if __name__ == "__main__":
  
    parser = argparse.ArgumentParser(description='MLP 5-CV training and predict for m6ATEpre')
    parser.add_argument('--embeedding_feat', 
                        default="./m6ADF1_feat.csv")
    parser.add_argument('--samples_label', 
                        default="./labels.csv")
    parser.add_argument('--output_path', 
                        default="./output/")
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--n_splits', type=int, default=5)
    parser.add_argument('--epochs', type=int, default=200)
    parser.add_argument('--lr', type=float, default=0.001)
    parser.add_argument('--weight_decay', type=float, default=5e-4)

    args = parser.parse_args()
    
    try:
        train_mlp_kfold(
            feat_path=args.embeedding_feat,
            label_path=args.samples_label,
            output_path=args.output_path,
            seed=args.seed,
            n_splits=args.n_splits,
            epochs=args.epochs,
            lr=args.lr,
            weight_decay=args.weight_decay
        )
