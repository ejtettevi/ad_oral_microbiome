# microbiome_ml.py
import pandas as pd
import numpy as np
import os
import sys
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV, cross_val_score
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.feature_selection import SelectFromModel
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import (
    accuracy_score, classification_report, roc_auc_score, roc_curve, confusion_matrix, balanced_accuracy_score, f1_score
)
import matplotlib.pyplot as plt
import seaborn as sns
import platform

# -------------------- Configuration --------------------
output_dir = "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/ml"
os.makedirs(output_dir, exist_ok=True)
random_seed = 42

# -------------------- Data Loading --------------------
def safe_read_csv(path, **kwargs):
    try:
        df = pd.read_csv(path, **kwargs)
        return df
    except FileNotFoundError:
        print(f"ERROR: File not found: {path}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR reading {path}: {e}")
        sys.exit(1)

asv_table_path = "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/dada2/asv_table.tsv"
asv_table_df = safe_read_csv(asv_table_path, sep='\t')
if 'SampleID' in asv_table_df.columns:
    asv_table_df = asv_table_df.set_index('SampleID')
else:
    raise ValueError("Expected 'SampleID' column in ASV table.")

if asv_table_df.index.duplicated().any():
    print("WARNING: Duplicate SampleIDs found in ASV table.")
if asv_table_df.isnull().any().any():
    print("WARNING: Missing values found in ASV table.")

asv_taxonomy_path = "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/dada2/asv_taxonomy.tsv"
tmp_tax = safe_read_csv(asv_taxonomy_path, sep='\t')
tax_index_col = None
for col in ['ASV', 'ASV_ID', 'Feature']:
    if col in tmp_tax.columns:
        tax_index_col = col
        break
if tax_index_col:
    asv_taxonomy_df = safe_read_csv(asv_taxonomy_path, sep='\t', index_col=tax_index_col)
else:
    asv_taxonomy_df = tmp_tax

metadata_path = "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/metadata.tsv"
metadata_df = safe_read_csv(metadata_path, sep='\t')
if 'Sample' in metadata_df.columns:
    metadata_df = metadata_df.set_index('Sample')
else:
    raise ValueError("Expected 'Sample' column in metadata.")

if metadata_df.index.duplicated().any():
    print("WARNING: Duplicate Sample names found in metadata.")
if metadata_df.isnull().any().any():
    print("WARNING: Missing values found in metadata.")

# -------------------- Data Preparation --------------------
zero_sum_samples = asv_table_df.sum(axis=1) == 0
if zero_sum_samples.any():
    print(f"WARNING: {zero_sum_samples.sum()} samples with all zero counts will be dropped.")
    asv_table_df = asv_table_df.loc[~zero_sum_samples]

asv_table_normalized = asv_table_df.div(asv_table_df.sum(axis=1), axis=0)

common_samples = list(set(asv_table_normalized.index) & set(metadata_df.index))
dropped_samples = set(asv_table_normalized.index).symmetric_difference(metadata_df.index)
print(f"Retained {len(common_samples)} samples. Dropped {len(dropped_samples)} samples.")
if len(common_samples) == 0:
    raise ValueError("No common samples between ASV table and metadata.")
asv_table_aligned = asv_table_normalized.loc[common_samples]
metadata_aligned = metadata_df.loc[common_samples]

with open(os.path.join(output_dir, "dropped_samples.txt"), "w") as f:
    for s in dropped_samples:
        f.write(f"{s}\n")

if 'Group' not in metadata_aligned.columns:
    raise ValueError("Expected 'Group' column in metadata.")
label_encoder = LabelEncoder()
encoded_labels = label_encoder.fit_transform(metadata_aligned['Group'])
label_map = dict(zip(label_encoder.classes_, range(len(label_encoder.classes_))))
print("Label mapping:", label_map)

# -------------------- Data Splitting --------------------
X = asv_table_aligned
y = encoded_labels
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=random_seed, stratify=y
)

cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=random_seed)

# -------------------- Define Models and Pipelines --------------------
models = {
    "RandomForest": RandomForestClassifier(random_state=random_seed, class_weight='balanced'),
    "GradientBoosting": GradientBoostingClassifier(random_state=random_seed),
    "LogisticRegression": LogisticRegression(
        random_state=random_seed, solver='liblinear', class_weight='balanced'
    ),
    "SVM": SVC(random_state=random_seed, probability=True, class_weight='balanced'),
}

feature_selectors = {
    "RandomForest": SelectFromModel(RandomForestClassifier(n_estimators=100, random_state=random_seed, class_weight='balanced')),
    "GradientBoosting": SelectFromModel(RandomForestClassifier(n_estimators=100, random_state=random_seed, class_weight='balanced')),
    "LogisticRegression": SelectFromModel(LogisticRegression(penalty='l1', solver='liblinear', random_state=random_seed, class_weight='balanced')),
    "SVM": SelectFromModel(RandomForestClassifier(n_estimators=100, random_state=random_seed, class_weight='balanced')),
}

pipelines = {}
for name, model in models.items():
    pipelines[name] = Pipeline([
        ('scaler', StandardScaler()),
        ('feature_selection', feature_selectors[name]),
        ('classifier', model)
    ])

# -------------------- Nested Cross-Validation for Hyperparameter Tuning --------------------
# Example for RandomForest (can be extended to other models)
param_grid_rf = {
    'classifier__n_estimators': [100, 200],
    'classifier__max_depth': [None, 5, 10]
}
grid_search_rf = GridSearchCV(
    pipelines['RandomForest'], param_grid_rf, cv=cv, scoring='accuracy', n_jobs=-1
)
grid_search_rf.fit(X_train, y_train)
print("RandomForest best params:", grid_search_rf.best_params_)
pipelines['RandomForest'] = grid_search_rf.best_estimator_

# -------------------- Model Training and Evaluation (No Data Leakage) --------------------
for name, pipeline in pipelines.items():
    print(f"\n=== {name} ===")
    # Cross-validation accuracy (on training set only)
    cv_scores = cross_val_score(
        pipeline, X_train, y_train, cv=cv, scoring='accuracy'
    )
    print(f"Cross-validation accuracy (train only): {cv_scores.mean():.4f} +/- {cv_scores.std():.4f}")

    # Fit on training set only
    pipeline.fit(X_train, y_train)
    y_pred = pipeline.predict(X_test)

    # Classification report
    report = classification_report(y_test, y_pred, target_names=label_encoder.classes_)
    print("Classification Report:\n", report)

    # Balanced accuracy and F1-score
    bal_acc = balanced_accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred, average='weighted')
    print(f"Balanced Accuracy: {bal_acc:.4f}")
    print(f"Weighted F1-score: {f1:.4f}")

    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    cm_df = pd.DataFrame(cm, index=label_encoder.classes_, columns=label_encoder.classes_)

    plt.figure(figsize=(7, 5))
    sns.heatmap(cm_df, annot=True, fmt="d", cmap="Blues")
    plt.title(f"Confusion Matrix ({name})")
    plt.ylabel("Actual")
    plt.xlabel("Predicted")
    cm_path = os.path.join(output_dir, f"{name}_confusion_matrix.png")
    plt.savefig(cm_path)
    plt.close()
    print(f"Confusion matrix saved to {cm_path}")

    # ROC and AUC (macro for multiclass)
    try:
        y_pred_proba = pipeline.predict_proba(X_test)
        if y_pred_proba.shape[1] == 2:
            auc = roc_auc_score(y_test, y_pred_proba[:, 1])
            fpr, tpr, _ = roc_curve(y_test, y_pred_proba[:, 1])
            plt.figure(figsize=(7, 6))
            plt.plot(fpr, tpr, label=f"{name} (AUC = {auc:.2f})")
            plt.plot([0, 1], [0, 1], 'k--', label="Random Chance")
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel("False Positive Rate")
            plt.ylabel("True Positive Rate")
            plt.title(f"ROC Curve ({name})")
            plt.legend(loc="lower right")
            roc_path = os.path.join(output_dir, f"{name}_roc_curve.png")
            plt.savefig(roc_path)
            plt.close()
            print(f"ROC curve saved to {roc_path}")
            print(f"AUC: {auc:.4f}")
            with open(os.path.join(output_dir, f"{name}_auc.txt"), "w") as f:
                f.write(f"AUC: {auc:.4f}\n")
        else:
            # Multiclass ROC/AUC
            auc = roc_auc_score(y_test, y_pred_proba, multi_class='ovr', average='macro')
            print(f"Multiclass macro AUC: {auc:.4f}")
            with open(os.path.join(output_dir, f"{name}_auc.txt"), "w") as f:
                f.write(f"Multiclass macro AUC: {auc:.4f}\n")
    except AttributeError:
        print(f"{name} does not have predict_proba, skipping ROC/AUC.")

    # Save classification report to a file
    report_path = os.path.join(output_dir, f"{name}_classification_report.txt")
    with open(report_path, "w") as f:
        f.write(report)
        f.write(f"\nBalanced Accuracy: {bal_acc:.4f}\nWeighted F1-score: {f1:.4f}\n")
    print(f"Classification report saved to {report_path}")

    # Feature importance (from feature selection step)
    feature_selector = pipeline.named_steps['feature_selection']
    selected_features = X_train.columns[feature_selector.get_support()]
    print("\nSelected Features:")
    print(selected_features)
    selected_features_path = os.path.join(output_dir, f"{name}_selected_features.txt")
    with open(selected_features_path, "w") as f:
        for feature in selected_features:
            f.write(f"{feature}\n")
    print(f"Selected features saved to {selected_features_path}")

# -------------------- Save Environment Info --------------------
env_path = os.path.join(output_dir, "python_env_info.txt")
with open(env_path, "w") as f:
    f.write("Python version: " + platform.python_version() + "\n")
    f.write("Platform: " + platform.platform() + "\n")
    f.write("Random seed: 42\n")
    f.write("Packages:\n")
    for pkg in ['sklearn', 'pandas', 'seaborn', 'matplotlib', 'numpy']:
        try:
            mod = __import__(pkg)
            f.write(f"{pkg}: {mod.__version__}\n")
        except Exception:
            f.write(f"{pkg}: not found\n")
print(f"\nEnvironment info saved to {env_path}")
print("\nAll model results saved in:", output_dir)
