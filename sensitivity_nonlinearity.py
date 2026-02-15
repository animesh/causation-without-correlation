import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr

def load_geo_data():
    base_dir = os.path.dirname(__file__)
    filepath = os.path.join(base_dir, 'GSE3431_series_matrix.txt')
    with open(filepath, 'r') as f:
        lines = f.readlines()
    data_start = None
    for i, line in enumerate(lines):
        if '!series_matrix_table_begin' in line:
            data_start = i + 1
            break
    data = []
    probe_ids = []
    n_samples = None
    for line in lines[data_start + 1:]:
        if line.startswith('!'):
            break
        parts = line.strip().split('\t')
        if len(parts) > 1:
            try:
                probe_ids.append(parts[0].strip('"'))
                values = [float(x.strip('"')) for x in parts[1:]]
                if n_samples is None:
                    n_samples = len(values)
                if len(values) == n_samples:
                    data.append(values)
            except ValueError:
                continue
    df = pd.DataFrame(data, columns=[f'T{i+1}' for i in range(len(data[0]))])
    df['ProbeID'] = probe_ids
    return df

def create_shadow_attractor(time_series, embedding_dim=3, time_delay=1):
    ts = np.asarray(time_series).flatten()
    n = len(ts)
    max_idx = n - (embedding_dim - 1) * time_delay
    attractor = np.zeros((max_idx, embedding_dim))
    for i in range(embedding_dim):
        idx_start = i * time_delay
        idx_end = max_idx + i * time_delay
        attractor[:, i] = ts[idx_start:idx_end]
    return attractor

def smap_skill_for_theta(cause, effect, embedding_dim=3, theta=0.0):
    cause = np.asarray(cause).flatten()
    effect = np.asarray(effect).flatten()
    attractor = create_shadow_attractor(effect, embedding_dim)
    if attractor.shape[0] < 5:
        return np.nan
    linear_actuals = []
    linear_predictions = []
    n = attractor.shape[0]
    for i in range(n):
        point = attractor[i:i+1, :]
        distances = cdist(point, attractor, metric='euclidean').flatten()
        avg_dist = np.mean(distances) if np.mean(distances)>0 else 1.0
        if theta == 0.0:
            weights = np.ones_like(distances)
        else:
            weights = np.exp(-theta * distances / avg_dist)
        weights /= weights.sum()
        n_neighbors = min(attractor.shape[1] + 1, len(cause))
        nn_idx = np.argsort(distances)[:n_neighbors]
        pred = np.average(cause[nn_idx], weights=weights[nn_idx])
        linear_predictions.append(pred)
        if i < len(cause):
            linear_actuals.append(cause[i])
    preds = np.array(linear_predictions)
    acts = np.array(linear_actuals[:len(preds)])
    if len(preds) > 2:
        return pearsonr(preds, acts)[0]
    return np.nan

def convergent_cross_mapping(cause, effect, embedding_dim=3, library_sizes=None):
    cause = np.asarray(cause).flatten()
    effect = np.asarray(effect).flatten()
    attractor = create_shadow_attractor(effect, embedding_dim)
    if library_sizes is None:
        max_lib = len(cause) - (embedding_dim - 1) - 1
        library_sizes = np.round(np.linspace(5, max_lib, 8)).astype(int)
    ccm_skills = []
    for lib_size in library_sizes:
        if lib_size > len(cause) or lib_size > attractor.shape[0]:
            ccm_skills.append(np.nan)
            continue
        valid_idx = min(len(cause), attractor.shape[0])
        lib_indices = np.random.choice(valid_idx, min(lib_size, valid_idx), replace=False)
        lib_cause = cause[lib_indices]
        lib_attractor = attractor[lib_indices, :]
        predictions = []
        actuals = []
        for i in range(valid_idx):
            point = attractor[i:i+1, :]
            distances = cdist(point, lib_attractor, metric='euclidean').flatten()
            if len(distances) == 0:
                continue
            nn_idx = np.argmin(distances)
            predictions.append(lib_cause[nn_idx])
            if i < len(cause):
                actuals.append(cause[i])
        predictions = np.array(predictions)
        actuals = np.array(actuals[:len(predictions)])
        if len(predictions) > 2:
            try:
                skill = pearsonr(predictions, actuals)[0]
            except:
                skill = np.nan
            ccm_skills.append(skill)
        else:
            ccm_skills.append(np.nan)
    return library_sizes, np.array(ccm_skills)

def main():
    print('\nRunning sensitivity analysis for nonlinearity and CCM...')
    df = load_geo_data()
    gene_map = {
        'CLN3': ['11369_at'],
        'SWI4': ['5596_at'],
        'CLB2': ['7651_at'],
        'WHI5': ['8464_at'],
        'YHP1': ['6010_at'],
        'MBP1': ['6542_at'],
        'CLN2': ['7993_at'],
        'CLB5': ['7652_at'],
    }
    genes = {}
    for g, probes in gene_map.items():
        for p in probes:
            if p in df['ProbeID'].values:
                ts = df[df['ProbeID']==p].iloc[0, :-1].values.astype(float)
                ts_norm = (ts - np.mean(ts)) / (np.std(ts) + 1e-10)
                genes[g] = ts_norm
                break
    interaction_pairs = [
        ('WHI5', 'SWI4'),
        ('SWI4', 'CLN3'),
        ('CLB2', 'CLN3'),
        ('MBP1', 'SWI4'),
        ('CLN2', 'CLB2'),
    ]
    Es = [2,3,4]
    thetas = [0.0, 0.5, 1.0, 1.5, 2.0]
    libsets = [np.array([5,9,13,17,21,25,29,33]), np.array([10,20,30])]

    for pair in interaction_pairs:
        a,b = pair
        if a not in genes or b not in genes:
            continue
        print(f"\nInteraction: {a} -> {b}")
        cause = genes[a]
        effect = genes[b]
        for E in Es:
            print(f" E={E}")
            lin = smap_skill_for_theta(cause, effect, embedding_dim=E, theta=0.0)
            print(f"  Linear skill (theta=0): {lin:.4f}")
            for th in thetas:
                skill = smap_skill_for_theta(cause, effect, embedding_dim=E, theta=th)
                delta = skill - lin if (not np.isnan(skill) and not np.isnan(lin)) else np.nan
                print(f"   theta={th:.1f}: skill={skill:.4f}, delta={delta:.4f}")
            # CCM sensitivity
            for libset in libsets:
                libs, skills = convergent_cross_mapping(cause, effect, embedding_dim=E, library_sizes=libset)
                print("   CCM libs:", libs.tolist())
                print("   CCM skills:", [f"{s:.4f}" if not np.isnan(s) else 'nan' for s in skills])

if __name__ == '__main__':
    main()
