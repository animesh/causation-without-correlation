# EDM and CCM: Practical Implementation Guide

## For Medical Professionals and Data Scientists

This guide provides step-by-step explanations of the algorithms with practical code examples suitable for understanding and implementing dynamical systems analysis.

---

## Part 1: Understanding the Concepts

### The Problem: Why Correlation Fails

Imagine two genes, **A** and **B**, where:
- When A is **low:** B increases as A increases (positive effect)
- When A is **high:** B decreases as A increases (negative effect)

**Result:** Overall Pearson correlation ≈ 0 (the effects cancel)

**But A clearly affects B!** ← This is what CCM detects.

### The Solution: Takens Embedding Theorem

**Key Insight:** In a deterministic system, the past contains information about the future.

```
Original system: Multiple genes (high-dimensional) interact
Time series: We only measure one gene (one-dimensional projection)
Embedding: Use time lags to reconstruct the full system's behavior
```

**Mathematical Idea:**
```
If gene X(t) is influenced by genes Y and Z:
X(t) = f(X, Y, Z)

Then from X alone, we can estimate:
X(t), X(t-τ), X(t-2τ) ≈ (X(t), Y(t), Z(t))

This reconstructed 3D space behaves like the original system!
```

### Why This Works: Nearest Neighbor Prediction

```
Today:    We observe gene expression state
Looking around: Find similar past states (nearest neighbors)
Prediction: Where did those past states go next?
Weighted average of their futures = Our prediction for today
```

If prediction is accurate → dynamics are deterministic ✓
If prediction improves with more data → true causality ✓

---

## Part 2: Step-by-Step Implementation

### Step 1: Create Shadow Attractor

```python
import numpy as np

def embed_time_series(x, embedding_dim=3, delay=1):
    """
    Convert 1D time series to high-dimensional phase space
    
    Example with x = [1, 2, 3, 4, 5, 6] and dim=3, delay=1:
    Point 1: [1, 2, 3]
    Point 2: [2, 3, 4]
    Point 3: [3, 4, 5]
    Point 4: [4, 5, 6]
    
    Each row = system state at one time point
    Columns = time-delayed coordinates
    """
    n = len(x)
    embedded = np.zeros((n - (embedding_dim - 1) * delay, embedding_dim))
    
    for i in range(embedding_dim):
        start = i * delay
        end = n - (embedding_dim - 1 - i) * delay
        embedded[:, i] = x[start:end]
    
    return embedded

# Example
time_series = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 0.8, 0.6, 0.4, 0.2])
attractor = embed_time_series(time_series, embedding_dim=3, delay=1)
print("Attractor shape:", attractor.shape)
print("First few points:\n", attractor[:3])
```

**Output:**
```
Attractor shape: (7, 3)
First few points:
 [[0.1 0.3 0.5]
  [0.3 0.5 0.7]
  [0.5 0.7 0.9]]
```

**Interpretation:**
Each row is a "snapshot" of the system in 3D space. Similar rows should evolve similarly.

### Step 2: Find Nearest Neighbors

```python
from scipy.spatial.distance import cdist

def find_nearest_neighbors(attractor, point_index, n_neighbors=4):
    """
    Find the nearest neighbors of a given point on the attractor
    """
    point = attractor[point_index:point_index+1, :]
    
    # Calculate distances to all other points
    distances = cdist(point, attractor, metric='euclidean').flatten()
    
    # Sort by distance and get indices
    sorted_indices = np.argsort(distances)
    
    # Return the n nearest (excluding the point itself)
    neighbors = sorted_indices[1:n_neighbors+1]
    neighbor_distances = distances[neighbors]
    
    return neighbors, neighbor_distances

# Example
neighbors, distances = find_nearest_neighbors(attractor, point_index=3, n_neighbors=3)
print("Neighbors of point 3:", neighbors)
print("Distances:", distances)
```

### Step 3: Weighted Average Prediction

```python
def predict_next_value(attractor, time_series, point_index, n_neighbors=4, tp=1):
    """
    Predict where the system goes next using nearest neighbor analogy
    
    tp = time ahead to predict
    """
    neighbors, distances = find_nearest_neighbors(attractor, point_index, n_neighbors)
    
    # Weight inversely by distance (closer points get more weight)
    weights = np.exp(-distances / np.mean(distances))
    weights = weights / np.sum(weights)
    
    # Where do the neighbors go at time tp?
    target_times = neighbors + tp
    target_values = time_series[target_times[target_times < len(time_series)]]
    
    # Weighted average
    prediction = np.average(target_values, weights=weights[:len(target_values)])
    
    return prediction

# Example: Predict value at t=5 based on neighbors of t=3
prediction = predict_next_value(attractor, time_series, point_index=3, tp=1)
actual = time_series[4]  # Actual value one step ahead
print(f"Predicted: {prediction:.3f}, Actual: {actual:.3f}")
```

### Step 4: Calculate Prediction Skill

```python
from scipy.stats import pearsonr

def simplex_projection(attractor, time_series, tp=1):
    """
    Test predictability by making predictions for all points
    and calculating correlation with actual values
    """
    n_points = attractor.shape[0]
    predictions = []
    actuals = []
    
    for i in range(n_points - tp):
        try:
            pred = predict_next_value(attractor, time_series, i, n_neighbors=4, tp=tp)
            actual = time_series[i + tp]
            predictions.append(pred)
            actuals.append(actual)
        except:
            continue
    
    predictions = np.array(predictions)
    actuals = np.array(actuals)
    
    # Correlation = prediction skill
    skill, pvalue = pearsonr(predictions, actuals)
    
    return {
        'predictions': predictions,
        'actuals': actuals,
        'skill': skill,
        'pvalue': pvalue,
        'is_predictable': skill > 0.1
    }

# Example
result = simplex_projection(attractor, time_series)
print(f"Prediction skill: {result['skill']:.3f}")
print(f"P-value: {result['pvalue']:.4f}")
print(f"Predictable: {result['is_predictable']}")
```

### Step 5: Test for Nonlinearity (S-map)

```python
def smap_local_regression(attractor, time_series, theta=0, tp=1):
    """
    S-map = Sequentially Weighted Global Linear Map
    
    theta = nonlinearity parameter
      theta=0 → global linear regression
      theta>0 → increasingly local (nonlinear)
    """
    n_points = attractor.shape[0]
    predictions = []
    actuals = []
    
    for i in range(n_points - tp):
        point = attractor[i:i+1, :]
        
        # Distance-based weights
        distances = cdist(point, attractor, metric='euclidean').flatten()
        avg_dist = np.mean(distances)
        
        # S-map weights: exponential with theta
        if theta == 0:
            weights = np.ones_like(distances)
        else:
            weights = np.exp(-theta * distances / avg_dist)
        
        weights = weights / np.sum(weights)
        
        # Weighted prediction
        target_idx = i + tp
        if target_idx < len(time_series):
            # Simple weighted average (could be more sophisticated)
            prediction = np.average(time_series[max(0, i-4):i+tp+4], 
                                   weights=weights[max(0, i-4):i+tp+4])
            predictions.append(prediction)
            actuals.append(time_series[target_idx])
    
    predictions = np.array(predictions)
    actuals = np.array(actuals)
    
    if len(predictions) > 2:
        skill, _ = pearsonr(predictions, actuals)
    else:
        skill = np.nan
    
    return skill

# Test different theta values
thetas = np.linspace(0, 2, 11)
skills = [smap_local_regression(attractor, time_series, theta=t) for t in thetas]

print("Theta values:", thetas)
print("Skills:      ", [f"{s:.3f}" for s in skills])

# Nonlinearity = max skill > linear skill
linear_skill = skills[0]  # theta=0
nonlinear_skill = max(skills)
is_nonlinear = nonlinear_skill > (linear_skill + 0.05)
print(f"\nLinear skill: {linear_skill:.3f}")
print(f"Nonlinear skill: {nonlinear_skill:.3f}")
print(f"Is nonlinear: {is_nonlinear}")
```

### Step 6: Convergent Cross Mapping (CCM)

```python
def convergent_cross_mapping(cause, effect, embedding_dim=3, library_sizes=None):
    """
    Test if effect series contains information about cause series
    
    If cause → effect (causal), then:
    - Effect's attractor encodes cause's behavior
    - We can predict cause from effect's state space
    - Skill improves as we have more data (convergence)
    """
    
    # Reconstruct effect's attractor
    attractor_effect = embed_time_series(effect, embedding_dim)
    
    if library_sizes is None:
        library_sizes = np.round(np.linspace(10, len(cause)//2, 8)).astype(int)
    
    ccm_skills = []
    
    for lib_size in library_sizes:
        # Sample a library of points
        lib_indices = np.random.choice(len(cause), min(lib_size, len(cause)), 
                                      replace=False)
        lib_cause = cause[lib_indices]
        lib_attractor = attractor_effect[lib_indices, :]
        
        # Try to predict cause from effect's state space
        predictions = []
        actuals = []
        
        for i in range(attractor_effect.shape[0]):
            point = attractor_effect[i:i+1, :]
            distances = cdist(point, lib_attractor, metric='euclidean').flatten()
            
            # Find nearest neighbor
            nn_idx = np.argmin(distances)
            predictions.append(lib_cause[nn_idx])
            if i < len(cause):
                actuals.append(cause[i])
        
        predictions = np.array(predictions)
        actuals = np.array(actuals[:len(predictions)])
        
        # Calculate skill
        if len(predictions) > 2:
            skill, _ = pearsonr(predictions, actuals)
            ccm_skills.append(skill)
        else:
            ccm_skills.append(np.nan)
    
    return {
        'library_sizes': library_sizes,
        'ccm_skills': np.array(ccm_skills),
        'convergence': np.nanmean(np.diff(ccm_skills)) > 0
    }

# Example: Two time series where A influences B
time_series_a = np.sin(np.linspace(0, 4*np.pi, 50)) + 0.1*np.random.randn(50)
time_series_b = np.sin(np.linspace(0, 4*np.pi, 50) + 0.5) + 0.1*np.random.randn(50)

result = convergent_cross_mapping(time_series_a, time_series_b)
print("CCM Analysis: Can we predict A from B's state space?")
print("Library sizes:", result['library_sizes'])
print("CCM skills:  ", [f"{s:.3f}" for s in result['ccm_skills']])
print("Convergence: ", result['convergence'])
```

---

## Part 3: Full Workflow Example

```python
# 1. GENERATE DATA
np.random.seed(42)
t = np.linspace(0, 10, 100)
gene_x = np.sin(t) + 0.1*np.random.randn(100)
gene_y = np.sin(t + 0.3) * (1 - 0.5*np.sin(t)) + 0.1*np.random.randn(100)

# Linear correlation might miss this relationship
corr = np.corrcoef(gene_x, gene_y)[0, 1]
print(f"Linear correlation: {corr:.3f}")

# 2. EMBEDDING
E = 3
attractor_x = embed_time_series(gene_x, E)
attractor_y = embed_time_series(gene_y, E)

# 3. PREDICTABILITY TEST
pred_x = simplex_projection(attractor_x, gene_x)
pred_y = simplex_projection(attractor_y, gene_y)
print(f"Gene X skill: {pred_x['skill']:.3f}, Predictable: {pred_x['is_predictable']}")
print(f"Gene Y skill: {pred_y['skill']:.3f}, Predictable: {pred_y['is_predictable']}")

# 4. NONLINEARITY TEST
smap_result_x = [smap_local_regression(attractor_x, gene_x, theta) 
                 for theta in np.linspace(0, 2, 9)]
smap_result_y = [smap_local_regression(attractor_y, gene_y, theta) 
                 for theta in np.linspace(0, 2, 9)]

print(f"Gene X nonlinearity: {max(smap_result_x) - smap_result_x[0]:.3f}")
print(f"Gene Y nonlinearity: {max(smap_result_y) - smap_result_y[0]:.3f}")

# 5. CAUSAL ANALYSIS (Does X cause Y?)
ccm_result = convergent_cross_mapping(gene_x, gene_y, embedding_dim=3)
print(f"CCM (X→Y): {ccm_result['ccm_skills']}")
print(f"Convergence: {ccm_result['convergence']}")
```

---

## Part 4: Biological Interpretation

### Example: Checkpoint Control

```
Scenario: G1/S Checkpoint with signal integration

Input genes (uncorrelated):
  A: Nutrient sensor (Pearson r = -0.02 with checkpoint)
  B: Stress signal (Pearson r = 0.03)
  C: DNA damage (Pearson r = -0.05)

Output: Checkpoint protein (WHI5/RB)

Problem: Correlation-based approach says nothing related
Solution: CCM detects that A, B, C all influence checkpoint

Why this works:
  - When A=high AND B=low: Response (nonlinear)
  - When A=low AND B=high: No response
  - Different combinations give different outputs
  - Overall correlation can be near zero due to averaging
  - But causality is real and testable by CCM
```

### Validation Strategy

```
1. Predict causal network structure using CCM
   └─ Identify WHI5 as hub with 8 inputs

2. Perturb the system (overexpression or knockout)
   └─ WHI5 overexpression

3. Measure downstream effects (co-prediction)
   └─ Compare wildtype vs. overexpression attractors
   └─ If CCM was correct, dynamics should change significantly

4. Statistical test
   └─ 71% of predicted targets show significant dynamic change
   └─ Only 52% of non-predicted genes show change (random)
   └─ p < 10^-8 (highly significant)

Result: CCM predictions validated!
```

---

## Part 5: Practical Considerations

### Choosing Embedding Dimension

```
Too small (E=1): Linear approximation, bad
E=2: Might miss complexity
E=3-5: Usually optimal for gene networks
E>5: Curse of dimensionality, need more data

Rule of thumb:
  Required data points ≥ 2^E
  For E=3: need ≥ 8 points
  For E=5: need ≥ 32 points
  For E=8: need ≥ 256 points ← Paper uses this
```

### Time Series Length

```
Minimum: (E + 1) × 10
Recommended: (E + 1) × 50+

For gene expression:
  Cell cycle: 50-60 points per cycle
  Long-term: 100+ points
  Short pulses: 20-30 points
```

### Statistical Significance

```
Simplex (predictability):
  p < 0.05 threshold
  skill > 0.1 indicates determinism

S-map (nonlinearity):
  Δskill > 0.05 (5% improvement)
  Or: skill(θ*) > skill(0) with p < 0.05

CCM (causality):
  Convergence: increasing skill with library size
  Significance: Wilcoxon test (p < 0.05)
  Validation: Experimental manipulation
```

---

## Part 6: Troubleshooting

### Problem: Low Predictability

**Causes:**
- Series too short
- Embedding dimension wrong
- Stochastic noise dominates
- Non-stationary dynamics

**Solutions:**
- Increase time series length
- Optimize embedding dimension
- Test different E values
- Check for temporal trends

### Problem: Spurious Causality

**Causes:**
- Both series driven by common third variable
- Autocorrelation artifact
- Too small library size

**Solutions:**
- Test both directions (X→Y and Y→X)
- True causality: asymmetric CCM
- Partialing out common drivers
- Use larger library sizes

### Problem: Nonlinearity Not Detected

**Causes:**
- Dynamics are actually linear
- Time series too short
- Noise obscures nonlinearity
- Threshold too stringent

**Solutions:**
- Reduce threshold (Δskill > 0.01)
- Increase time series length
- Smooth data lightly
- Check underlying biology

---

## References

1. **Takens Theorem:** Takens, F. (1981)
2. **Convergent Cross Mapping:** Sugihara et al. (2012)
3. **S-map:** Sugihara (1994)
4. **Yeast Application:** Pao et al. (2026)

---

## Summary

**EDM and CCM provide:**
- Detection of causation in nonlinear systems
- Works when correlation fails
- Rooted in dynamical systems theory
- Experimentally validatable

**Key advantages:**
- No a priori model needed
- Detects hidden causal links
- Asymmetric (directional)
- Convergence tests causality (not just correlation)

**Key limitations:**
- Requires time series data
- Computational cost O(n²) or O(n³)
- Parameters need optimization
- Needs experimental validation
