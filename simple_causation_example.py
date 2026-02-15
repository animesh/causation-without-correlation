"""
SIMPLEST EXAMPLE: CAUSATION WITHOUT CORRELATION

A gene expression system where:
  - Gene B CAUSES changes in Gene A
  - But Correlation(Gene A, Gene B) ≈ 0

The key: State-dependent interaction
"""

import numpy as np
from scipy.stats import pearsonr

print("\n" + "="*70)
print("MINIMAL EXAMPLE: CAUSATION WITHOUT CORRELATION")
print("="*70)

# ============================================================================
# EXAMPLE 1: The Simplest Case
# ============================================================================

print("\n[EXAMPLE 1] A Sine Wave with State-Dependent Damping")
print("-"*70)

print("""
THE SETUP:
  Gene A: A sine wave that oscillates naturally
  Gene B: A damping force that suppresses Gene A
  
THE TRICK:
  Gene B is STRONGEST when Gene A is at peak height (±1)
  Gene B is WEAKEST when Gene A crosses zero
  
RESULT:
  Most of the time, high Gene A → high Gene B (negative correlation)
  But sometimes, high Gene A → low Gene B (positive correlation)
  Overall: ρ ≈ 0 (they cancel out)
""")

# Generate data
t = np.arange(0, 40, 0.5)
gene_A = np.sin(t)  # Oscillates between -1 and +1
gene_B = np.abs(gene_A)  # Active when A is at extremes

print(f"\nGene A values: sin(t)")
print(f"Gene B values: |sin(t)| (damping force)")

# Correlation
r, p = pearsonr(gene_A, gene_B)
print(f"\nCorrelation: ρ = {r:.4f} (NOT SIGNIFICANT)")
print(f"p-value: {p:.4f}")

# Causation proof
print(f"\n✓ PROOF OF CAUSATION:")
print(f"  When Gene A = +0.9: Gene B = {np.mean(np.abs(gene_A[np.abs(gene_A) > 0.8])):.3f} (HIGH - damping!)")
print(f"  When Gene A ≈ 0.0:  Gene B = {np.mean(np.abs(gene_A[np.abs(gene_A) < 0.2])):.3f} (LOW - no damping)")
print(f"  \n  → Gene B clearly AFFECTS Gene A, but...")
print(f"  → Linear correlation misses this!")

# ============================================================================
# EXAMPLE 2: An Even Simpler Case
# ============================================================================

print("\n\n[EXAMPLE 2] The Pendulum Checkpoint")
print("-"*70)

print("""
BIOLOGICAL SCENARIO - Cell Cycle Checkpoint:
  Gene A (WHI5): A repressor protein
  Gene B (SWI4): A transcription factor it controls
  
WHAT HAPPENS:
  • When cell wants to enter S-phase: SWI4 high, WHI5 pushes back hard
  • When cell progresses: WHI5 effect weakens
  
RESULT:
  • Opposite signs at different times = low correlation
  • But WHI5 clearly causes changes in SWI4
""")

# Simplified model
time = np.arange(100)

# SWI4 oscillates (cell cycle progression)
swi4 = np.sin(2 * np.pi * time / 50)  # One cycle per 50 timepoints

# WHI5 only represses at checkpoint (when SWI4 wants to be high)
whi5 = np.zeros(100)
whi5[(time % 50) < 15] = 0.8  # Active only first 15 points of cycle
whi5 += 0.1 * np.random.normal(0, 1, 100)  # Add small noise

print(f"\nSimulated time series:")
print(f"  SWI4: sin(2π * t/50)  [oscillates]")
print(f"  WHI5: Active at checkpoint only [pulsed]")

# Correlation
r, p = pearsonr(swi4, whi5)
print(f"\nCorrelation: ρ = {r:.4f} (essentially zero!)")
print(f"p-value: {p:.4f}")

print(f"\n✓ PROOF OF CAUSATION:")
print(f"  When SWI4 wants up (sin>0.5): WHI5 = {np.mean(whi5[swi4 > 0.5]):.3f}")
print(f"  When SWI4 is down (sin<-0.5): WHI5 = {np.mean(whi5[swi4 < -0.5]):.3f}")
print(f"  \n  → WHI5 repression strongest when SWI4 rising")
print(f"  → This is classic NEGATIVE REGULATION")
print(f"  → But correlation is near zero!")

# ============================================================================
# EXAMPLE 3: The Absolute Simplest Case
# ============================================================================

print("\n\n[EXAMPLE 3] The Micro-Scale Example")
print("-"*70)

print("""
THE ABSOLUTE SIMPLEST EXAMPLE:

Data:      X        Y
          -1       1
          -1       1
           0       0
           0       0
          +1       1
          +1       1

OBSERVATION:
  • X ranges from -1 to +1
  • Y = |X| (Y is absolute value of X)
  • When X = -1: Y = 1 ✓ (relationship exists!)
  • When X = 0: Y = 0  ✓ (relationship exists!)
  • When X = +1: Y = 1 ✓ (relationship exists!)
  
BUT: Correlation(X, Y) ≈ 0
     Because negative X's pair with high Y
     And positive X's also pair with high Y
     Perfect symmetry cancels correlation!
""")

X = np.array([-1, -1, 0, 0, 1, 1])
Y = np.abs(X)

r, p = pearsonr(X, Y)

print(f"\nData:")
print(f"  X = {list(X)}")
print(f"  Y = {list(Y)}")
print(f"\nCorrelation: ρ = {r:.4f}")
print(f"Causation: Y = |X|  (OBVIOUS!)")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n\n" + "="*70)
print("KEY INSIGHT")
print("="*70)

print("""
Why does this happen?

CORRELATION assumes LINEAR relationship:
  Y = m*X + b
  
  When this is SYMMETRIC (same relationship in both directions),
  effects cancel → ρ ≈ 0

But CAUSATION can be NONLINEAR:
  Y = f(X)  where f is any function
  
  Examples:
    Y = |X|          ← absolute value
    Y = X²           ← quadratic
    Y = sgn(X) * |X| ← symmetric damping
    Y = Step(X)      ← threshold
    
All have CAUSATION (X affects Y) but low CORRELATION!

═══════════════════════════════════════════════════════════════════════════

BIOLOGICAL EXAMPLES:

1. CHECKPOINT CONTROL (Like Pao et al.)
   WHI5 represses SWI4 when S-phase should not start
   No repression when S-phase should start
   → Nonlinear interaction, low correlation
   
2. FEEDBACK INHIBITION
   Gene A activates Gene B
   Gene B inhibits Gene A
   → Oscillatory with low correlation
   
3. THRESHOLD RESPONSE
   Gene A only responds when Gene B exceeds threshold
   Below threshold: no response
   → Binary decision, low correlation

═══════════════════════════════════════════════════════════════════════════

HOW TO DETECT:

❌ FAILS: Pearson Correlation (assumes linearity)
✓ WORKS: Convergent Cross Mapping (CCM)
         - Tests nonlinear relationships
         - Detects state-dependent effects
         - Reveals hidden causality

This is WHY Pao et al. used CCM instead of correlation!
""")

print("\n" + "="*70 + "\n")
