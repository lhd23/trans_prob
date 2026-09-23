<p align="center">
  <img src="docs/transitions_N5000_R200_L1000_transition.gif"
       alt="Transition probability animation"
       width="700">
</p>

Code to compute transition probabilities from a
1D cosmological N-body simulation with 1LPT initial conditions.

## Run

```sh
python run.py transitions
python run.py peak
```

* **transitions:** runs a set of simulations to estimate the transition histograms and
  the evolution of the first three separation cumulants.
* **peak:** runs another set of simulations to measure the evolution of the correlation
  function at five redshifts. Does not save transition pairs or full snapshots.

Both use the supplied Planck 2018 power spectrum table, growing-mode momenta,
the exact periodic rank-ordered force, and the same scale-factor-integrated
kick-drift-kick integrator. Positions are comoving; each equal-mass particle
represents a planar slab. No force softening in this one-dimensional simulation.

Change the corresponding file in `configs/` first
if a smaller test is wanted. Default are production runs that take tens of minutes on a laptop.

To remake figures without simulating:

```sh
python run.py transitions --plot-only
python run.py transitions --plot-only --linear-y
python run.py peak --plot-only
```

To change histogram binning without simulating, reanalyse the saved pairs:

```sh
python src/analyse.py data/transitions_N5000_R200_L1000 --histogram-bins 135
python plots/plot_transitions.py data/transitions_N5000_R200_L1000
```

## What each file does

```
configs/transitions.json       transition simulation settings
configs/peak.json              peak simulation settings
src/settings.py               paths, basic settings checks, numerical metadata
src/initial_conditions.py     spectrum, growth, Gaussian displacements and momenta
src/dynamics.py               periodic force and kick-drift-kick evolution
src/transitions.py            select displaced initial pairs and follow labels
src/simulate_transitions.py   evolve and save selected-pair separations in batches
src/analyse.py                normalized probabilities, cumulants, sampled paths
src/simulate_peak.py          evolve and accumulate correlation means and errors
plots/plot_transitions.py     export transition/cumulant tables and invoke Gnuplot
plots/histograms.gp           four rows, with trajectories and four histograms
plots/cumulants.gp            two-panel cumulant evolution
plots/plot_peak.py            spectrum-based models, peak tables, and Gnuplot call
plots/peak.gp                 five-panel peak evolution
run.py                       short command sequence for either workflow
check.py                     small numerical checks
data/Pk_Planck18_large.dat    supplied redshift-zero input spectrum
```

## Measurements and stored data

The transition initial cells have q/L = 0.10, 0.05, 0.02, 0.01 and width 0.001L.
Selection uses the actual initial Eulerian positions **after displacement**.
Every oriented non-self pair in `[q-width/2, q+width/2)` is followed by its
original particle labels through all crossings. Later separations use the
signed periodic interval `[-L/2, L/2)`.

There are 20 logarithmic snapshots between a=0.02 and a=10, plus a=0.1, 0.3,
0.6, 1.0: 24 in total. Only those four extra times enter the histograms. The
cumulant figure omits their markers and shows magnitudes of the first three
population cumulants of `(r-q)/q` on logarithmic axes, for the first and last
q rows. The signed cumulants remain available in `analysis.npz` and the tables.
There are no cumulant model curves.