<p align="center">
  <img src="docs/transitions_N5000_R200_L1000_S72_transition.gif"
       alt="Transition probability animation">
</p>

Code to estimate the two-particle transition probabilities from a
1D cosmological N-body simulation with 1LPT initial conditions.
This code was used to obtain the results described in
[arXiv:2609.30251](https://arxiv.org/abs/2609.30251).

## Run

Two options:

```sh
python run.py transitions
python run.py peak
```

* **transitions:** runs a set of simulations to estimate the transition histograms
  and the evolution of the first three separation cumulants.
  Stores 20 snapshots logspaced between a=0.02 and a=10,
  with four further snapshots at a=0.1, 0.3, 0.6, 1.0. 
* **peak:** runs another set of simulations to measure the evolution of the correlation
  function at five redshifts. Does not save transition pairs or full snapshots.

Change the settings file in `configs/` for different separations, bins, times, etc.

For the smoother transition animation, save 72 logarithmically spaced snapshots
between a=0.02 and a=10 in a separate run directory:

```sh
python src/simulate_transitions.py --config configs/transitions_72.json --output-directory data/transitions_N5000_R200_L1000_S72
python src/analyse.py data/transitions_N5000_R200_L1000_S72 --trajectory-max-scale-factor 10
python plots/animate_transitions.py data/transitions_N5000_R200_L1000_S72 docs/transitions_N5000_R200_L1000_S72_transition.gif
```

The animation retains the original dimensions and a four-second loop, with an
average of 18 frames per second. Rendering additionally requires Matplotlib
and Pillow. The 72-snapshot settings omit the four extra histogram times, so
this run is intended for animation rather than the fixed-time histogram figure.

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

## Notes

The transition initial cells have q/L = 0.10, 0.05, 0.02, 0.01 and width 0.001L.
Selection uses the actual initial Eulerian positions **after displacement**.
Every oriented non-self pair in `[q-width/2, q+width/2)` is followed by its
original particle labels through all crossings. Later separations use the
signed periodic interval `[-L/2, L/2)`.
