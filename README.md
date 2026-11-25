# Biharmonic source reconstruction

MATLAB code for the smooth and piecewise-constant source examples in *Quantitative Direct Sampling Methods for Source Reconstruction in the Biharmonic Wave Equation* by Xiaodong Liu, Qingxiang Shi, and Jing Wang. It includes the reconstructions in Figures 2 and 5, the smooth-source errors in Table 1, and the Lanczos-window comparison for the bear-shaped source.

## Running the examples

The scripts were run with MATLAB R2025a. They use base MATLAB functions and generate the synthetic data locally.

1. Download or clone the repository. If you downloaded a ZIP file, unzip it.
2. Open MATLAB and set **Current Folder** to the directory containing this README.
3. Run the following in the Command Window.

    ~~~matlab
    run_all
    ~~~

    You can also run `run_smooth` or `run_bear` individually. To repeat the bear example with the Lanczos window, run `run_bear('lanczos')` after the sharp-cutoff run.

The first run saves the noiseless scattered field and its Laplacian to `data/`. Subsequent runs load those matrices. To regenerate the data after editing a source function, use `for_save_data('smooth',true)` or `for_save_data('bear',true)`.

## Settings

| Example | Sensors | Wavenumber interval | Step |
| --- | ---: | --- | ---: |
| Smooth source, Table 1 | 30 and 60 | [0.5, 30] | 0.5 |
| Smooth source, Table 1 | 30 and 60 | [0.1, 50] | 0.1 |
| Bear, coarse setting | 30 | [0.5, 30] | 0.5 |
| Bear, fine setting | 60 | [0.1, 50] | 0.1 |

The sensors are equally spaced on the circle of radius 3. The forward solve uses composite midpoint quadrature on [-3, 3]², with cell side length 0.01. The reconstruction uses the 401-by-401 vertex grid on [-2, 2]². The integrals in the inverse solve use left rectangles. For the first indicator, the auxiliary wavenumber has the same step as the measurement wavenumber and an upper limit of 30. The radial nodes are 0.01, 0.02, ..., 6.

For the noisy runs, each real and imaginary component is multiplied by `1 + 0.2*xi`. The four sets of uniform draws on [-1, 1] are independent, with a fixed MATLAB random seed. The sharp-cutoff and windowed runs share the same data and noise samples. The Lanczos window is applied to the outer auxiliary wavenumber in the first indicator and to the measurement wavenumber in the second indicator.

## Output

| Directory | Contents |
| --- | --- |
| `results/smooth/` | The six panels for Fig. 2, the reconstruction matrices, and the relative errors for Table 1 |
| `results/bear/sharp/` | The six panels for Fig. 5, the reconstruction matrices, and the error table |
| `results/bear/lanczos/` | The corresponding panels and errors from the Lanczos-window run |

The PNGs are saved in each output directory's `figures/` subdirectory. The error panels show the signed difference between the reconstruction and the source. RE is the discrete relative Euclidean error as a fraction. Near-jump RE uses the same definition on points within 0.05 of the source jump interfaces. The two `reference_*.csv` files give reference values for comparison.

## Code guide

- `parameters.m` sets the common grid sizes, auxiliary integration parameters, and noise seed. The experiment-specific frequency ranges are listed in `run_smooth.m` and `run_bear.m`.
- `general_source_function.m` and `bear_head.m` evaluate the two sources.
- `u_s.m` evaluates the forward field and its Laplacian. `for_save_data.m` saves and reloads the data, and `add_noise.m` applies the independent component noise.
- `Simple_f_yy.m` evaluates both indicators. `left_weights.m` gives the quadrature weights and `spectral_window.m` supplies the optional taper.
- `bear_boundary_masks.m` identifies the reconstruction points near each jump. `export_paper_figures.m` writes the individual panels.
