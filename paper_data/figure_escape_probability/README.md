# Line escape probability figure

This figure uses the final outer iteration (011) of run `64e9e118`. Its input is
`line_escape_selected_64e9e118_iter011.dat`; the run settings are saved in
`params_64e9e118.json`. The four selected transitions are He II lines at
303.8, 243.0, 256.3, and 1640.4 Å. They illustrate the bound–bound escape
calculation; this dataset does not test Fe K fluorescence.

Generate the PDF and PNG in this folder and copy both to `../../Figure/`, where
`../../draft.tex` reads the PDF:

```bash
python3 paper_data/figure_escape_probability/plot_escape_probability.py
```

To generate fresh diagnostic data, rebuild `maindaocl` and run the desired
model. The final outer iteration automatically writes
`results/<run_hash>/line_escape_lines_iterNNN.dat` and
`results/<run_hash>/line_escape_selected_iterNNN.dat`, including with parallel
Cloudy workers. Copy its `line_escape_selected_iterNNN.dat` into this folder,
and update `SELECTED_FILE` in the plotting script. The per-line summary is not
needed to reproduce this figure.

The plotted mean survival fraction uses the local thin line power times each
cell's physical width. The lower panels split the survival fraction into line
escape and electron-scattering contributions; their sum is the black curve.
