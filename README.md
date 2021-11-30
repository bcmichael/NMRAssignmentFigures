This code draws figures showing which atomic sites within a peptide sequence have had their chemical shifts assigned via
NMR experiments. The drawings are created using Luxor.jl and data is read from resonance lists saved from the sparky NMR
assignment program to identify assigned sites. The example scripts show how to produce a figure showing proton, carbon,
and nitrogen assignments:

![HCN_figure](https://github.com/bcmichael/NMRAssignmentFigures/blob/master/example_figure.svg)

or with only carbon and nitrogen, but not protons:

![CN_figure](https://github.com/bcmichael/NMRAssignmentFigures/blob/master/example_figure_no_protons.svg)
