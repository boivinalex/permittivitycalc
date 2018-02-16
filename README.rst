permittivity-calc
=================
Scripts to calculate and plot the complex permittivity from S-parameter data acquired with transmission-line measurements

Overview
--------
permittivity-calc is a Python package made to take S-parameter data output from METAS VNA Tools II (https://www.metas.ch/metas/en/home/fabe/hochfrequenz/vna-tools.html) and process it to calculate and plot the complex permittivity of a material measured in a coaxial transmission line.

Currently, permittivity-calc uses the New Non-iterative Method for permittivity calculation from S-parameters from [Boughriet,1997]_ which assumes that the material is non-magnetic (i.e. \mu = 1).

References
----------
.. [Boughriet,1997] Boughriet, A. H., Legrand, C., & Chapoton, A. (1997). Noniterative stable transmission/reflection method for low-loss material complex permittivity determination. IEEE Transactions on Microwave Theory and Techniques, 45(1), 52–57. http://doi.org/10.1109/22.552032
