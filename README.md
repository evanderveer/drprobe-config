# drprobe-config

A user-friendly front-end for the command line tools associated with the DrProbe software package for multislice simulations, written by Dr. Juri Barthel. `drprobe-config` can be used either as a standalone CLI command or as a Julia module. In either case, the parameters for the calculation are supplied through a configuration file in YAML format. An example configuration file is available in the `test` folder. A major benefit of `drprobe-config` is that it supports parallel computation which greatly speeds up any calculation when many cores are available (e.g. on a computing cluster). 

## Authors and Copyright
`drprobe-config` was written by:

Ewout van der Veer \
University of Groningen, Groningen, The Netherlands

DrProbe command-line tools (MSA, CELSLC) were written by:

Juri Barthel, \
Forschungszentrum Jülich GmbH, 52425 Jülich, Germany

Copyright (c) 2008 - 2023 - Forschungszentrum Jülich GmbH\
Published under the GNU General Public License, version 3.

