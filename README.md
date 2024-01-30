# drprobe-config

A user-friendly front-end for the command line tools associated with the DrProbe software package for multislice simulations, written by Dr. Juri Barthel. `drprobe-config` can be used either as a standalone CLI command or as a Julia module. In either case, the parameters for the calculation are supplied through a configuration file in YAML format. An example configuration file is available in the `example` folder. A major benefit of `drprobe-config` is that it supports parallel computation which greatly speeds up any calculation when many cores are available (e.g. on a computing cluster). 

## Installation
Clone the GitHub repo:
```
git clone https://github.com/evanderveer/drprobe-config
```

`drprobe-config` requires that the DrProbe command line tools be available on the system. If necessary they can be installed by running:

```
> julia ./drprobe-config/deps/install_drprobe_cli.jl
```
which will place them in `~/.local/bin`.

To make the CLI command available, run:

```
julia --project=./drprobe-config ./drprobe-config/deps/build.jl install
```


## Usage
`drprobe-config` can be used in two ways: as a CLI command or as a Julia module.

### CLI command
```
> drprobe-config ./config.yml
```

### Julia module
In file `run-drprobe.jl`:
```
using DrProbeConfig
run_drprobe("config.yml")
```
or
```
using DrProbeConfig
config = load_config("config.yml")
run_drprobe(config)
```
then run
```
> julia run-drprobe.jl
```

### Multithreading
To run drprobe-config multithreaded, make sure to set the JULIA_NUM_THREADS environment variable before starting Julia. In this case, the calculation is split up by scan line, so increasing the number of threads beyond the number of scan lines is useless.

```
> export JULIA_NUM_THREADS=10
> julia run-drprobe.jl

> export JULIA_NUM_THREADS=10
> drprobe-config ./config.yml
```

## Authors and Copyright
`drprobe-config` was written by:

Ewout van der Veer \
University of Groningen, Groningen, The Netherlands

DrProbe command-line tools (MSA, CELSLC) were written by:

Juri Barthel, \
Forschungszentrum Jülich GmbH, 52425 Jülich, Germany

Copyright (c) 2008 - 2023 - Forschungszentrum Jülich GmbH\
Published under the GNU General Public License, version 3.

