# `navier-2D`

![](.github/car-sim.gif)

Navier-2D is a 2D CFD simulator. Read the [article about it](https://kyletennison.com/articles/navier-2d) on my website.

## Install

```sh
cargo build --release && alias navier-2d=target/release/navier-2d
```

## Usage


```sh
Usage: navier-2d [OPTIONS] [MASK_PATH]

Arguments:
  [MASK_PATH]  The path to a PNG image to use as the solid object.

Options:
      --input-json <INPUT_JSON>
          An input file with pre-loaded parameters.
      --mode <MODE>
          The mode to run the simulation in: `stream` or `video` [default: video]
      --input-json-savepath <INPUT_JSON_SAVEPATH>
          Optional directiory to save input file to.
      --frames-dir <FRAMES_DIR>
          An optional directory pointing to where frames should be saved.
  -r, --retain-frames
          Whether or not frames should be retained after saving. Always false if streaming.
  -d, --display-video
          Whether the frame animation should play after solving. Disabled if streaming.
      --inflow-x <INFLOW_X>
          Inflow x velocity. [default: 8.0]
      --inflow-y <INFLOW_Y>
          Inflow y velocity. [default: 0.0]
      --length-x <LENGTH_X>
          Domain length in x axis. [default: 5.0]
      --length-y <LENGTH_Y>
          Domain length in y axis. [default: 5.0]
  -s, --simtime <SIMTIME>
          Simulation time in seconds. [default: 10]
      --density <DENSITY>
          Fluid density in g/ml [default: 1.0]
      --viscosity <VISCOSITY>
          Shear viscosity [default: 0.002]
      --cfl <CFL>
          Maximum CFL [default: 1.0]
  -h, --help
          Print help
  -V, --version
          Print version
```

You can display this window with `navier-2d --help`.