# Lid driven cavity using Vorticity-Stream function formulation

## Stream function solver

Set the grid size and test the stream function solver

```bash
julia test_stream.jl
```

We obtain a plot like this

<img align="center" src="output/test_stream.svg">

## Lid-driven cavity

To run the lid-driver cavity code

```bash
julia main.jl
```

The solution at the end is shown below.

<img align="center" src="output/stream_vort.svg">

<img align="center" src="output/velocity1.svg">

<img align="center" src="output/velocity2.svg">
