# CentroSymmetryParameter

`CentroSymmetryParameter` computes centrosymmetry values for each atom.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- CentroSymmetryParameter
```

## CLI

Usage:

```bash
centrosymmetry <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--numNeighbors <int>` | No | Even number of neighbors, up to `32`. | `12` |
| `--mode <conventional\|matching>` | No | CSP evaluation mode. | `conventional` |
| `--threads <int>` | No | Maximum worker threads. | auto |
| `--help` | No | Print CLI help. | |
