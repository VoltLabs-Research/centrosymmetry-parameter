# Centrosymmetry Parameter

Computes the per-atom centrosymmetry parameter (CSP) to detect defects and local symmetry.

## Install

```bash
vpm install @voltlabs/centrosymmetry-parameter
```

## CLI

```bash
centrosymmetry <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--num_neighbors <int>` | no | `12` | Even number of neighbors, up to `32`. |
| `--mode <conventional\|matching>` | no | `conventional` | CSP evaluation mode. |
| `--threads <int>` | no | auto | Maximum worker threads. |
| `--help` | no | — | Print CLI help. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_centrosymmetry.parquet` | Centrosymmetry | — |
| `{output_base}_atoms.parquet` | Centrosymmetry Model | AtomisticExporter → glb |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins/centrosymmetry-parameter
