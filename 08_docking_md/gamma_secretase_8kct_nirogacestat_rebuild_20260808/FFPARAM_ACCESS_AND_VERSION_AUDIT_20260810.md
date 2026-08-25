# FFParam access and version audit

## Result

The registered academic account was successfully authenticated without retaining credentials in the project. The official authenticated download page exposed `versions/ffparam_v1.2.0.tar.gz` as the latest downloadable archive and did not expose an FFParam-v2 source archive.

The user-supplied archive was compared with a fresh authenticated official download route and retained under its actual software identity:

- filename: `ffparam_v1.2.0.tar.gz`;
- size: 30,611,745 bytes;
- SHA-256: `d9508f3a1590ba9fbfb1d048e832d6726ef6131fd40a70572f5662c5bdc2cbdb`;
- package metadata version: `1.2.0`;
- license: University of Maryland, Baltimore FFParam academic/research license;
- server installation: isolated prefix `/root/autodl-tmp/envs/ssc_ffparam_py310`;
- source storage: user-only restricted directory outside the project tree.

## Scientific boundary

The 2024 FFParam-v2.0 paper describes capabilities that are not present in the downloadable 1.2.0 package, including a CLI-oriented workflow and expanded normal-mode/PED tooling. Therefore:

1. FFParam 1.2.0 may generate and compare only the standard CHARMM additive targets it actually supports.
2. Psi4 1.9.1 and OpenMM 8.1.2 will generate independently retained Hessian, normal-mode/PED, and headless validation artifacts.
3. Every output records the tool and version that produced it.
4. The manuscript and reproducibility package must not say that FFParam-v2 was executed.
5. No ligand parameter is approved merely because FFParam 1.2.0 installed or generated an input file.

## Installation qualification

The source tree compiled under Python 3.10.20. Core modules for topology parsing, Psi4 input/output, OpenMM calculations, molecular coordinates, and CGenFF handling imported successfully. GUI-only dependencies are absent on the headless server by design; this does not convert an untested GUI action into a headless calculation. Formal O6U work proceeds through retained scripts and raw engine outputs.
