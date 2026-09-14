# Data policy

No study data are distributed in this repository.

The phenotype schema accepted by the analysis is strictly:

```text
ENV,GEN,REP,GY
```

Only grain yield (`GY`) is within scope. Keep real phenotype and marker files outside version control; the repository's `.gitignore` blocks the expected input filenames and private-data directories.
