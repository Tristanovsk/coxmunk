<!-- README wrapper for the docs: rewrites image paths relative to docs/source
     and skips the title line (the page already has one); the installation
     section has its own page (installation.md) -->
```{include} ../../README.md
:relative-images:
:start-line: 1
:end-before: <!-- installation-start -->
```

See {doc}`installation` to install `coxmunk` with conda or pip.

```{include} ../../README.md
:relative-images:
:start-after: <!-- installation-end -->
```
