# ETL process for TDC

> [!IMPORTANT]
> This code only works on Linux-based systems due to the use of the `pytdc` library. Specifically, `tiledbsoma` has no builds for Windows[^1].



## Requirements


## `pkg_resources` module error

It is possible that during the run of the code you may encounter the following error:

```python
import pkg_resources
ModuleNotFoundError: No module named 'pkg_resources'
```

In this case, please follow the following instructions inside the .venv environment, or call the python from your environment you want to use:


```bash
.venv/bin/python -m pip uninstall -y setuptools
.venv/bin/python -m pip install setuptools
```


[^1]: https://github.com/single-cell-data/TileDB-SOMA/issues/2631