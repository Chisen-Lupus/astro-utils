# astro-utils

> **"Why do I need to write a plotting function every time I start a new project?"**

`astro-utils` is a Python package designed to streamline the process of visualizing and analyzing astronomical data. It includes various utility functions to help you explore and interpret your data more efficiently, eliminating repetitive setup tasks.

## Local Installation

This package is not published yet. To install `astro-utils` locally for development:

```sh
git clone https://github.com/Chisen-Lupus/astro-utils.git
cd astro-utils
pip install -e .
```

> **Note**: This package is under active development. Features and functionalities may be subject to change.

## Finished Utilities

### "Manual" DS9 Viewer

All inside a Jupyter Notebook block, the class `autils.jupyter.InteractivePlot` provides a simplified widget for displaying astronomical images and adjusting visualization parameters, without the need to open a new DS9 window.

```python
from autils import InteractivePlot
InteractivePlot(image_data, wcs=image_wcs)
```

![Manual DS9](fig/meome12.gif)

### Incantations

`incantations.ipynb` contains common packages and settings that would show up in most astronomy projects. I prepared this for copy-and-paste, to save time of typing these code everytime.

### Additional Utilities

_(More utilities and descriptions will be added as features are developed)_

## Contributing

If you would like to contribute to `astro-utils`, feel free to submit pull requests or open issues. Contributions of new utilities, documentation, or examples are always welcome!