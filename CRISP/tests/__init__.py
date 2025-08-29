[tool.setuptools.packages.find]
namespaces = false
where = ["."]
include = ["CRISP*"]

[tool.setuptools.package-data]
CRISP = ["py.typed"]
"CRISP.tests" = ["*.py"]
