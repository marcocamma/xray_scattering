"""Python setup.py for xray_scattering package"""
import io
import os
from setuptools import find_packages, setup


def read(*paths, **kwargs):
    """Read the contents of a text file safely.
    >>> read("xray_scattering", "VERSION")
    '0.1.0'
    >>> read("README.md")
    ...
    """

    content = ""
    with io.open(
        os.path.join(os.path.dirname(__file__), *paths),
        encoding=kwargs.get("encoding", "utf8"),
    ) as open_file:
        content = open_file.read().strip()
    return content


def read_requirements(path):
    return [
        line.strip()
        for line in read(path).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


setup(
    name="xray_scattering",
    version=read("xray_scattering", "VERSION"),
    description="xray_scattering created by marco cammarata",
    url="https://github.com/marcocamma/xray_scattering/",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    author="marcocamma",
    packages=find_packages(where="."),
    package_dir={"": "."},
    package_data={
        "xray_scattering.database": ["*.h5", "*.txt"],
        "xray_scattering.pdb": ["*.pdb",],
        "xray_scattering.xyz": ["*.xyz",],
        },
    include_package_data=True,
    install_requires=read_requirements("requirements.txt"),
)
