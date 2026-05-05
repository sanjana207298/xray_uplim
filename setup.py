from setuptools import setup, find_packages

setup(
    name             = "xray_uplim",
    version          = "2.0.0",
    author           = "Sanjana Gupta",
    description      = "Unified X-ray non-detection upper limit calculator (NuSTAR, XMM-Newton, Swift, Chandra)",
    long_description = open("README.md").read(),
    long_description_content_type = "text/markdown",
    packages         = find_packages(),
    package_data     = {
        # Bundle the Swift PSF coefficient file shipped with the package
        "xray_uplim.swift": ["data/*.fits", "data/*.FTZ"],
    },
    python_requires  = ">=3.8",
    install_requires = [
        "numpy>=1.21",
        "scipy>=1.7",
        "astropy>=5.0",
        "matplotlib>=3.4",
        "openpyxl>=3.0",      # Excel output
        "pyyaml>=6.0",        # YAML config files (xray_uplim-cli)
    ],
    extras_require = {
        # pip install xray_uplim[gui]
        "gui":  ["PySide6>=6.4"],
        # pip install xray_uplim[full]  — GUI + everything
        "full": ["PySide6>=6.4"],
    },
    entry_points = {
        "console_scripts": [
            # GUI launcher (falls back to CLI if PySide6 is missing)
            "xray_uplim     = xray_uplim.__main__:main",
            # CLI-only launcher (YAML/JSON config file)
            "xray_uplim-cli = xray_uplim.cli:main",
        ],
    },
    classifiers = [
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
    ],
)
