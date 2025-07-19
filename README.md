# Method Development Template Repository

This repository provides a standardized template for building, testing, documenting, and releasing Python-based methods. It includes recommended layouts, tooling integrations, and CI/CD workflows to streamline development.

---

## 📁 Repository Structure
```bash
├── src/                   # Source code
│   └── my_package/        # Your Python package modules
├── test/                  # Pytest test suite
├── doc/                   # Sphinx documentation source
├── .github/               # GitHub workflows and configuration
│   ├── workflows/         # CI/CD pipelines
│   │   ├── publish-package.yml
│   │   └── docker-publish.yml
│   └── dependabot.yml     # Dependabot configuration
├── env.yml                # Conda environment spec
├── requirements.txt       # pip dependencies
├── lint.sh                # Run Flake8
├── build_doc.sh           # Build Sphinx HTML docs
├── pyproject.toml         # Package metadata & build config
└── README.md              # This file

```


## ⚙️ Installation

**Conda** (preferred for RDKit compatibility):
```bash
conda env create -f env.yml
conda activate template-env
```
**Other**
```bash
python -m venv .venv
source .venv/bin/activate  # macOS/Linux
.\venv\Scripts\activate # Windows
pip install --upgrade pip
pip install -r requirements.txt
```

## 🛠 Development Workflow

1. **Source Layout**  
   Store all library code under `src/` (rename folder as desired).

2. **Unit Tests**  
   Keep tests in `tests/` and run via pytest:
   ```bash
   pytest --maxfail=1 --disable-warnings -q

3. **Linting & Formatting**  
   - Follow PEP 8 style and PEP 257 docstrings.  
   - Use `black` for automatic formatting and `flake8` to enforce quality.  
   - Run both via:
     ```bash
     ./lint.sh
     ```
   - Edit `lint.sh` to adjust rules or exclude files.

4. **Documentation**  
   - Write docs in `docs/` using Sphinx.  
   - Build locally with:
     ```bash
     ./build_doc.sh
     ```
   - Automate publishing via ReadTheDocs with `.readthedocs.yml`

5. **Dependency Management**  
   - Use `env.yml` for Conda (pins `rdkit>=2025.3.1`, `pytest`, `black`, `flake8`).  
   - Alternatively, manage with `requirements.txt` for pip users.

6. **Release & Packaging**  
   - Define package metadata in `pyproject.toml`.  
   - Local install:
     ```bash
     pip install .
     ```
   - CI/CD workflows:
     - **PyPI**: `.github/workflows/publish-package.yml`  
     - **Docker**: `.github/workflows/docker-publish.yml`

7. **Automated Dependency Updates**  
   - `.github/dependabot.yml` configured for weekly checks and PRs.

*Feel free to adapt any naming conventions or tooling versions to fit your project needs.*

## License

This project is licensed under MIT License - see the [License](LICENSE) file for details.

## Acknowledgments

This project has received funding from the European Unions Horizon Europe Doctoral Network programme under the Marie-Skłodowska-Curie grant agreement No 101072930 ([TACsy](https://tacsy.eu/) -- Training Alliance for Computational)