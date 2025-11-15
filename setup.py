from setuptools import setup

with open("requirements.txt") as f:
    requirements = []
    for line in f:
        line = line.strip()
        if line and not line.startswith("#"):
            requirements.append(line)

setup(
    name="parkinsons-ec3-model",
    version="0.1.0",
    description="Parkinson's energetic collapse project with EC3 bistable model",
    author="Max Anfilofyev",
    python_requires=">=3.7",
    install_requires=requirements,
    packages=[],  # No packages to install, just standalone scripts
)
