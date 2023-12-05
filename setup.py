from setuptools import setup, find_packages

if __name__ == "__main__":
    setup(version="1.0", packages=find_packages(where="./src", exclude=("./tests",)))
