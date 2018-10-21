default: test

build:
	python setup.py build
install:
	python setup.py install
test:
	python -m pytest -v --profile-svg --cov=two --cov-report="html" tests 
test:
	python -m pytest -v --cov=two --cov-report="html" tests 
profile:
	python -m pytest -v --profile-svg profiling
